////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2013 Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify it either under the
//  terms of the GNU General Public License version 3 as published by the Free Software
//  Foundation (the "GPL") or, at your option, under the terms of the MIT License.
//  If you do not alter this notice, a recipient may use your version of this
//  file under either the GPL or the MIT License.
//
//  You should have received a copy of the GPL along with this program in a
//  file LICENSE.GPL.txt.  You should have received a copy of the MIT License along
//  with this program in a file LICENSE.MIT.txt
//
//  This software is distributed on an "AS IS" basis, WITHOUT WARRANTY OF ANY KIND,
//  either express or implied. See the GPL or the MIT License for the specific language
//  governing rights and limitations.
//
////////////////////////////////////////////////////////////////////////////////////////

#include <ovito/gui/GUI.h>
#include <ovito/core/dataset/scene/SceneNode.h>
#include <ovito/core/dataset/scene/RootSceneNode.h>
#include <ovito/core/dataset/DataSetContainer.h>
#include <ovito/core/dataset/DataSet.h>
#include "SceneNodesListModel.h"

namespace Ovito { OVITO_BEGIN_INLINE_NAMESPACE(Gui) OVITO_BEGIN_INLINE_NAMESPACE(Internal)

/******************************************************************************
* Constructs the model.
******************************************************************************/
SceneNodesListModel::SceneNodesListModel(DataSetContainer& datasetContainer, QWidget* parent) : QAbstractListModel(parent),
		_datasetContainer(datasetContainer)
{
	// Listen for changes of the data set.
	connect(&datasetContainer, &DataSetContainer::dataSetChanged, this, &SceneNodesListModel::onDataSetChanged);

	// Listen for events of the root node.
	connect(&_rootNodeListener, &RefTargetListener<RootSceneNode>::notificationEvent, this, &SceneNodesListModel::onRootNodeNotificationEvent);

	// Listen for events of the other scene nodes.
	connect(&_nodeListener, &VectorRefTargetListener<SceneNode>::notificationEvent, this, &SceneNodesListModel::onNodeNotificationEvent);
}

/******************************************************************************
* Returns the number of rows of the model.
******************************************************************************/
int SceneNodesListModel::rowCount(const QModelIndex& parent) const
{
	return _nodeListener.targets().size();
}

/******************************************************************************
* Returns the model's data stored under the given role for the item referred to by the index.
******************************************************************************/
QVariant SceneNodesListModel::data(const QModelIndex& index, int role) const
{
	if(role == Qt::DisplayRole) {
		return QVariant::fromValue(_nodeListener.targets()[index.row()]->objectTitle());
	}
	else if(role == Qt::UserRole) {
		return QVariant::fromValue(_nodeListener.targets()[index.row()]);
	}
	return QVariant();
}

/******************************************************************************
* This is called when a new dataset has been loaded.
******************************************************************************/
void SceneNodesListModel::onDataSetChanged(DataSet* newDataSet)
{
	beginResetModel();
	_nodeListener.clear();
	_rootNodeListener.setTarget(nullptr);
	if(newDataSet) {
		_rootNodeListener.setTarget(newDataSet->sceneRoot());
		newDataSet->sceneRoot()->visitChildren([this](SceneNode* node) -> bool {
			_nodeListener.push_back(node);
			return true;
		});
	}
	endResetModel();
}

/******************************************************************************
* This handles reference events generated by the root node.
******************************************************************************/
void SceneNodesListModel::onRootNodeNotificationEvent(const ReferenceEvent& event)
{
	onNodeNotificationEvent(_rootNodeListener.target(), event);
}

/******************************************************************************
* This handles reference events generated by the scene nodes.
******************************************************************************/
void SceneNodesListModel::onNodeNotificationEvent(RefTarget* source, const ReferenceEvent& event)
{
	// Whenever a new node is being inserted into the scene, add it to our internal list.
	if(event.type() == ReferenceEvent::ReferenceAdded) {
		const ReferenceFieldEvent& refEvent = static_cast<const ReferenceFieldEvent&>(event);
		if(refEvent.field() == &PROPERTY_FIELD(SceneNode::children)) {
			if(SceneNode* node = dynamic_object_cast<SceneNode>(refEvent.newTarget())) {
				beginInsertRows(QModelIndex(), _nodeListener.targets().size(), _nodeListener.targets().size());
				_nodeListener.push_back(node);
				endInsertRows();
				// Add all child nodes too.
				node->visitChildren([this](SceneNode* node) -> bool {
					beginInsertRows(QModelIndex(), _nodeListener.targets().size(), _nodeListener.targets().size());
					_nodeListener.push_back(node);
					endInsertRows();
					return true;
				});
			}
		}
	}

	// If a node is being removed from the scene, remove it from our internal list.
	if(event.type() == ReferenceEvent::ReferenceRemoved) {
		// Don't know how else to do this in a safe manner. Rebuild the entire model from scratch.
		onDataSetChanged(_datasetContainer.currentSet());
	}

	// If a node is being renamed, let the model emit an update signal.
	if(event.type() == ReferenceEvent::TitleChanged) {
		int index = _nodeListener.targets().indexOf(static_cast<SceneNode*>(source));
		if(index >= 0) {
			QModelIndex modelIndex = createIndex(index, 0, source);
			dataChanged(modelIndex, modelIndex);
		}
	}
}

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace
