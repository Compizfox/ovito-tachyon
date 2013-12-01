///////////////////////////////////////////////////////////////////////////////
// 
//  Copyright (2013) Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  OVITO is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include <core/Core.h>
#include <core/scene/SelectionSet.h>

namespace Ovito {

IMPLEMENT_SERIALIZABLE_OVITO_OBJECT(Core, SelectionSet, RefTarget)
DEFINE_FLAGS_VECTOR_REFERENCE_FIELD(SelectionSet, _selection, "SelectedNodes", SceneNode, PROPERTY_FIELD_NEVER_CLONE_TARGET)
SET_PROPERTY_FIELD_LABEL(SelectionSet, _selection, "Nodes")

/******************************************************************************
* Default constructor.
******************************************************************************/
SelectionSet::SelectionSet(DataSet* dataset) : RefTarget(dataset)
{
	INIT_PROPERTY_FIELD(SelectionSet::_selection);
}

/******************************************************************************
* Adds a scene node to this selection set. 
******************************************************************************/
void SelectionSet::add(SceneNode* node)
{
	OVITO_CHECK_OBJECT_POINTER(node);
	if(contains(node)) return;
			
	// Insert into children array.
	_selection.push_back(node);
	OVITO_ASSERT(contains(node));
}

/******************************************************************************
* Adds multiple scene nodes to this selection set. 
******************************************************************************/
void SelectionSet::addAll(const QVector<SceneNode*>& nodes)
{
	for(SceneNode* node : nodes)
		add(node);
}

/******************************************************************************
* Completely replaces the contents of the selection set.  
******************************************************************************/
void SelectionSet::setNodes(const QVector<SceneNode*>& nodes)
{
	// Remove all nodes from the selection set that are not in the new list of selected nodes.
	for(int i = _selection.size(); i--; ) {
		if(!nodes.contains(_selection[i]))
			_selection.remove(i);
	}
	addAll(nodes);
}

/******************************************************************************
* Resets the selection set to contain only a single node.  
******************************************************************************/
void SelectionSet::setNode(SceneNode* node)
{
	OVITO_CHECK_POINTER(node);
	if(!_selection.contains(node)) {
		clear();
		add(node);
	}
	else {
		// Remove all other nodes from the selection set.
		for(int i = _selection.size(); i--; ) {
			if(node != _selection[i])
				_selection.remove(i);
		}
	}
}

/******************************************************************************
* Removes a scene node from this selection set. 
******************************************************************************/
void SelectionSet::remove(SceneNode* node)
{
	int index = _selection.indexOf(node);
	if(index == -1) return;	
	_selection.remove(index);
	OVITO_ASSERT(!contains(node));
}

/******************************************************************************
* Clears the selection.
******************************************************************************/
void SelectionSet::clear()
{
	_selection.clear();
}

/******************************************************************************
* From RefMaker.
******************************************************************************/
bool SelectionSet::referenceEvent(RefTarget* source, ReferenceEvent* event)
{
	if(event->type() == ReferenceEvent::TargetChanged) {
        // Convert changes messages from scene nodes.
        SceneNode* sourceNode = dynamic_object_cast<SceneNode>(source);
        OVITO_CHECK_OBJECT_POINTER(sourceNode);
        if(sourceNode) {
        	NodeInSelectionSetChangedEvent e(this, sourceNode, event);
        	notifyDependents(e);
        }
	}
	// Do not propagate events from selected nodes.
	return false;
}

/******************************************************************************
* Returns the bounding box that includes all selected nodes.
******************************************************************************/
Box3 SelectionSet::boundingBox(TimePoint time)
{
	Box3 bb;
	for(SceneNode* node : nodes()) {
		// Get node's world bounding box
		// and add it to global box.
		bb.addBox(node->worldBoundingBox(time));
	}
	return bb;
}

};
