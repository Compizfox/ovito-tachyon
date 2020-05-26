////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2020 Alexander Stukowski
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

#pragma once


#include <ovito/gui/desktop/GUI.h>
#include <ovito/core/oo/OvitoObject.h>
#include <ovito/core/dataset/pipeline/PipelineFlowState.h>
#include <ovito/core/dataset/scene/PipelineSceneNode.h>

namespace Ovito {

/**
 * \brief Abstract base class for applets shown in the data inspector.
 */
class OVITO_GUI_EXPORT DataInspectionApplet : public OvitoObject
{
	Q_OBJECT
	OVITO_CLASS(DataInspectionApplet)

public:

	/// Returns the key value for this applet that is used for ordering the applet tabs.
	virtual int orderingKey() const { return std::numeric_limits<int>::max(); }

	/// Determines whether the given pipeline data contains data that can be displayed by this applet.
	virtual bool appliesTo(const DataCollection& data);

	/// Determines the list of data objects that are displayed by the applet.
	virtual std::vector<ConstDataObjectPath> getDataObjectPaths() {
		return currentState().getObjectsRecursive(_dataObjectClass); 
	}

	/// Lets the applet create the UI widget that is to be placed into the data inspector panel.
	virtual QWidget* createWidget(MainWindow* mainWindow) = 0;

	/// Creates and returns the list widget displaying the list of data object objects.
	QListWidget* objectSelectionWidget();

	/// Lets the applet update the contents displayed in the inspector.
	virtual void updateDisplay(const PipelineFlowState& state, PipelineSceneNode* pipeline);

	/// This is called when the applet is no longer visible.
	virtual void deactivate(MainWindow* mainWindow) {}

	/// Selects a specific data object in this applet.
	virtual bool selectDataObject(PipelineObject* dataSource, const QString& objectIdentifierHint, const QVariant& modeHint);

	/// Returns the currently selected data pipeline in the scene.
	PipelineSceneNode* currentPipeline() const { return _pipelineNode.data(); }

	/// Returns the current output of the data pipeline displayed in the applet.
	const PipelineFlowState& currentState() const { return _pipelineState; }

	/// Returns the data object that is currently selected.
	const DataObject* selectedDataObject() const { return _selectedDataObject; }

	/// Returns the data collection path of the currently selected data object.
	const ConstDataObjectPath& selectedDataObjectPath() const { return _selectedDataObjectPath; }

protected:

	/// Constructor.
	DataInspectionApplet(const DataObject::OOMetaClass& dataObjectClass) : _dataObjectClass(dataObjectClass) {}

	/// Updates the list of data objects displayed in the inspector.
	void updateDataObjectList();

Q_SIGNALS:

	/// This signal is emitted when the user selects a different data object in the list.
	void currentObjectChanged(const DataObject* dataObject);

public:

	/// A specialized QTableView widget, which allows copying the selected contents of the
	/// table to the clipboard.
	class OVITO_GUI_EXPORT TableView : public QTableView
	{
	public:

		/// Constructor.
		TableView(QWidget* parent = nullptr) : QTableView(parent) {
			setWordWrap(false);
		}

	protected:

		/// Handles key press events for this widget.
		virtual void keyPressEvent(QKeyEvent* event) override;
	};

private:

	/// The type of data object displayed by this applet.
	const DataObject::OOMetaClass& _dataObjectClass;

	/// The widget for selecting the current data object.
	QListWidget* _objectSelectionWidget = nullptr;

	/// The path of the currently selected data object.
	ConstDataObjectPath _selectedDataObjectPath;

	/// The identifier path of the currently selected data object.
	QString _selectedDataObjectPathString;

	/// Pointer to the currently selected data object.
	const DataObject* _selectedDataObject = nullptr;

	/// The currently selected pipeline in the scene.
	QPointer<PipelineSceneNode> _pipelineNode;

	/// The pipeline output data being displayed.
	PipelineFlowState _pipelineState;
};

}	// End of namespace
