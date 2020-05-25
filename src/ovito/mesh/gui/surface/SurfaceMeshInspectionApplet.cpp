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

#include <ovito/gui/desktop/GUI.h>
#include <ovito/gui/desktop/widgets/general/AutocompleteLineEdit.h>
#include "SurfaceMeshInspectionApplet.h"

namespace Ovito { namespace Mesh {

IMPLEMENT_OVITO_CLASS(SurfaceMeshInspectionApplet);

/******************************************************************************
* Lets the applet create the UI widget that is to be placed into the data
* inspector panel.
******************************************************************************/
QWidget* SurfaceMeshInspectionApplet::createWidget(MainWindow* mainWindow)
{
	QSplitter* splitter = new QSplitter();
	splitter->addWidget(objectSelectionWidget());

	QWidget* rightContainer = new QWidget();
	splitter->addWidget(rightContainer);
	splitter->setStretchFactor(0, 1);
	splitter->setStretchFactor(1, 3);

	QHBoxLayout* rightLayout = new QHBoxLayout(rightContainer);
	rightLayout->setContentsMargins(0,0,0,0);
	rightLayout->setSpacing(0);

	QToolBar* toolbar = new QToolBar();
	toolbar->setOrientation(Qt::Vertical);
	toolbar->setToolButtonStyle(Qt::ToolButtonIconOnly);
	toolbar->setIconSize(QSize(22,22));
	toolbar->setStyleSheet("QToolBar { padding: 0px; margin: 0px; border: 0px none black; spacing: 0px; }");

	QActionGroup* subobjectActionGroup = new QActionGroup(this);
	_switchToVerticesAction = subobjectActionGroup->addAction(QIcon(":/mesh/icons/vertices_view.svg"), tr("Vertices"));
	_switchToFacesAction = subobjectActionGroup->addAction(QIcon(":/mesh/icons/faces_view.svg"), tr("Faces"));
	_switchToRegionsAction = subobjectActionGroup->addAction(QIcon(":/mesh/icons/regions_view.svg"), tr("Regions"));
	toolbar->addAction(_switchToVerticesAction);
	toolbar->addAction(_switchToFacesAction);
	toolbar->addAction(_switchToRegionsAction);
	_switchToVerticesAction->setCheckable(true);
	_switchToFacesAction->setCheckable(true);
	_switchToRegionsAction->setCheckable(true);
	_switchToVerticesAction->setChecked(true);

	_stackedWidget = new QStackedWidget();
	rightLayout->addWidget(_stackedWidget, 1);
	rightLayout->addSpacing(6);
	rightLayout->addWidget(toolbar, 0);

	_verticesApplet = new SurfaceMeshVertexInspectionApplet(this);
	_verticesApplet->setParent(this);
	_stackedWidget->addWidget(_verticesApplet->createWidget(mainWindow));

	_facesApplet = new SurfaceMeshFaceInspectionApplet(this);
	_facesApplet->setParent(this);
	_stackedWidget->addWidget(_facesApplet->createWidget(mainWindow));

	_regionsApplet = new SurfaceMeshRegionInspectionApplet(this);
	_regionsApplet->setParent(this);
	_stackedWidget->addWidget(_regionsApplet->createWidget(mainWindow));

	connect(_switchToVerticesAction, &QAction::triggered, this, [this]() {
		_stackedWidget->setCurrentIndex(0);
	});
	connect(_switchToFacesAction, &QAction::triggered, this, [this]() {
		_stackedWidget->setCurrentIndex(1);
	});
	connect(_switchToRegionsAction, &QAction::triggered, this, [this]() {
		_stackedWidget->setCurrentIndex(2);
	});

	connect(this, &DataInspectionApplet::currentObjectChanged, this, &SurfaceMeshInspectionApplet::onCurrentDataObjectChanged);

	return splitter;
}

/******************************************************************************
* Is called when the user selects a different container object from the list.
******************************************************************************/
void SurfaceMeshInspectionApplet::onCurrentDataObjectChanged()
{
	_verticesApplet->updateDisplay(currentState(), currentPipeline());
	_facesApplet->updateDisplay(currentState(), currentPipeline());
	_regionsApplet->updateDisplay(currentState(), currentPipeline());
}

/******************************************************************************
* Selects a specific data object in this applet.
******************************************************************************/
bool SurfaceMeshInspectionApplet::selectDataObject(PipelineObject* dataSource, const QString& objectIdentifierHint, const QVariant& modeHint)
{
	// Let the base class switch to the right data object. 
	bool result = DataInspectionApplet::selectDataObject(dataSource, objectIdentifierHint, modeHint);
	
	if(result) {
		// The mode hint is used to switch between vertex/face/region views.
		bool ok;
		int mode = modeHint.toInt(&ok);
		if(ok) {
			if(mode == 0)
				_switchToVerticesAction->trigger();	// Vertex list view
			else if(mode == 1)
				_switchToFacesAction->trigger(); 	// Face list view
			else if(mode == 2)
				_switchToRegionsAction->trigger();	// Region list view
		}
	}

	return result;
}

/******************************************************************************
* Determines the list of data objects that are displayed by the applet.
******************************************************************************/
std::vector<ConstDataObjectPath> SurfaceMeshVertexInspectionApplet::getDataObjectPaths()
{
	ConstDataObjectPath path = _parentApplet->selectedDataObjectPath();
	if(path.empty()) return {};
	path.push_back(static_object_cast<SurfaceMesh>(path.back())->vertices());
	return { std::move(path) };
}

/******************************************************************************
* Determines the list of data objects that are displayed by the applet.
******************************************************************************/
std::vector<ConstDataObjectPath> SurfaceMeshFaceInspectionApplet::getDataObjectPaths()
{
	ConstDataObjectPath path = _parentApplet->selectedDataObjectPath();
	if(path.empty()) return {};
	path.push_back(static_object_cast<SurfaceMesh>(path.back())->faces());
	return { std::move(path) };
}

/******************************************************************************
* Determines the list of data objects that are displayed by the applet.
******************************************************************************/
std::vector<ConstDataObjectPath> SurfaceMeshRegionInspectionApplet::getDataObjectPaths()
{
	ConstDataObjectPath path = _parentApplet->selectedDataObjectPath();
	if(path.empty()) return {};
	path.push_back(static_object_cast<SurfaceMesh>(path.back())->regions());
	return { std::move(path) };
}

/******************************************************************************
* Lets the applet create the UI widget that is to be placed into the data
* inspector panel.
******************************************************************************/
QWidget* SurfaceMeshVertexInspectionApplet::createWidget(MainWindow* mainWindow)
{
	createBaseWidgets();

	QWidget* panel = new QWidget();
	QGridLayout* layout = new QGridLayout(panel);
	layout->setContentsMargins(0,0,0,0);
	layout->setSpacing(0);

	QToolBar* toolbar = new QToolBar();
	toolbar->setOrientation(Qt::Horizontal);
	toolbar->setToolButtonStyle(Qt::ToolButtonIconOnly);
	toolbar->setIconSize(QSize(18,18));
	toolbar->setStyleSheet("QToolBar { padding: 0px; margin: 0px; border: 0px none black; spacing: 0px; }");
	toolbar->addAction(resetFilterAction());
	layout->addWidget(toolbar, 0, 0);

	filterExpressionEdit()->setPlaceholderText(tr("Filter vertices list..."));
	layout->addWidget(filterExpressionEdit(), 0, 1);
	QSplitter* splitter = new QSplitter();
	splitter->setChildrenCollapsible(false);
	splitter->addWidget(tableView());
	layout->addWidget(splitter, 1, 0, 1, 2);
	layout->setRowStretch(1, 1);

	return panel;
}

/******************************************************************************
* Lets the applet create the UI widget that is to be placed into the data
* inspector panel.
******************************************************************************/
QWidget* SurfaceMeshFaceInspectionApplet::createWidget(MainWindow* mainWindow)
{
	createBaseWidgets();

	QWidget* panel = new QWidget();
	QGridLayout* layout = new QGridLayout(panel);
	layout->setContentsMargins(0,0,0,0);
	layout->setSpacing(0);

	QToolBar* toolbar = new QToolBar();
	toolbar->setOrientation(Qt::Horizontal);
	toolbar->setToolButtonStyle(Qt::ToolButtonIconOnly);
	toolbar->setIconSize(QSize(18,18));
	toolbar->setStyleSheet("QToolBar { padding: 0px; margin: 0px; border: 0px none black; spacing: 0px; }");
	toolbar->addAction(resetFilterAction());
	layout->addWidget(toolbar, 0, 0);

	filterExpressionEdit()->setPlaceholderText(tr("Filter faces list..."));
	layout->addWidget(filterExpressionEdit(), 0, 1);
	QSplitter* splitter = new QSplitter();
	splitter->setChildrenCollapsible(false);
	splitter->addWidget(tableView());
	layout->addWidget(splitter, 1, 0, 1, 2);
	layout->setRowStretch(1, 1);

	return panel;
}

/******************************************************************************
* Lets the applet create the UI widget that is to be placed into the data
* inspector panel.
******************************************************************************/
QWidget* SurfaceMeshRegionInspectionApplet::createWidget(MainWindow* mainWindow)
{
	createBaseWidgets();

	QWidget* panel = new QWidget();
	QGridLayout* layout = new QGridLayout(panel);
	layout->setContentsMargins(0,0,0,0);
	layout->setSpacing(0);

	QToolBar* toolbar = new QToolBar();
	toolbar->setOrientation(Qt::Horizontal);
	toolbar->setToolButtonStyle(Qt::ToolButtonIconOnly);
	toolbar->setIconSize(QSize(18,18));
	toolbar->setStyleSheet("QToolBar { padding: 0px; margin: 0px; border: 0px none black; spacing: 0px; }");
	toolbar->addAction(resetFilterAction());
	layout->addWidget(toolbar, 0, 0);

	filterExpressionEdit()->setPlaceholderText(tr("Filter regions list..."));
	layout->addWidget(filterExpressionEdit(), 0, 1);
	QSplitter* splitter = new QSplitter();
	splitter->setChildrenCollapsible(false);
	splitter->addWidget(tableView());
	layout->addWidget(splitter, 1, 0, 1, 2);
	layout->setRowStretch(1, 1);

	return panel;
}

}	// End of namespace
}	// End of namespace
