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
#include <ovito/mesh/surface/SurfaceMesh.h>
#include <ovito/stdobj/gui/properties/PropertyInspectionApplet.h>

namespace Ovito { namespace Mesh {

class SurfaceMeshInspectionApplet;  // defined below

/**
 * \brief Data inspector page for surface meshes vertices.
 */
class SurfaceMeshVertexInspectionApplet : public PropertyInspectionApplet
{
	Q_OBJECT

public:

	/// Constructor.
	SurfaceMeshVertexInspectionApplet(SurfaceMeshInspectionApplet* parentApplet) : PropertyInspectionApplet(SurfaceMeshVertices::OOClass()), _parentApplet(parentApplet) {}

	/// Lets the applet create the UI widget that is to be placed into the data inspector panel.
	virtual QWidget* createWidget(MainWindow* mainWindow) override;

	/// Determines the list of data objects that are displayed by the applet.
	virtual std::vector<ConstDataObjectPath> getDataObjectPaths() override;

private:

	/// The parent applet for the SurfaceMesh, which hosts this sub-object applet.
	SurfaceMeshInspectionApplet* _parentApplet;
};

/**
 * \brief Data inspector page for surface meshes faces.
 */
class SurfaceMeshFaceInspectionApplet : public PropertyInspectionApplet
{
	Q_OBJECT

public:

	/// Constructor.
	SurfaceMeshFaceInspectionApplet(SurfaceMeshInspectionApplet* parentApplet) : PropertyInspectionApplet(SurfaceMeshFaces::OOClass()), _parentApplet(parentApplet) {}

	/// Lets the applet create the UI widget that is to be placed into the data inspector panel.
	virtual QWidget* createWidget(MainWindow* mainWindow) override;

	/// Determines the list of data objects that are displayed by the applet.
	virtual std::vector<ConstDataObjectPath> getDataObjectPaths() override;

private:

	/// The parent applet for the SurfaceMesh, which hosts this sub-object applet.
	SurfaceMeshInspectionApplet* _parentApplet;
};

/**
 * \brief Data inspector page for surface meshes regions.
 */
class SurfaceMeshRegionInspectionApplet : public PropertyInspectionApplet
{
	Q_OBJECT

public:

	/// Constructor.
	SurfaceMeshRegionInspectionApplet(SurfaceMeshInspectionApplet* parentApplet) : PropertyInspectionApplet(SurfaceMeshRegions::OOClass()), _parentApplet(parentApplet) {}

	/// Lets the applet create the UI widget that is to be placed into the data inspector panel.
	virtual QWidget* createWidget(MainWindow* mainWindow) override;

	/// Determines the list of data objects that are displayed by the applet.
	virtual std::vector<ConstDataObjectPath> getDataObjectPaths() override;

private:

	/// The parent applet for the SurfaceMesh, which hosts this sub-object applet.
	SurfaceMeshInspectionApplet* _parentApplet;
};

/**
 * \brief Data inspector page for surface meshes.
 */
class SurfaceMeshInspectionApplet : public DataInspectionApplet
{
	Q_OBJECT
	OVITO_CLASS(SurfaceMeshInspectionApplet)
	Q_CLASSINFO("DisplayName", "Surfaces");

public:

	/// Constructor.
	Q_INVOKABLE SurfaceMeshInspectionApplet() : DataInspectionApplet(SurfaceMesh::OOClass()) {}

	/// Returns the key value for this applet that is used for ordering the applet tabs.
	virtual int orderingKey() const override { return 220; }

	/// Lets the applet create the UI widget that is to be placed into the data inspector panel.
	virtual QWidget* createWidget(MainWindow* mainWindow) override;

private Q_SLOTS:

	/// Is called when the user selects a different data object from the list.
	void onCurrentDataObjectChanged();

private:

	SurfaceMeshVertexInspectionApplet* _verticesApplet;
	SurfaceMeshFaceInspectionApplet* _facesApplet;
	SurfaceMeshRegionInspectionApplet* _regionsApplet;
	QStackedWidget* _stackedWidget;
	QAction* _switchToVerticesAction;
	QAction* _switchToFacesAction;
	QAction* _switchToRegionsAction;
};

}	// End of namespace
}	// End of namespace
