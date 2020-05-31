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


#include <ovito/stdobj/gui/StdObjGui.h>
#include <ovito/stdobj/gui/properties/PropertyInspectionApplet.h>
#include <ovito/grid/objects/VoxelGrid.h>

namespace Ovito { namespace Grid {

/**
 * \brief Data inspector page for voxel grid objects.
 */
class VoxelGridInspectionApplet : public PropertyInspectionApplet
{
	Q_OBJECT
	OVITO_CLASS(VoxelGridInspectionApplet)
	Q_CLASSINFO("DisplayName", "Voxel Grids");

public:

	/// Constructor.
	Q_INVOKABLE VoxelGridInspectionApplet() : PropertyInspectionApplet(VoxelGrid::OOClass()) {}

	/// Returns the key value for this applet that is used for ordering the applet tabs.
	virtual int orderingKey() const override { return 210; }

	/// Lets the applet create the UI widget that is to be placed into the data inspector panel.
	virtual QWidget* createWidget(MainWindow* mainWindow) override;

protected:

	/// Determines the text shown in cells of the vertical header column.
	virtual QVariant headerColumnText(int section) override;

private Q_SLOTS:

	/// Is called when the user selects a different property container object in the list.
	void onCurrentContainerChanged(const DataObject* dataObject);

private:

	MainWindow* _mainWindow;
	QLabel* _gridInfoLabel;
};

}	// End of namespace
}	// End of namespace
