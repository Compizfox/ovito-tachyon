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

#include <ovito/stdobj/gui/StdObjGui.h>
#include <ovito/gui/desktop/mainwin/MainWindow.h>
#include "VoxelGridInspectionApplet.h"

namespace Ovito { namespace Grid {

IMPLEMENT_OVITO_CLASS(VoxelGridInspectionApplet);

/******************************************************************************
* Lets the applet create the UI widget that is to be placed into the data
* inspector panel.
******************************************************************************/
QWidget* VoxelGridInspectionApplet::createWidget(MainWindow* mainWindow)
{
	createBaseWidgets();
	_mainWindow = mainWindow;

	QSplitter* splitter = new QSplitter();
	splitter->addWidget(objectSelectionWidget());

	QWidget* rightContainer = new QWidget();
	splitter->addWidget(rightContainer);
	splitter->setStretchFactor(0, 1);
	splitter->setStretchFactor(1, 3);

	QHBoxLayout* rightLayout = new QHBoxLayout(rightContainer);
	rightLayout->setContentsMargins(0,0,0,0);
	rightLayout->setSpacing(4);

	_gridInfoLabel = new QLabel();
	_gridInfoLabel->setTextInteractionFlags(Qt::TextBrowserInteraction);
	_gridInfoLabel->setTextFormat(Qt::RichText);
	_gridInfoLabel->setMargin(3);
	_gridInfoLabel->setAlignment(Qt::AlignLeft | Qt::AlignTop);
	rightLayout->addWidget(tableView(), 1);
	rightLayout->addWidget(_gridInfoLabel);

	connect(this, &DataInspectionApplet::currentObjectChanged, this, &VoxelGridInspectionApplet::onCurrentContainerChanged);

	return splitter;
}

/******************************************************************************
* Is called when the user selects a different container object from the list.
******************************************************************************/
void VoxelGridInspectionApplet::onCurrentContainerChanged(const DataObject* dataObject)
{
	// Update the displayed information.
	if(const VoxelGrid* grid = static_object_cast<VoxelGrid>(dataObject)) {
		QString text = tr("<p><b>Grid cells:</b> ");
		if(grid->domain()->is2D() && grid->shape()[2] <= 1)
			text += tr("%1 x %2</p>").arg(grid->shape()[0]).arg(grid->shape()[1]);
		else
			text += tr("%1 x %2 x %3</p>").arg(grid->shape()[0]).arg(grid->shape()[1]).arg(grid->shape()[2]);
		const Vector3& v1 = grid->domain()->cellVector1();
		const Vector3& v2 = grid->domain()->cellVector2();
		const Vector3& v3 = grid->domain()->cellVector3();
		const Point3& origin = grid->domain()->cellOrigin();
		text += tr("<p><b>Grid vector 1:</b> (%1 %2 %3)</p>").arg(v1.x()).arg(v1.y()).arg(v1.z());
		text += tr("<p><b>Grid vector 2:</b> (%1 %2 %3)</p>").arg(v2.x()).arg(v2.y()).arg(v2.z());
		if(grid->domain()->is2D() && grid->shape()[2] <= 1)
			text += tr("<p><b>Grid vector 3:</b> -</p>");
		else
			text += tr("<p><b>Grid vector 3:</b> (%1 %2 %3)</p>").arg(v3.x()).arg(v3.y()).arg(v3.z());
		text += tr("<p><b>Grid origin:</b> (%1 %2 %3)</p>").arg(origin.x()).arg(origin.y()).arg(origin.z());
		_gridInfoLabel->setText(text);
	}
	else {
		_gridInfoLabel->setText({});
	}
}

/******************************************************************************
* Determines the text shown in cells of the vertical header column.
******************************************************************************/
QVariant VoxelGridInspectionApplet::headerColumnText(int section) 
{ 
	if(const VoxelGrid* grid = static_object_cast<VoxelGrid>(selectedContainerObject())) {
		std::array<size_t, 3> coords = grid->voxelCoords(section);
		if(grid->domain()->is2D() && grid->shape()[2] <= 1) {
			return QStringLiteral("(%1, %2)").arg(coords[0]).arg(coords[1]);
		}
		else {
			return QStringLiteral("(%1, %2, %3)").arg(coords[0]).arg(coords[1]).arg(coords[2]);
		}
	}
	return section;
}

}	// End of namespace
}	// End of namespace
