///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2019) Alexander Stukowski
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

#include <ovito/gui/GUI.h>
#include <ovito/gui/properties/FloatParameterUI.h>
#include <ovito/gui/properties/BooleanParameterUI.h>
#include <ovito/grid/objects/VoxelGridVis.h>
#include "VoxelGridVisEditor.h"

namespace Ovito { namespace Grid { OVITO_BEGIN_INLINE_NAMESPACE(Internal)

IMPLEMENT_OVITO_CLASS(VoxelGridVisEditor);
SET_OVITO_OBJECT_EDITOR(VoxelGridVis, VoxelGridVisEditor);

/******************************************************************************
* Sets up the UI widgets of the editor.
******************************************************************************/
void VoxelGridVisEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	// Create a rollout.
	QWidget* rollout = createRollout(tr("Grid display"), rolloutParams);

    // Create the rollout contents.
	QGridLayout* layout = new QGridLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(4);
	layout->setColumnStretch(1, 1);

	FloatParameterUI* transparencyUI = new FloatParameterUI(this, PROPERTY_FIELD(VoxelGridVis::transparencyController));
	layout->addWidget(transparencyUI->label(), 1, 0);
	layout->addLayout(transparencyUI->createFieldLayout(), 1, 1);

	BooleanParameterUI* highlightLinesUI = new BooleanParameterUI(this, PROPERTY_FIELD(VoxelGridVis::highlightGridLines));
	layout->addWidget(highlightLinesUI->checkBox(), 2, 0, 1, 2);
}

OVITO_END_INLINE_NAMESPACE
}	// End of namespace
}	// End of namespace