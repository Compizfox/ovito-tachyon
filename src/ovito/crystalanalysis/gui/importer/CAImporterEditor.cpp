////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2017 Alexander Stukowski
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

#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/crystalanalysis/importer/CAImporter.h>
#include <ovito/gui/desktop/properties/BooleanParameterUI.h>
#include "CAImporterEditor.h"

namespace Ovito { namespace CrystalAnalysis {

IMPLEMENT_OVITO_CLASS(CAImporterEditor);
SET_OVITO_OBJECT_EDITOR(CAImporter, CAImporterEditor);

/******************************************************************************
* Sets up the UI widgets of the editor.
******************************************************************************/
void CAImporterEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	// Create a rollout.
	QWidget* rollout = createRollout(tr("Crystal analysis file"), rolloutParams);

    // Create the rollout contents.
	QVBoxLayout* layout = new QVBoxLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(4);

	// Multi-timestep file
	BooleanParameterUI* multitimestepUI = new BooleanParameterUI(this, PROPERTY_FIELD(FileSourceImporter::isMultiTimestepFile));
	// The following signal handler updates the parameter UI whenever the isMultiTimestepFile parameter of the current file source importer changes.
	// It is needed, because target-changed messages are surpressed for this property field and the normal update mechanism for the parameter UI doesn't work.
	connect(this, &PropertiesEditor::contentsReplaced, this, [con = QMetaObject::Connection(), multitimestepUI = multitimestepUI](RefTarget* editObject) mutable {
#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4573)
#endif
		disconnect(con);
		con = editObject ? connect(static_object_cast<FileSourceImporter>(editObject), &FileSourceImporter::isMultiTimestepFileChanged, multitimestepUI, &ParameterUI::updateUI) : QMetaObject::Connection();
#ifdef _MSC_VER
	#pragma warning(pop)
#endif
	});
	layout->addWidget(multitimestepUI->checkBox());
}

}	// End of namespace
}	// End of namespace
