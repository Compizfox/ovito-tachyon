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

///////////////////////////////////////////////////////////////////////////////
//
//  This module implements import of AMBER-style NetCDF trajectory files.
//  For specification documents see <http://ambermd.org/netcdf/>.
//
//  Extensions to this specification are supported through OVITO's manual
//  column mappings.
//
//  A LAMMPS dump style for this file format can be found at
//  <https://github.com/pastewka/lammps-netcdf>.
//
//  An ASE trajectory container is found in ase.io.netcdftrajectory.
//  <https://wiki.fysik.dtu.dk/ase/epydoc/ase.io.netcdftrajectory-module.html>.
//
//  Please contact Lars Pastewka <lars.pastewka@iwm.fraunhofer.de> for
//  questions and suggestions.
//
///////////////////////////////////////////////////////////////////////////////

#include <ovito/particles/gui/ParticlesGui.h>
#include <ovito/netcdf/AMBERNetCDFImporter.h>
#include <ovito/stdobj/gui/properties/InputColumnMappingDialog.h>
#include <ovito/gui/desktop/properties/BooleanParameterUI.h>
#include <ovito/gui/desktop/properties/BooleanRadioButtonParameterUI.h>
#include <ovito/gui/desktop/mainwin/MainWindow.h>
#include <ovito/core/dataset/io/FileSource.h>
#include "AMBERNetCDFImporterEditor.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(AMBERNetCDFImporterEditor);
SET_OVITO_OBJECT_EDITOR(AMBERNetCDFImporter, AMBERNetCDFImporterEditor);

/******************************************************************************
 * Displays a dialog box that allows the user to edit the custom file column to particle
 * property mapping.
 *****************************************************************************/
bool AMBERNetCDFImporterEditor::showEditColumnMappingDialog(AMBERNetCDFImporter* importer, const FileSourceImporter::Frame& frame, MainWindow* mainWindow)
{
	Future<ParticleInputColumnMapping> inspectFuture = importer->inspectFileHeader(frame);
	if(!importer->dataset()->taskManager().waitForFuture(inspectFuture))
		return false;
	ParticleInputColumnMapping mapping = inspectFuture.result();

	if(!importer->customColumnMapping().empty()) {
		ParticleInputColumnMapping customMapping = importer->customColumnMapping();
		customMapping.resize(mapping.size());
		for(size_t i = 0; i < customMapping.size(); i++)
			customMapping[i].columnName = mapping[i].columnName;
		mapping = customMapping;
	}

	InputColumnMappingDialog dialog(mapping, mainWindow, importer->dataset()->taskManager());
	if(dialog.exec() == QDialog::Accepted) {
		importer->setCustomColumnMapping(dialog.mapping());
		importer->setUseCustomColumnMapping(true);
		return true;
	}
	else {
		return false;
	}
}

/******************************************************************************
* Sets up the UI widgets of the editor.
******************************************************************************/
void AMBERNetCDFImporterEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	// Create a rollout.
	QWidget* rollout = createRollout(tr("NetCDF file"), rolloutParams);

    // Create the rollout contents.
	QVBoxLayout* layout = new QVBoxLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(4);

	QGroupBox* optionsBox = new QGroupBox(tr("Options"), rollout);
	QVBoxLayout* sublayout = new QVBoxLayout(optionsBox);
	sublayout->setContentsMargins(4,4,4,4);
	layout->addWidget(optionsBox);

	// Sort particles
	BooleanParameterUI* sortParticlesUI = new BooleanParameterUI(this, PROPERTY_FIELD(ParticleImporter::sortParticles));
	sublayout->addWidget(sortParticlesUI->checkBox());

	QGroupBox* columnMappingBox = new QGroupBox(tr("File columns"), rollout);
	sublayout = new QVBoxLayout(columnMappingBox);
	sublayout->setContentsMargins(4,4,4,4);
	layout->addWidget(columnMappingBox);

	BooleanRadioButtonParameterUI* useCustomMappingUI = new BooleanRadioButtonParameterUI(this, PROPERTY_FIELD(AMBERNetCDFImporter::useCustomColumnMapping));
	useCustomMappingUI->buttonFalse()->setText(tr("Automatic mapping"));
	sublayout->addWidget(useCustomMappingUI->buttonFalse());
	useCustomMappingUI->buttonTrue()->setText(tr("User-defined mapping to particle properties"));
	sublayout->addWidget(useCustomMappingUI->buttonTrue());
	connect(useCustomMappingUI->buttonFalse(), &QRadioButton::clicked, this, [this]() {
		if(AMBERNetCDFImporter* importer = static_object_cast<AMBERNetCDFImporter>(editObject()))
			importer->requestReload();
	}, Qt::QueuedConnection);

	QPushButton* editMappingButton = new QPushButton(tr("Edit column mapping..."));
	sublayout->addWidget(editMappingButton);
	connect(editMappingButton, &QPushButton::clicked, this, &AMBERNetCDFImporterEditor::onEditColumnMapping);
}

/******************************************************************************
* Is called when the user pressed the "Edit column mapping" button.
******************************************************************************/
void AMBERNetCDFImporterEditor::onEditColumnMapping()
{
	if(AMBERNetCDFImporter* importer = static_object_cast<AMBERNetCDFImporter>(editObject())) {
		UndoableTransaction::handleExceptions(importer->dataset()->undoStack(), tr("Change file column mapping"), [this, importer]() {

			// Determine the currently loaded data file of the FileSource.
			FileSource* fileSource = importer->fileSource();
			if(!fileSource || fileSource->frames().empty()) return;
			int frameIndex = qBound(0, fileSource->dataCollectionFrame(), fileSource->frames().size()-1);

			// Show the dialog box, which lets the user modify the file column mapping.
			if(showEditColumnMappingDialog(importer, fileSource->frames()[frameIndex], mainWindow())) {
				importer->requestReload();
			}
		});
	}
}

}	// End of namespace
}	// End of namespace
