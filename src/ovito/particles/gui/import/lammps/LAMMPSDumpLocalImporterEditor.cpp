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

#include <ovito/particles/gui/ParticlesGui.h>
#include <ovito/particles/import/lammps/LAMMPSDumpLocalImporter.h>
#include <ovito/stdobj/gui/properties/InputColumnMappingDialog.h>
#include <ovito/gui/desktop/properties/BooleanParameterUI.h>
#include <ovito/gui/desktop/mainwin/MainWindow.h>
#include <ovito/core/dataset/DataSetContainer.h>
#include <ovito/core/dataset/io/FileSource.h>
#include <ovito/core/utilities/concurrent/TaskManager.h>
#include "LAMMPSDumpLocalImporterEditor.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(LAMMPSDumpLocalImporterEditor);
SET_OVITO_OBJECT_EDITOR(LAMMPSDumpLocalImporter, LAMMPSDumpLocalImporterEditor);

/******************************************************************************
* This is called by the system when the user has selected a new file to import.
******************************************************************************/
bool LAMMPSDumpLocalImporterEditor::inspectNewFile(FileImporter* importer, const QUrl& sourceFile, QWidget* parent)
{
	// Retrieve column information of input file.
	LAMMPSDumpLocalImporter* lammpsImporter = static_object_cast<LAMMPSDumpLocalImporter>(importer);
	Future<InputColumnMapping> inspectFuture = lammpsImporter->inspectFileHeader(FileSourceImporter::Frame(sourceFile));
	if(!importer->dataset()->taskManager().waitForFuture(inspectFuture))
		return false;
	InputColumnMapping mapping = inspectFuture.result();

	// If this is a newly created file importer, load old mapping from application settings store.
	if(lammpsImporter->columnMapping().empty()) {
		QSettings settings;
		settings.beginGroup("viz/importer/lammps_dump_local/");
		if(settings.contains("colmapping")) {
			try {
				InputColumnMapping storedMapping;
				storedMapping.fromByteArray(settings.value("colmapping").toByteArray(), importer->dataset()->taskManager());
				std::copy_n(storedMapping.begin(), std::min(storedMapping.size(), mapping.size()), mapping.begin());
			}
			catch(Exception& ex) {
				ex.prependGeneralMessage(tr("Failed to load last used column-to-property mapping from application settings store."));
				ex.logError();
			}
		}
	}

	InputColumnMappingDialog dialog(mapping, parent, importer->dataset()->taskManager());
	if(dialog.exec() == QDialog::Accepted) {
		lammpsImporter->setColumnMapping(dialog.mapping());
		return true;
	}

	return false;
}

/******************************************************************************
 * Displays a dialog box that allows the user to edit the file column to property mapping.
 *****************************************************************************/
bool LAMMPSDumpLocalImporterEditor::showEditColumnMappingDialog(LAMMPSDumpLocalImporter* importer, const FileSourceImporter::Frame& frame, MainWindow* mainWindow)
{
	Future<InputColumnMapping> inspectFuture = importer->inspectFileHeader(frame);
	if(!importer->dataset()->taskManager().waitForFuture(inspectFuture))
		return false;
	InputColumnMapping mapping = inspectFuture.result();

	if(!importer->columnMapping().empty()) {
		InputColumnMapping newMapping = importer->columnMapping();
		newMapping.resize(mapping.size());
		for(size_t i = 0; i < newMapping.size(); i++)
			newMapping[i].columnName = mapping[i].columnName;
		mapping = std::move(newMapping);
	}

	InputColumnMappingDialog dialog(mapping, mainWindow, importer->dataset()->taskManager());
	if(dialog.exec() == QDialog::Accepted) {
		importer->setColumnMapping(dialog.mapping());
		// Remember the user-defined mapping for the next time.
		QSettings settings;
		settings.beginGroup("viz/importer/lammps_dump_local/");
		settings.setValue("colmapping", dialog.mapping().toByteArray(importer->dataset()->taskManager()));
		settings.endGroup();
		return true;
	}
	return false;
}

/******************************************************************************
* Sets up the UI widgets of the editor.
******************************************************************************/
void LAMMPSDumpLocalImporterEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	// Create a rollout.
	QWidget* rollout = createRollout(tr("LAMMPS dump local reader"), rolloutParams);

    // Create the rollout contents.
	QVBoxLayout* layout = new QVBoxLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(4);

	QGroupBox* optionsBox = new QGroupBox(tr("Options"), rollout);
	QVBoxLayout* sublayout = new QVBoxLayout(optionsBox);
	sublayout->setContentsMargins(4,4,4,4);
	layout->addWidget(optionsBox);

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
	sublayout->addWidget(multitimestepUI->checkBox());

	QGroupBox* columnMappingBox = new QGroupBox(tr("File columns"), rollout);
	sublayout = new QVBoxLayout(columnMappingBox);
	sublayout->setContentsMargins(4,4,4,4);
	layout->addWidget(columnMappingBox);

	QPushButton* editMappingButton = new QPushButton(tr("Edit column mapping..."));
	sublayout->addWidget(editMappingButton);
	connect(editMappingButton, &QPushButton::clicked, this, &LAMMPSDumpLocalImporterEditor::onEditColumnMapping);
}

/******************************************************************************
* Is called when the user pressed the "Edit column mapping" button.
******************************************************************************/
void LAMMPSDumpLocalImporterEditor::onEditColumnMapping()
{
	if(LAMMPSDumpLocalImporter* importer = static_object_cast<LAMMPSDumpLocalImporter>(editObject())) {
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
