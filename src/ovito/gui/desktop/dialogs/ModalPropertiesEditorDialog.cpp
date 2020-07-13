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
#include <ovito/gui/desktop/properties/PropertiesPanel.h>
#include <ovito/gui/desktop/mainwin/MainWindow.h>
#include <ovito/core/dataset/UndoStack.h>
#include <ovito/core/dataset/DataSet.h>
#include "ModalPropertiesEditorDialog.h"

namespace Ovito {

/******************************************************************************
* The constructor of the dialog.
******************************************************************************/
ModalPropertiesEditorDialog::ModalPropertiesEditorDialog(RefTarget* object, OORef<PropertiesEditor> editor, QWidget* parent, MainWindow* mainWindow, const QString& dialogTitle, const QString& undoString, const QString& helpTopic) :
		QDialog(parent), _editor(std::move(editor)), UndoableTransaction(object->dataset()->undoStack(), undoString)
{
	setWindowTitle(dialogTitle);

	QVBoxLayout* layout = new QVBoxLayout(this);

	PropertiesPanel* propertiesPanel = new PropertiesPanel(this, mainWindow);
	propertiesPanel->setVisible(false);
	layout->addWidget(propertiesPanel, 1);
	_editor->initialize(propertiesPanel, mainWindow, RolloutInsertionParameters().insertInto(this), nullptr);
	_editor->setEditObject(object);

	QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel | QDialogButtonBox::Help, Qt::Horizontal, this);
	layout->addWidget(buttonBox);

	connect(buttonBox, &QDialogButtonBox::rejected, this, &ModalPropertiesEditorDialog::reject);
	connect(buttonBox, &QDialogButtonBox::accepted, this, [&]() {
		commit();
		accept();
	});
	connect(buttonBox, &QDialogButtonBox::helpRequested, [helpTopic]() {
		MainWindow::openHelpTopic(helpTopic);
	});
}

}	// End of namespace
