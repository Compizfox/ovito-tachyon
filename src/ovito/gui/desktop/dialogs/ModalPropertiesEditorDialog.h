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
#include <ovito/core/dataset/UndoStack.h>

namespace Ovito {

/**
 * A dialog box displaying a PropertiesEditor for an object.
 */
class ModalPropertiesEditorDialog : public QDialog, private UndoableTransaction
{
	Q_OBJECT

public:

	/// Constructor.
	ModalPropertiesEditorDialog(RefTarget* object, OORef<PropertiesEditor> editor, QWidget* parentWindow, MainWindow* mainWindow, const QString& dialogTitle, const QString& undoString, const QString& helpTopic);

private:

	OORef<PropertiesEditor> _editor;
};

}	// End of namespace
