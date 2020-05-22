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

namespace Ovito {

/******************************************************************************
* A push button that allows the user to open a modifier's output in the data 
* inspector panel.
******************************************************************************/
class OVITO_GUI_EXPORT OpenDataInspectorButton : public QPushButton
{
	Q_OBJECT

public:

	/// Constructor.
	OpenDataInspectorButton(ModifierPropertiesEditor* editor, const QString& buttonTitle, const QString& objectNameHint = {}, const QVariant& modeHint = {});

	/// Returns the properties editor hosting this button.
	ModifierPropertiesEditor* editor() const { return _editor; }

private:

	/// The properties editor hosting this button.
	ModifierPropertiesEditor* _editor;

	/// Data object name hint to be passed to the data inspector when the button is clicked.
	QString _objectNameHint;
	
	/// Mode hint to be passed to the data inspector when the button is clicked.
	QVariant _modeHint;
};

}	// End of namespace
