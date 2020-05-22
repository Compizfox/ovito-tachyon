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
#include <ovito/gui/desktop/properties/ModifierPropertiesEditor.h>
#include <ovito/gui/desktop/mainwin/MainWindow.h>
#include "OpenDataInspectorButton.h"

namespace Ovito {

/******************************************************************************
* Constructor.
******************************************************************************/
OpenDataInspectorButton::OpenDataInspectorButton(ModifierPropertiesEditor* editor, const QString& buttonTitle, const QString& objectNameHint, const QVariant& modeHint)
	: QPushButton(buttonTitle), _editor(editor), _objectNameHint(objectNameHint), _modeHint(modeHint)
{
	connect(this, &QPushButton::clicked, [this]() {
		if(!this->editor()->modifierApplication() || !this->editor()->modifierApplication()->modifier() || !this->editor()->modifierApplication()->modifier()->isEnabled()) {
			QToolTip::showText(mapToGlobal(QPoint(0, height()/2)), tr("No results available, because modifier is turned off."), 
				this->editor()->container(), this->editor()->container()->rect(), 3000);
		}
		else if(!this->editor()->mainWindow()->openDataInspector(this->editor()->modifierApplication(), _objectNameHint, _modeHint)) {
			QToolTip::showText(mapToGlobal(QPoint(0, height()/2)), tr("Results not available yet. Try again later."), 
				this->editor()->container(), this->editor()->container()->rect(), 3000);
		}
		else {
			QToolTip::hideText();
		}
	});
}

}	// End of namespace
