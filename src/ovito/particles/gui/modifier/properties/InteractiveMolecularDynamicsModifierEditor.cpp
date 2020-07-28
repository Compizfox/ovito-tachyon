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
#include <ovito/particles/modifier/properties/InteractiveMolecularDynamicsModifier.h>
#include <ovito/gui/desktop/properties/StringParameterUI.h>
#include <ovito/gui/desktop/properties/IntegerParameterUI.h>
#include "InteractiveMolecularDynamicsModifierEditor.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(InteractiveMolecularDynamicsModifierEditor);
SET_OVITO_OBJECT_EDITOR(InteractiveMolecularDynamicsModifier, InteractiveMolecularDynamicsModifierEditor);

/******************************************************************************
* Sets up the UI widgets of the editor.
******************************************************************************/
void InteractiveMolecularDynamicsModifierEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	QWidget* rollout = createRollout(tr("Interactive molecular dynamics (IMD)"), rolloutParams, "particles.modifiers.interactive_molecular_dynamics.html");

    // Create the rollout contents.
	QGridLayout* layout = new QGridLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(4);
	layout->setColumnStretch(1, 1);

	QGridLayout* sublayout = new QGridLayout();
	sublayout->setContentsMargins(0,0,0,0);
	sublayout->setHorizontalSpacing(4);
	sublayout->setVerticalSpacing(2);
	sublayout->setColumnStretch(0, 3);
	sublayout->setColumnStretch(1, 1);
	layout->addLayout(sublayout, 0, 0, 1, 2);

	// Server hostname.
	_serverHostnameUI = new StringParameterUI(this, PROPERTY_FIELD(InteractiveMolecularDynamicsModifier::hostName));
	sublayout->addWidget(new QLabel(tr("IMD Server:")), 0, 0);
	sublayout->addWidget(_serverHostnameUI->textBox(), 1, 0);

	// Server port.
	_serverPortUI = new IntegerParameterUI(this, PROPERTY_FIELD(InteractiveMolecularDynamicsModifier::port));
	sublayout->addWidget(_serverPortUI->label(), 0, 1);
	sublayout->addLayout(_serverPortUI->createFieldLayout(), 1, 1);

	// Connect button.
	_connectButton = new QPushButton();
	_connectButton->setEnabled(false);
	layout->addWidget(_connectButton, 1, 0, 1, 2);

	// Transmission interval.
	layout->setRowMinimumHeight(2, 10);
	IntegerParameterUI* transmissionIntervalUI = new IntegerParameterUI(this, PROPERTY_FIELD(InteractiveMolecularDynamicsModifier::transmissionInterval));
	layout->addWidget(transmissionIntervalUI->label(), 3, 0);
	layout->addLayout(transmissionIntervalUI->createFieldLayout(), 3, 1);

	// Status label.
	layout->setRowMinimumHeight(4, 10);
	layout->addWidget(statusLabel(), 5, 0, 1, 2);

	// Handle connect button.
	connect(_connectButton, &QPushButton::clicked, this, [&]() {
		InteractiveMolecularDynamicsModifier* modifier = static_object_cast<InteractiveMolecularDynamicsModifier>(editObject());
		if(!modifier) return;

		if(modifier->socket().state() == QAbstractSocket::UnconnectedState)
			modifier->connectToServer();
		else
			modifier->disconnectFromServer();
	});

	// Whenever a new modifier gets loaded into the editor:
	connect(this, &PropertiesEditor::contentsReplaced, this, [this, con = QMetaObject::Connection()](RefTarget* editObject) mutable {
		disconnect(con);

		// Update displayed information.
		connectionStateChanged();

		// Update UI when IMD server connection changes.
		con = editObject ? connect(&static_object_cast<InteractiveMolecularDynamicsModifier>(editObject)->socket(), &QAbstractSocket::stateChanged, this, &InteractiveMolecularDynamicsModifierEditor::connectionStateChanged) : QMetaObject::Connection();
	});

	// Handle parameter change.
	connect(this, &PropertiesEditor::contentsChanged, this, &InteractiveMolecularDynamicsModifierEditor::connectionStateChanged);
}

/******************************************************************************
*  Is called whenever the IMD connection state of the modifier changes.
******************************************************************************/
void InteractiveMolecularDynamicsModifierEditor::connectionStateChanged()
{
	InteractiveMolecularDynamicsModifier* modifier = static_object_cast<InteractiveMolecularDynamicsModifier>(editObject());

	_connectButton->setEnabled(modifier && (modifier->socket().state() != QAbstractSocket::UnconnectedState || !modifier->hostName().trimmed().isEmpty()));
	_connectButton->setText(!modifier || modifier->socket().state() == QAbstractSocket::UnconnectedState ?
		tr("Connect") : tr("Disconnect"));
	_serverHostnameUI->setEnabled(modifier && modifier->socket().state() == QAbstractSocket::UnconnectedState);
	_serverPortUI->setEnabled(modifier && modifier->socket().state() == QAbstractSocket::UnconnectedState);
}

}	// End of namespace
}	// End of namespace
