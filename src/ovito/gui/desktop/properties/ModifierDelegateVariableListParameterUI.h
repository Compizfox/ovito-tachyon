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
#include "ParameterUI.h"

namespace Ovito {

/******************************************************************************
* A list view that shows the fixed set of delegates of a MultiDelegatingModifier,
* which can each be enabled or disabled by the user.
******************************************************************************/
class OVITO_GUI_EXPORT ModifierDelegateVariableListParameterUI : public ParameterUI
{
	Q_OBJECT
	OVITO_CLASS(ModifierDelegateVariableListParameterUI)

public:

	/// Constructor.
	ModifierDelegateVariableListParameterUI(QObject* parentEditor, const OvitoClass& delegateType);

	/// Destructor.
	virtual ~ModifierDelegateVariableListParameterUI();

	/// This method is called when a new editable object has been assigned to the properties owner this
	/// parameter UI belongs to. The parameter UI should react to this change appropriately and
	/// show the properties value for the new edit object in the UI.
	virtual void resetUI() override;

	/// This method updates the displayed value of the parameter UI.
	virtual void updateUI() override;

	/// This returns the container widget managed by this class.
	QWidget* containerWidget() const { return _containerWidget; }

	/// Sets the enabled state of the UI.
	virtual void setEnabled(bool enabled) override;

private Q_SLOTS:

	/// Is called when the user requested the addition of a new delegate to the modifier.
	void onAddDelegate();

	/// Is called when the user requested the removal of an existing delegate from the modifier.
	void onRemoveDelegate();

	/// Is called when the user selects an entry in a combo box.
	void onDelegateSelected(int index);

protected:

	/// This method is called when a reference target changes.
	virtual bool referenceEvent(RefTarget* source, const ReferenceEvent& event) override;

	/// Is called when a RefTarget has been added to a VectorReferenceField of this RefMaker.
	virtual void referenceInserted(const PropertyFieldDescriptor& field, RefTarget* newTarget, int listIndex) override;

	/// Is called when a RefTarget has been removed from a VectorReferenceField of this RefMaker.
	virtual void referenceRemoved(const PropertyFieldDescriptor& field, RefTarget* oldTarget, int listIndex) override;

private:

	/// The type of modifier delegates, which the user can choose from.
	const OvitoClass& _delegateType;

	/// The container widget managed by this parameter UI.
	QPointer<QWidget> _containerWidget;

	/// The QAction for each delegate that removes it. 
	QVector<QAction*> _removeDelegateActions;

	/// The QComboBox for each delegate. 
	QVector<QComboBox*> _delegateBoxes;

	/// The current list of delegates.
	DECLARE_VECTOR_REFERENCE_FIELD_FLAGS(ModifierDelegate, delegates, PROPERTY_FIELD_NO_UNDO | PROPERTY_FIELD_WEAK_REF | PROPERTY_FIELD_NO_CHANGE_MESSAGE);
};

}	// End of namespace
