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
#include <ovito/core/dataset/pipeline/DelegatingModifier.h>
#include <ovito/gui/desktop/properties/ModifierPropertiesEditor.h>
#include "ModifierDelegateVariableListParameterUI.h"
#include "ModifierDelegateParameterUI.h"

namespace Ovito {

IMPLEMENT_OVITO_CLASS(ModifierDelegateVariableListParameterUI);
DEFINE_REFERENCE_FIELD(ModifierDelegateVariableListParameterUI, delegates);

/******************************************************************************
* Constructor.
******************************************************************************/
ModifierDelegateVariableListParameterUI::ModifierDelegateVariableListParameterUI(QObject* parentEditor, const OvitoClass& delegateType) :
	ParameterUI(parentEditor), 
	_delegateType(delegateType),
	_containerWidget(new QWidget())
{
	QVBoxLayout* layout = new QVBoxLayout(_containerWidget);
	layout->setContentsMargins(0,0,0,0);

	QToolBar* toolbar = new QToolBar();
	toolbar->setStyleSheet("QToolBar { padding: 0px; margin: 0px; border: 0px none black; }");
	toolbar->setFloatable(false);
	toolbar->setIconSize(QSize(16,16));
	QAction* addDelegateAction = toolbar->addAction(QIcon(":/gui/actions/animation/add_animation_key.bw.svg"), tr("Add entry"));
	connect(addDelegateAction, &QAction::triggered, this, &ModifierDelegateVariableListParameterUI::onAddDelegate);
	layout->addWidget(toolbar, 0, Qt::AlignRight | Qt::AlignTop);
}

/******************************************************************************
* Destructor.
******************************************************************************/
ModifierDelegateVariableListParameterUI::~ModifierDelegateVariableListParameterUI()
{
	clearAllReferences();

	// Release widget.
	delete containerWidget();
}

/******************************************************************************
* This method is called when a new editable object has been assigned to the properties owner this
* parameter UI belongs to. The parameter UI should react to this change appropriately and
* show the properties value for the new edit object in the UI.
******************************************************************************/
void ModifierDelegateVariableListParameterUI::resetUI()
{
	// Create our own copy of the list of delegates of the modifier.
	if(editObject())
		_delegates.set(this, PROPERTY_FIELD(delegates), static_object_cast<MultiDelegatingModifier>(editObject())->delegates());
	else
		_delegates.clear(this, PROPERTY_FIELD(delegates));

	// Update enabled state of UI widget.
	if(containerWidget())
		containerWidget()->setEnabled(editObject() && isEnabled());

	ParameterUI::resetUI();
}

/******************************************************************************
* This method is called when a new editable object has been assigned to the
* properties owner this parameter UI belongs to.
******************************************************************************/
void ModifierDelegateVariableListParameterUI::updateUI()
{
	ParameterUI::updateUI();

	MultiDelegatingModifier* modifier = dynamic_object_cast<MultiDelegatingModifier>(editObject());
	OVITO_ASSERT(!modifier || modifier->delegates() == delegates());
	OVITO_ASSERT(modifier || delegates().empty());
	OVITO_ASSERT(_delegateBoxes.size() == delegates().size());

	// Update the contents of the combo boxes.
	for(int index = 0; index < _delegateBoxes.size(); index++) {
		ModifierDelegateParameterUI::populateComboBox(_delegateBoxes[index], modifier, delegates()[index], 
			delegates()[index] ? delegates()[index]->inputDataObject() : DataObjectReference(), _delegateType);
	}
}

/******************************************************************************
* This method is called when a reference target changes.
******************************************************************************/
bool ModifierDelegateVariableListParameterUI::referenceEvent(RefTarget* source, const ReferenceEvent& event)
{
	if(source == editObject()) {
		if(event.type() == ReferenceEvent::ReferenceAdded) {
			const ReferenceFieldEvent& refevent = static_cast<const ReferenceFieldEvent&>(event);
			if(refevent.field() == &PROPERTY_FIELD(MultiDelegatingModifier::delegates)) {
				_delegates.insert(this, PROPERTY_FIELD(delegates), refevent.index(), static_object_cast<ModifierDelegate>(refevent.newTarget()));
			}
		}
		else if(event.type() == ReferenceEvent::ReferenceRemoved) {
			const ReferenceFieldEvent& refevent = static_cast<const ReferenceFieldEvent&>(event);
			if(refevent.field() == &PROPERTY_FIELD(MultiDelegatingModifier::delegates)) {
				OVITO_ASSERT(refevent.oldTarget() == _delegates[refevent.index()]);
				_delegates.remove(this, PROPERTY_FIELD(delegates), refevent.index());
			}
		}
		else if(event.type() == ReferenceEvent::ModifierInputChanged) {
			// The modifier's input from the pipeline has changed -> update list of available delegates.
			updateUI();
		}
	}
	return ParameterUI::referenceEvent(source, event);
}

/******************************************************************************
* Is called when a RefTarget has been added to a VectorReferenceField of this RefMaker.
******************************************************************************/
void ModifierDelegateVariableListParameterUI::referenceInserted(const PropertyFieldDescriptor& field, RefTarget* newTarget, int listIndex)
{
	if(field == PROPERTY_FIELD(delegates) && containerWidget()) {
		QHBoxLayout* sublayout = new QHBoxLayout();
		sublayout->setContentsMargins(0,0,0,0);
		sublayout->setSpacing(2);
		QComboBox* comboBox = new QComboBox();
		connect(comboBox, qOverload<int>(&QComboBox::activated), this, &ModifierDelegateVariableListParameterUI::onDelegateSelected);
		sublayout->addWidget(comboBox, 1);
		QToolBar* toolbar = new QToolBar();
		toolbar->setStyleSheet("QToolBar { padding: 0px; margin: 0px; border: 0px none black; }");
		toolbar->setFloatable(false);
		toolbar->setIconSize(QSize(16,16));
		QAction* removeDelegateAction = toolbar->addAction(QIcon(":/gui/actions/animation/delete_animation_key.bw.svg"), tr("Remove entry"));
		connect(removeDelegateAction, &QAction::triggered, this, &ModifierDelegateVariableListParameterUI::onRemoveDelegate);
		_removeDelegateActions.insert(listIndex, removeDelegateAction);
		_delegateBoxes.insert(listIndex, comboBox);
		sublayout->addWidget(toolbar, 0, Qt::AlignRight | Qt::AlignVCenter);
		QVBoxLayout* layout = static_cast<QVBoxLayout*>(containerWidget()->layout());
		layout->insertLayout(listIndex, sublayout);
		OVITO_ASSERT(layout->count() == delegates().size() + 1);
		ModifierDelegateParameterUI::populateComboBox(comboBox, static_object_cast<MultiDelegatingModifier>(editObject()), newTarget, 
			newTarget ? static_object_cast<ModifierDelegate>(newTarget)->inputDataObject() : DataObjectReference(), _delegateType);
		editor()->container()->updateRolloutsLater();
	}
	ParameterUI::referenceInserted(field, newTarget, listIndex);
}

/******************************************************************************
* Is called when a RefTarget has been removed from a VectorReferenceField of this RefMaker.
******************************************************************************/
void ModifierDelegateVariableListParameterUI::referenceRemoved(const PropertyFieldDescriptor& field, RefTarget* oldTarget, int listIndex)
{
	if(field == PROPERTY_FIELD(delegates) && containerWidget()) {
		QVBoxLayout* layout = static_cast<QVBoxLayout*>(containerWidget()->layout());
		QLayoutItem* item = layout->takeAt(listIndex);
		while(QLayoutItem* child = item->layout()->takeAt(0)) {
			child->widget()->deleteLater();
			delete child;
		}
		delete item;
		OVITO_ASSERT(layout->count() == delegates().size() + 1);
		_removeDelegateActions.remove(listIndex);
		_delegateBoxes.remove(listIndex);
		editor()->container()->updateRolloutsLater();
	}
	ParameterUI::referenceRemoved(field, oldTarget, listIndex);
}

/******************************************************************************
* Sets the enabled state of the UI.
******************************************************************************/
void ModifierDelegateVariableListParameterUI::setEnabled(bool enabled)
{
	if(enabled == isEnabled()) 
		return;
	ParameterUI::setEnabled(enabled);
	if(containerWidget()) 
		containerWidget()->setEnabled(editObject() && isEnabled());
}

/******************************************************************************
* Is called when the user requested the addition of a new delegate to the modifier.
******************************************************************************/
void ModifierDelegateVariableListParameterUI::onAddDelegate()
{
	if(!editObject()) return;

	undoableTransaction(tr("Add modifier input"), [&]() {
		// Add an empty delegate slot to the modifier.
		MultiDelegatingModifier* modifier = static_object_cast<MultiDelegatingModifier>(editObject());
		modifier->_delegates.push_back(modifier, PROPERTY_FIELD(MultiDelegatingModifier::delegates), nullptr);
	});
}

/******************************************************************************
* Is called when the user requested the removal of an existing delegate from the modifier.
******************************************************************************/
void ModifierDelegateVariableListParameterUI::onRemoveDelegate()
{
	// Get the QAction from which the signal originated.
	QAction* action = qobject_cast<QAction*>(sender());
	OVITO_ASSERT(action);
	if(!action || !editObject()) return;

	// Determine the list index of the delegate to be deleted.
	int index = _removeDelegateActions.indexOf(action);
	OVITO_ASSERT(index >= 0);

	undoableTransaction(tr("Remove modifier input"), [&]() {
		MultiDelegatingModifier* modifier = static_object_cast<MultiDelegatingModifier>(editObject());
		modifier->_delegates.remove(modifier, PROPERTY_FIELD(MultiDelegatingModifier::delegates), index);
	});
}

/******************************************************************************
* Is called when the user selects an entry in a combo box.
******************************************************************************/
void ModifierDelegateVariableListParameterUI::onDelegateSelected(int index)
{
	// Get the QComboBox from which the signal originated.
	QComboBox* comboBox = qobject_cast<QComboBox*>(sender());
	OVITO_ASSERT(comboBox);
	if(!comboBox || !editObject()) return;

	// Determine the list index of the delegate to be changed.
	int delegateIndex = _delegateBoxes.indexOf(comboBox);
	OVITO_ASSERT(delegateIndex >= 0);

	// Get the selected delegate type.
	OvitoClassPtr delegateType = comboBox->currentData().value<OvitoClassPtr>();
	if(!delegateType) return;

	// Get the selected data object.
	DataObjectReference ref = comboBox->currentData(Qt::UserRole + 1).value<DataObjectReference>();

	undoableTransaction(tr("Select modifier input"), [&]() {
		MultiDelegatingModifier* modifier = static_object_cast<MultiDelegatingModifier>(editObject());
		OVITO_ASSERT(modifier->delegates() == delegates());

		ModifierDelegate* oldDelegate = delegates()[delegateIndex];
		if(!oldDelegate || &oldDelegate->getOOClass() != delegateType || oldDelegate->inputDataObject() != ref) {
			// Create the new delegate object.
			OORef<ModifierDelegate> delegate = static_object_cast<ModifierDelegate>(delegateType->createInstance(modifier->dataset()));
			// Set which input data object the delegate should operate on.
			delegate->setInputDataObject(ref);
			// Activate the new delegate.
			modifier->_delegates.set(modifier, PROPERTY_FIELD(MultiDelegatingModifier::delegates), delegateIndex, std::move(delegate));
			OVITO_ASSERT(delegates()[delegateIndex] == modifier->delegates()[delegateIndex]);
		}
		Q_EMIT valueEntered();
	});
}

}	// End of namespace
