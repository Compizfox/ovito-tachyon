///////////////////////////////////////////////////////////////////////////////
// 
//  Copyright (2013) Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  OVITO is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

/** 
 * \file ParameterUI.h 
 * \brief Contains the definition of the Ovito::ParameterUI class and some derived classes.
 */

#ifndef __OVITO_PARAMETER_UI_H
#define __OVITO_PARAMETER_UI_H

#include <core/Core.h>
#include <core/reference/RefTarget.h>
#include "PropertiesEditor.h"

namespace Ovito {

/**
 * \brief Base class for UI components that allow the user to edit a parameter
 *        of a RefTarget derived object in the PropertiesEditor.
 */
class ParameterUI : public RefMaker
{
public:

	/// \brief Constructor.
	/// \param parentEditor The editor in which this parameter UI is used. This becomes the parent of this object.
	///
	/// The parameter UI is automatically deleted when the editor is deleted. 
	ParameterUI(QObject* parent);

	/// \brief Destructor.
	/// \note The destructor must clear all references since this UI object
	///       might have been deleted via the delete operator and not via
	///       the OvitoObject::autoDeleteObject() method which should normally be used for PluginClass derived classes.
	virtual ~ParameterUI() { clearAllReferences(); }	
	
	/// \brief Gets the object whose parameter is being edited/shown in this parameter UI.
	/// \return The current object being edited.
	/// \sa setEditObject()
	RefTarget* editObject() const { return _editObject; }
	
	/// \brief Returns a pointer to the properties editor this parameter UI belongs to.
	/// \return The editor in which this parameter UI is used or NULL if the parameter UI is used outside of a PropertiesEditor.
	PropertiesEditor* editor() const { return qobject_cast<PropertiesEditor*>(this->parent()); }

	/// \brief Returns the enabled state of the UI.
	/// \return \c true if this parameter's value can be changed by the user;
	///         \c false otherwise.
	/// \sa setEnabled()
	bool isEnabled() const { return _enabled; }

	/// \brief Returns the disabled state of the UI. This is just the inverse of the enabled state.
	/// \return \c false if this parameter's value can be changed by the user;
	///         \c true otherwise.
	/// \sa isEnabled()
	bool isDisabled() const { return !isEnabled(); }

public:	
	
	Q_PROPERTY(RefTarget editObject READ editObject)		
	Q_PROPERTY(bool isEnabled READ isEnabled WRITE setEnabled)
	Q_PROPERTY(bool isDisabled READ isDisabled WRITE setDisabled)
		
public Q_SLOTS:	
		
	/// \brief This method is called when a new editable object has been assigned to the properties owner 
	///       this parameter UI belongs to.
	/// 
	/// The parameter UI should react to this change appropriately and
	/// show the properties value for the new edit object in the UI. The default implementation
	/// of this method just calls updateUI() to reflect the change.
	///
	/// \sa setEditObject()
	virtual void resetUI() { updateUI(); }
	
	/// \brief This method updates the displayed value of the parameter UI.
	/// 
	/// This method should be overridden by derived classes.
	virtual void updateUI() {}

	/// \brief Sets the enabled state of the UI.
	/// \param enabled Controls whether may change the parameter's value or not.
	/// \sa isEnabled()
	virtual void setEnabled(bool enabled) { _enabled = enabled; }

	/// \brief Sets the enabled state of the UI. This is just the reverse of setEnabled().
	/// \param disabled Controls whether may change the parameter's value or not.
	/// \sa setEnabled()
	/// \sa isDisabled()
	void setDisabled(bool disabled) { setEnabled(!disabled); }

	/// \brief Sets the object whose property is being displayed in this parameter UI.
	/// \sa editObject()
	void setEditObject(RefTarget* newObject) {
		_editObject = newObject;
		resetUI();
	}
	
private:

	/// The object whose parameter is being edited.
	ReferenceField<RefTarget> _editObject;

	/// Stores whether this UI is enabled.
	bool _enabled;

	Q_OBJECT
	OVITO_OBJECT

	DECLARE_REFERENCE_FIELD(_editObject);
};

/**
 * \brief Base class for UI components that allow the user to edit a property of
 *        an object that is stored in a reference field, a property field, or a Qt property.
 */
class PropertyParameterUI : public ParameterUI
{
public:
	
	/// \brief Constructor for a Qt property.
	/// \param parent The editor in which this parameter UI is used. This becomes the parent of this object.
	/// \param propertyName The name of the property that has been defined using the \c Q_PROPERTY macro. 
	PropertyParameterUI(QObject* parent, const char* propertyName);

	/// \brief Constructor for a PropertyField or ReferenceField.
	/// \param parent The editor in which this parameter UI is used. This becomes the parent of this object.
	/// \param propField The property or reference field.
	PropertyParameterUI(QObject* parent, const PropertyFieldDescriptor& propField);
	
	/// \brief Destructor.
	/// \note The destructor must clear all references since this UI object
	///       might have been deleted via the delete operator and not via
	///       the OvitoObject::autoDeleteObject() method which should normally be used for OvitoObject-derived classes.
	virtual ~PropertyParameterUI() { clearAllReferences(); }

	/// \brief Returns the property being edited in this parameter UI.
	/// \return The name of the QObject property this UI is bound to or \c
	///         NULL if this PropertyUI is bound to a PropertyField.
	/// \sa propertyField()
	const char* propertyName() const { return _propertyName; }

	/// \brief Returns the property or reference field being edited.
	/// \return A pointer to the descriptor of the PropertyField or ReferenceField being edited or
	///         \c NULL if this PropertyUI is bound to a normal Qt property.
	/// \sa propertyName()
	const PropertyFieldDescriptor* propertyField() const { return _propField; }
	
	/// \brief Indicates whether this parameter UI is representing a sub-object property (e.g. an animation controller).
	bool isReferenceFieldUI() const { return (_propField && _propField->isReferenceField()); }

	/// \brief Indicates whether this parameter UI is representing a PropertyField based property.
	bool isPropertyFieldUI() const { return (_propField && !_propField->isReferenceField()); }

	/// \brief Indicates whether this parameter UI is representing a Qt property.
	bool isQtPropertyUI() const { return _propField == nullptr; }

	/// \brief Returns the sub-object that is bound to this parameter UI.
	/// \return The object stored in the reference field. This may be \c NULL either when there is no editable object selected in the parent editor
	///         or if the editable object's reference field is currently empty.
	RefTarget* parameterObject() const { return _parameterObject; }

	/// \brief This method is called when parameter object has been assigned to the reference field of the editable object
	/// this parameter UI is bound to.
	///
	/// It is also called when the editable object itself has
	/// been replaced in the editor. The parameter UI should react to this change appropriately and
	/// show the properties value for the new edit object in the UI. New implementations of this
	/// method must call the base implementation before any other action is taken.
	virtual void resetUI() override;

public:
	
	Q_PROPERTY(const char* propertyName READ propertyName)
	Q_PROPERTY(RefTarget parameterObject READ parameterObject)

protected:

	/// This method is called when a reference target changes.
	virtual bool referenceEvent(RefTarget* source, ReferenceEvent* event) override;

private:
	
	/// The controller or sub-object whose value is being edited.
	ReferenceField<RefTarget> _parameterObject;

	/// The property or reference field being edited or NULL if editing a Qt property.
	const PropertyFieldDescriptor* _propField;

	/// The name of the Qt property being edited or NULL.
	const char* _propertyName;
	
	Q_OBJECT
	OVITO_OBJECT

	DECLARE_REFERENCE_FIELD(_parameterObject);
};

};

#endif // __OVITO_PARAMETER_UI_H