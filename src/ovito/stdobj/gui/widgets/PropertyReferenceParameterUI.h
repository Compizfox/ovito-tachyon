////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2017 Alexander Stukowski
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


#include <ovito/stdobj/gui/StdObjGui.h>
#include <ovito/stdobj/gui/widgets/PropertySelectionComboBox.h>
#include <ovito/stdobj/properties/PropertyContainer.h>
#include <ovito/gui/desktop/properties/ParameterUI.h>

namespace Ovito { namespace StdObj {

/**
 * \brief This parameter UI lets the user select a property.
 */
class OVITO_STDOBJGUI_EXPORT PropertyReferenceParameterUI : public PropertyParameterUI
{
	Q_OBJECT
	OVITO_CLASS(PropertyReferenceParameterUI)

public:

	/// Constructor.
	PropertyReferenceParameterUI(QObject* parentEditor, const char* propertyName, PropertyContainerClassPtr containerClass, bool showComponents = true, bool inputProperty = true);

	/// Constructor.
	PropertyReferenceParameterUI(QObject* parentEditor, const PropertyFieldDescriptor& propField, PropertyContainerClassPtr containerClass, bool showComponents = true, bool inputProperty = true);

	/// Destructor.
	virtual ~PropertyReferenceParameterUI();

	/// This returns the combo box managed by this ParameterUI.
	QComboBox* comboBox() const { return _comboBox; }

	/// This method is called when a new editable object has been assigned to the properties owner this
	/// parameter UI belongs to.
	virtual void resetUI() override;

	/// This method updates the displayed value of the property UI.
	virtual void updateUI() override;

	/// Sets the enabled state of the UI.
	virtual void setEnabled(bool enabled) override;

	/// Sets the tooltip text for the combo box widget.
	void setToolTip(const QString& text) const {
		if(comboBox()) comboBox()->setToolTip(text);
	}

	/// Sets the What's This helper text for the combo box.
	void setWhatsThis(const QString& text) const {
		if(comboBox()) comboBox()->setWhatsThis(text);
	}

	/// Returns the property container from which the user can select a property.
	const PropertyContainerReference& containerRef() const { return _containerRef; }

	/// Sets the property container from which the user can select a property.
	void setContainerRef(const PropertyContainerReference& containerRef) {
		if(_containerRef != containerRef) {
			_containerRef = containerRef;
			_comboBox->setContainerClass(_containerRef.dataClass());
			updateUI();
		}
	}

	/// Installs optional callback function that allows clients to filter the displayed property list.
	void setPropertyFilter(std::function<bool(const PropertyObject*)> filter) {
		_propertyFilter = std::move(filter);
	}

public:

	Q_PROPERTY(QComboBox comboBox READ comboBox);

public Q_SLOTS:

	/// Takes the value entered by the user and stores it in the property field
	/// this property UI is bound to.
	void updatePropertyValue();

protected:

	/// This method is called when a reference target changes.
	virtual bool referenceEvent(RefTarget* source, const ReferenceEvent& event) override;

	/// Returns the value currently set for the property field.
	PropertyReference getPropertyReference();

	/// Populates the combox box with items.
	void addItemsToComboBox(const PipelineFlowState& state);

protected:

	/// The combo box of the UI component.
	QPointer<PropertySelectionComboBox> _comboBox;

	/// Controls whether the combo box should display a separate entry for each component of a property.
	bool _showComponents;

	/// Controls whether the combo box should list input or output properties.
	bool _inputProperty;

	/// The container from which properties that can be selected.
	PropertyContainerReference _containerRef;

	/// An optional callback function that allows clients to filter the displayed property list.
	std::function<bool(const PropertyObject*)> _propertyFilter;
};

}	// End of namespace
}	// End of namespace
