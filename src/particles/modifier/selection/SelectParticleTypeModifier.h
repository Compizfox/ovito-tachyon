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

#ifndef __OVITO_SELECT_PARTICLE_TYPE_MODIFIER_H
#define __OVITO_SELECT_PARTICLE_TYPE_MODIFIER_H

#include <core/Core.h>
#include "../ParticleModifier.h"

namespace Particles {

using namespace Ovito;

class ParticlePropertyComboBox;

/******************************************************************************
* Selects particles of of one or more types.
******************************************************************************/
class SelectParticleTypeModifier : public ParticleModifier
{
public:

	/// Default constructor.
	Q_INVOKABLE SelectParticleTypeModifier() : _inputPropertyRef(ParticleProperty::ParticleTypeProperty) {}

	/// This virtual method is called by the system when the modifier has been inserted into a PipelineObject.
	virtual void initializeModifier(PipelineObject* pipelineObject, ModifierApplication* modApp) override;

	/// Returns the particle type property that is used as source for the selection.
	const ParticlePropertyReference& sourceProperty() const { return _inputPropertyRef; }

	/// Sets the particle type property that is used as source for the selection.
	void setSourceProperty(const ParticlePropertyReference& prop);

	/// Returns the list of particle type identifiers that are selected.
	const QSet<int>& selectedParticleTypes() const { return _selectedParticleTypes; }

	/// Sets the list of particle type identifiers to select.
	void setSelectedParticleTypes(const QSet<int>& types);

	/// Sets a single particle type identifier to be selected.
	void setSelectedParticleType(int type) { QSet<int> list; list.insert(type); setSelectedParticleTypes(list); }

	/// Retrieves the selected input type property from the given modifier input state.
	ParticleTypeProperty* lookupInputProperty(const PipelineFlowState& inputState) const;

public:

	Q_PROPERTY(Particles::ParticlePropertyReference sourceProperty READ sourceProperty WRITE setSourceProperty)

protected:

	/// Saves the class' contents to the given stream.
	virtual void saveToStream(ObjectSaveStream& stream) override;

	/// Loads the class' contents from the given stream.
	virtual void loadFromStream(ObjectLoadStream& stream) override;

	/// Creates a copy of this object.
	virtual OORef<RefTarget> clone(bool deepCopy, CloneHelper& cloneHelper) override;

	/// Modifies the particle object.
	virtual ObjectStatus modifyParticles(TimePoint time, TimeInterval& validityInterval) override;

	/// The particle type property that is used as source for the selection.
	ParticlePropertyReference _inputPropertyRef;

	/// The identifiers of the particle types to select.
	QSet<int> _selectedParticleTypes;

private:

	Q_OBJECT
	OVITO_OBJECT

	Q_CLASSINFO("DisplayName", "Select particle type");
	Q_CLASSINFO("ModifierCategory", "Selection");
};

/******************************************************************************
* A properties editor for the SelectParticleTypeModifier class.
******************************************************************************/
class SelectParticleTypeModifierEditor : public ParticleModifierEditor
{
public:

	/// Default constructor
	Q_INVOKABLE SelectParticleTypeModifierEditor() {}

protected:

	/// Creates the user interface controls for the editor.
	virtual void createUI(const RolloutInsertionParameters& rolloutParams) override;

protected Q_SLOTS:

	/// Updates the contents of the property list combo box.
	void updatePropertyList();

	/// Updates the contents of the particle type list box.
	void updateParticleTypeList();

	/// This is called when the user has selected another item in the particle property list.
	void onPropertySelected(int index);

	/// This is called when the user has selected another particle type.
	void onParticleTypeSelected(QListWidgetItem* item);

protected:

	/// This method is called when a reference target changes.
	virtual bool referenceEvent(RefTarget* source, ReferenceEvent* event) override;

private:

	/// The list of particle type properties.
	ParticlePropertyComboBox* propertyListBox;

	/// The list of particle types.
	QListWidget* particleTypesBox;

	Q_OBJECT
	OVITO_OBJECT
};

};	// End of namespace

#endif // __OVITO_SELECT_PARTICLE_TYPE_MODIFIER_H