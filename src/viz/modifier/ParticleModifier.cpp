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

#include <core/Core.h>
#include <core/scene/pipeline/ModifierApplication.h>
#include "ParticleModifier.h"

namespace Viz {

IMPLEMENT_SERIALIZABLE_OVITO_OBJECT(Viz, ParticleModifier, Modifier)
IMPLEMENT_OVITO_OBJECT(Viz, ParticleModifierEditor, PropertiesEditor)

/******************************************************************************
* This modifies the input object.
******************************************************************************/
ObjectStatus ParticleModifier::modifyObject(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	// This method is not re-entrant. If this method is called while the modifier is already being
	// evaluated then we are not able to process the request.
	if(!_input.isEmpty())
		return ObjectStatus(ObjectStatus::Error, tr("Cannot handle re-entrant modifier calls."));

	// Prepare internal fields.
	_input = state;
	_output = state;
	_modApp = modApp;
	ObjectStatus status;

	try {
		ParticlePropertyObject* posProperty = inputStandardProperty(ParticleProperty::PositionProperty);
		if(!posProperty)
			throw Exception(tr("This modifier cannot be evaluated because the input does not contain any particles."));
		_particleCount = posProperty->size();

		// Let the derived class do the actual work.
		TimeInterval validityInterval = state.stateValidity();
		status = modifyParticles(time, validityInterval);

		// Put result into geometry pipeline.
		state = _output;
		state.intersectStateValidity(validityInterval);
	}
	catch(const Exception& ex) {
		ex.logError();
		// Transfer exception message to evaluation status.
		status = ObjectStatus(ObjectStatus::Error, ex.messages().join('\n'));
	}

	// Cleanup
	_cloneHelper.reset();
	_input.clear();
	_output.clear();
	_modApp = nullptr;

	return status;
}

/******************************************************************************
* Returns a standard particle property from the input state.
******************************************************************************/
ParticlePropertyObject* ParticleModifier::inputStandardProperty(ParticleProperty::Type which) const
{
	OVITO_ASSERT(which != ParticleProperty::UserProperty);
	for(const auto& o : _input.objects()) {
		ParticlePropertyObject* property = dynamic_object_cast<ParticlePropertyObject>(o.get());
		if(property && property->type() == which) return property;
	}
	return nullptr;
}

/******************************************************************************
* Returns the property with the given identifier from the input object.
******************************************************************************/
ParticlePropertyObject* ParticleModifier::expectCustomProperty(const QString& propertyName, int dataType, size_t componentCount) const
{
	for(const auto& o : _input.objects()) {
		ParticlePropertyObject* property = dynamic_object_cast<ParticlePropertyObject>(o.get());
		if(property && property->name() == propertyName) {
			if(property->dataType() != dataType)
				throw Exception(tr("The modifier cannot be evaluated because the particle property '%1' does not have the required data type.").arg(property->name()));
			if(property->componentCount() != componentCount)
				throw Exception(tr("The modifier cannot be evaluated because the particle property '%1' does not have the required number of components per particle.").arg(property->name()));
			return property;
		}
	}
	throw Exception(tr("The modifier cannot be evaluated because the input does not contain the required particle property (name: %1).").arg(propertyName));
}

/******************************************************************************
* Returns the given standard channel from the input object.
* The returned channel may not be modified. If they input object does
* not contain the standard channel then an exception is thrown.
******************************************************************************/
ParticlePropertyObject* ParticleModifier::expectStandardProperty(ParticleProperty::Type which) const
{
	ParticlePropertyObject* property = inputStandardProperty(which);
	if(!property)
		throw Exception(tr("The modifier cannot be evaluated because the input does not contain the required particle property '%1'.").arg(ParticleProperty::standardPropertyName(which)));
	return property;
}

/******************************************************************************
* Creates a standard particle in the modifier's output.
* If the particle property already exists in the input, its contents are copied to the
* output property by this method.
******************************************************************************/
ParticlePropertyObject* ParticleModifier::outputStandardProperty(ParticleProperty::Type which)
{
	// Check if property already exists in the input.
	OORef<ParticlePropertyObject> inputProperty = inputStandardProperty(which);

	// Check if property already exists in the output.
	OORef<ParticlePropertyObject> outputProperty;
	for(const auto& o : _output.objects()) {
		ParticlePropertyObject* property = dynamic_object_cast<ParticlePropertyObject>(o.get());
		if(property && property->type() == which) {
			outputProperty = property;
			break;
		}
	}

	if(outputProperty) {
		// Is the existing output property still a shallow copy of the input?
		if(outputProperty == inputProperty) {
			// Make a real copy of the property, which may be modified.
			outputProperty = cloneHelper()->cloneObject(inputProperty, false);
			_output.replaceObject(inputProperty.get(), inputProperty);
		}
	}
	else {
		// Create a new particle property in the output.
		outputProperty = ParticlePropertyObject::create(_particleCount, which);
		_output.addObject(outputProperty.get());
	}

	return outputProperty.get();
}


#if 0
/******************************************************************************
* Creates a shallow copy of the given input object on demand and returns it.
******************************************************************************/
AtomsObject* ParticleModifier::output()
{
	OVITO_ASSERT(inputAtoms != NULL);
	if(outputAtoms) return outputAtoms.get();

	// Make a shallow copy.
	outputAtoms = cloneHelper()->cloneObject(inputAtoms, false);
	return outputAtoms.get();
}


#endif

/******************************************************************************
* Saves the class' contents to the given stream.
******************************************************************************/
void ParticleModifier::saveToStream(ObjectSaveStream& stream)
{
	Modifier::saveToStream(stream);
	stream.beginChunk(0x01);
	// For future use...
	stream.endChunk();
}

/******************************************************************************
* Loads the class' contents from the given stream.
******************************************************************************/
void ParticleModifier::loadFromStream(ObjectLoadStream& stream)
{
	Modifier::loadFromStream(stream);
	stream.expectChunk(0x01);
	// For future use...
	stream.closeChunk();
}

/******************************************************************************
* This handler is called when a new edit object has been loaded into the editor.
******************************************************************************/
void ParticleModifierEditor::onContentsReplaced(RefTarget* newEditObject)
{
	ModifierApplication* modApp = nullptr;
	Modifier* mod = dynamic_object_cast<Modifier>(newEditObject);
	if(mod && mod->modifierApplications().empty() == false)
		modApp = mod->modifierApplications().front();
	_modAppListener.setTarget(modApp);

	updateStatusLabel(modApp);
}

/******************************************************************************
* This handler is called when the current ModifierApplication sends a
* notification event.
******************************************************************************/
void ParticleModifierEditor::onModAppNotificationEvent(ReferenceEvent* event)
{
	if(event->type() == ReferenceEvent::StatusChanged) {
		ModifierApplication* modApp = dynamic_object_cast<ModifierApplication>(event->sender());
		updateStatusLabel(modApp);
	}
}

/******************************************************************************
* Updates the text of the result label.
******************************************************************************/
void ParticleModifierEditor::updateStatusLabel(ModifierApplication* modApp)
{
	if(!_statusTextLabel || !_statusIconLabel) return;

	if(modApp != NULL) {
		_statusTextLabel->setText(modApp->status().longText());
		if(modApp->status().type() == ObjectStatus::Success) {
			if(modApp->status().longText().isEmpty() == false)
				_statusIconLabel->setPixmap(_modifierStatusInfoIcon);
			else
				_statusIconLabel->clear();
		}
		else if(modApp->status().type() == ObjectStatus::Warning)
			_statusIconLabel->setPixmap(_modifierStatusWarningIcon);
		else if(modApp->status().type() == ObjectStatus::Error)
			_statusIconLabel->setPixmap(_modifierStatusErrorIcon);
		else
			_statusIconLabel->clear();
	}
	else {
		_statusTextLabel->clear();
		_statusIconLabel->clear();
	}
}

/******************************************************************************
* Returns a widget that displays a message sent by the modifier that
* states the outcome of the modifier evaluation. Derived classes of this
* editor base class can add the widget to their user interface.
******************************************************************************/
QWidget* ParticleModifierEditor::statusLabel()
{
	if(_statusLabel) return _statusLabel;

	_statusLabel = new QWidget();
	QGridLayout* layout = new QGridLayout(_statusLabel);
	layout->setContentsMargins(0,0,0,0);
	layout->setColumnStretch(1, 1);
	//_statusLabel->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);

	_statusIconLabel = new QLabel(_statusLabel);
	_statusIconLabel->setAlignment(Qt::AlignTop);
	layout->addWidget(_statusIconLabel, 0, 0, Qt::AlignTop);

	_statusTextLabel = new QLabel(_statusLabel);
	_statusTextLabel->setAlignment(Qt::AlignTop);
	//_statusTextLabel->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
	_statusTextLabel->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard | Qt::LinksAccessibleByMouse | Qt::LinksAccessibleByKeyboard);
	_statusTextLabel->setWordWrap(true);
	layout->addWidget(_statusTextLabel, 0, 1);

	return _statusLabel;
}

};	// End of namespace
