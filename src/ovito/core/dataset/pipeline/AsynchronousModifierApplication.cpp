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

#include <ovito/core/Core.h>
#include <ovito/core/dataset/pipeline/AsynchronousModifierApplication.h>

namespace Ovito {

IMPLEMENT_OVITO_CLASS(AsynchronousModifierApplication);
SET_MODIFIER_APPLICATION_TYPE(AsynchronousModifier, AsynchronousModifierApplication);

/******************************************************************************
* Constructor.
******************************************************************************/
AsynchronousModifierApplication::AsynchronousModifierApplication(DataSet* dataset) : ModifierApplication(dataset)
{
}

/******************************************************************************
* Is called when a RefTarget referenced by this object has generated an event.
******************************************************************************/
bool AsynchronousModifierApplication::referenceEvent(RefTarget* source, const ReferenceEvent& event)
{
	if(event.type() == ReferenceEvent::TargetEnabledOrDisabled && source == modifier()) {
		// Throw away cached results when the modifier is being disabled.
		_validStages.clear();
		_completedEngine.reset();
	}
	else if(event.type() == ReferenceEvent::PreliminaryStateAvailable && source == input()) {
		// Throw away cached results when the modifier's input changes, unless the engine requests otherwise.
		if(_completedEngine && !_completedEngine->pipelineInputChanged())
			_completedEngine.reset();
	}
	else if((event.type() == ReferenceEvent::TargetChanged && source == input()) || (event.type() == ReferenceEvent::ModifierInputChanged && source == modifier())) {
		// Whenever the modifier's inputs change, invalidate the cached computation results hold on to any 
		// cached results needed for preliminary pipeline evaluation. 
		_validStages.clear();
		if(_completedEngine)
			_completedEngine->setValidityInterval(TimeInterval::empty());
	}
	else if(event.type() == ReferenceEvent::TargetChanged && source == modifier()) {
		// Whenever the modifier changes, invalidate the cached computation results
		// unless the engine requests otherwise.
		for(auto e = _validStages.begin(); e != _validStages.end(); ++e) {
			if(!(*e)->modifierChanged(static_cast<const PropertyFieldEvent&>(event))) {
				_validStages.erase(e, _validStages.end());
				if(_completedEngine)
					_completedEngine->setValidityInterval(TimeInterval::empty());
				break;
			}
		}
	}
	return ModifierApplication::referenceEvent(source, event);
}

/******************************************************************************
* Gets called when the data object of the node has been replaced.
******************************************************************************/
void AsynchronousModifierApplication::referenceReplaced(const PropertyFieldDescriptor& field, RefTarget* oldTarget, RefTarget* newTarget)
{
	// Throw away cached results when the modifier is being detached from this ModifierApplication.
	if(field == PROPERTY_FIELD(modifier)) {
		_validStages.clear();
		_completedEngine.reset();
	}
	ModifierApplication::referenceReplaced(field, oldTarget, newTarget);
}

}	// End of namespace
