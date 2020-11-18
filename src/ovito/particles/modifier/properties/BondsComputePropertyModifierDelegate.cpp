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

#include <ovito/particles/Particles.h>
#include <ovito/particles/objects/BondsObject.h>
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include <ovito/core/utilities/concurrent/ParallelFor.h>
#include <ovito/core/dataset/DataSet.h>
#include "BondsComputePropertyModifierDelegate.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(BondsComputePropertyModifierDelegate);

/******************************************************************************
* Indicates which data objects in the given input data collection the modifier
* delegate is able to operate on.
******************************************************************************/
QVector<DataObjectReference> BondsComputePropertyModifierDelegate::OOMetaClass::getApplicableObjects(const DataCollection& input) const
{
    if(const ParticlesObject* particles = input.getObject<ParticlesObject>()) {
        if(particles->bonds())
       		return { DataObjectReference(&ParticlesObject::OOClass()) };
    }
    return {};
}

/******************************************************************************
* Creates and initializes a computation engine that will compute the
* modifier's results.
******************************************************************************/
std::shared_ptr<ComputePropertyModifierDelegate::PropertyComputeEngine> BondsComputePropertyModifierDelegate::createEngine(
				TimePoint time,
				const PipelineFlowState& input,
				const ConstDataObjectPath& containerPath,
				PropertyPtr outputProperty,
				ConstPropertyPtr selectionProperty,
				QStringList expressions)
{
	// Create engine object. Pass all relevant modifier parameters to the engine as well as the input data.
	return std::make_shared<Engine>(
			input.stateValidity(),
			time,
			std::move(outputProperty),
			containerPath,
			std::move(selectionProperty),
			std::move(expressions),
			dataset()->animationSettings()->timeToFrame(time),
			input);
}

/******************************************************************************
* Constructor.
******************************************************************************/
BondsComputePropertyModifierDelegate::Engine::Engine(
		const TimeInterval& validityInterval,
		TimePoint time,
		PropertyPtr outputProperty,
		const ConstDataObjectPath& containerPath,
		ConstPropertyPtr selectionProperty,
		QStringList expressions,
		int frameNumber,
		const PipelineFlowState& input) :
	ComputePropertyModifierDelegate::PropertyComputeEngine(
			validityInterval,
			time,
			input,
			containerPath,
			std::move(outputProperty),
			std::move(selectionProperty),
			std::move(expressions),
			frameNumber,
			std::make_unique<BondExpressionEvaluator>()),
	_inputFingerprint(input.expectObject<ParticlesObject>())
{
}

/******************************************************************************
* Performs the actual computation. This method is executed in a worker thread.
******************************************************************************/
void BondsComputePropertyModifierDelegate::Engine::perform()
{
	setProgressText(tr("Computing property '%1'").arg(outputProperty()->name()));
	setProgressMaximum(outputProperty()->size());
	setProgressValue(0);

	// Parallelized loop over all bonds.
	parallelForChunks(outputProperty()->size(), *this, [this](size_t startIndex, size_t count, Task& promise) {
		BondExpressionEvaluator::Worker worker(*_evaluator);

		size_t endIndex = startIndex + count;
		size_t componentCount = outputProperty()->componentCount();
		for(size_t bondIndex = startIndex; bondIndex < endIndex; bondIndex++) {

			// Update progress indicator.
			if((bondIndex % 1024) == 0)
				promise.incrementProgressValue(1024);

			// Exit if operation was canceled.
			if(promise.isCanceled())
				return;

			// Skip unselected bonds if requested.
			if(selectionArray() && !selectionArray()[bondIndex])
				continue;

			for(size_t component = 0; component < componentCount; component++) {

				// Compute expression value.
				FloatType value = worker.evaluate(bondIndex, component);

				// Store results in output property.
				outputArray().set(bondIndex, component, value);
			}
		}
	});

	// Release data that is no longer needed to reduce memory footprint.
	releaseWorkingData();
}

/******************************************************************************
* Injects the computed results of the engine into the data pipeline.
******************************************************************************/
void BondsComputePropertyModifierDelegate::Engine::applyResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	if(_inputFingerprint.hasChanged(state.expectObject<ParticlesObject>()))
		modApp->throwException(tr("Cached modifier results are obsolete, because the number or the storage order of input particles has changed."));

	PropertyComputeEngine::applyResults(time, modApp, state);
}

}	// End of namespace
}	// End of namespace
