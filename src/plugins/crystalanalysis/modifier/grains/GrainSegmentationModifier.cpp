///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2018) Alexander Stukowski
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

#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <plugins/stdobj/simcell/SimulationCellObject.h>
#include <plugins/particles/objects/BondsVis.h>
//#include <plugins/particles/modifier/ParticleInputHelper.h>
//#include <plugins/particles/modifier/ParticleOutputHelper.h>
#include <core/utilities/concurrent/Promise.h>
#include <core/utilities/concurrent/TaskManager.h>
#include <core/utilities/units/UnitsManager.h>
#include "GrainSegmentationModifier.h"
#include "GrainSegmentationEngine.h"
//#include "GrainTrackingEngine.h"

#include <plugins/particles/Particles.h>
#include <plugins/stdobj/properties/PropertyStorage.h>
#include <plugins/stdobj/series/DataSeriesObject.h>
#include <core/dataset/pipeline/ModifierApplication.h>
#include <core/dataset/DataSet.h>

#include <ptm/ptm_functions.h>

namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

IMPLEMENT_OVITO_CLASS(GrainSegmentationModifier);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, rmsdCutoff);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, mergingThreshold);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, minGrainAtomCount);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, onlySelectedParticles);
DEFINE_REFERENCE_FIELD(GrainSegmentationModifier, bondsVis);
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, rmsdCutoff, "RMSD cutoff");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, minGrainAtomCount, "Minimum grain size (# of atoms)");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, mergingThreshold, "Merging threshold");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, onlySelectedParticles, "Use only selected particles");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, bondsVis, "Bonds display");
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(GrainSegmentationModifier, rmsdCutoff, FloatParameterUnit, 0);
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(GrainSegmentationModifier, minGrainAtomCount, IntegerParameterUnit, 1);
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(GrainSegmentationModifier, mergingThreshold, FloatParameterUnit, 0);

//IMPLEMENT_OVITO_CLASS(GrainSegmentationModifierApplication);
//SET_MODIFIER_APPLICATION_TYPE(GrainSegmentationModifier, GrainSegmentationModifierApplication);

/******************************************************************************
* Constructs the modifier object.
******************************************************************************/
GrainSegmentationModifier::GrainSegmentationModifier(DataSet* dataset) : StructureIdentificationModifier(dataset),
		_rmsdCutoff(0.1),
		_minGrainAtomCount(100),
		_onlySelectedParticles(false),
		_mergingThreshold(1.5)
{
	// Define the structure types.
	createStructureType(PTMAlgorithm::OTHER, ParticleType::PredefinedStructureType::OTHER);
	createStructureType(PTMAlgorithm::FCC, ParticleType::PredefinedStructureType::FCC);
	createStructureType(PTMAlgorithm::HCP, ParticleType::PredefinedStructureType::HCP);
	createStructureType(PTMAlgorithm::BCC, ParticleType::PredefinedStructureType::BCC);
	createStructureType(PTMAlgorithm::ICO, ParticleType::PredefinedStructureType::ICO)->setEnabled(false);
	createStructureType(PTMAlgorithm::SC, ParticleType::PredefinedStructureType::SC)->setEnabled(false);
	createStructureType(PTMAlgorithm::CUBIC_DIAMOND, ParticleType::PredefinedStructureType::CUBIC_DIAMOND)->setEnabled(false);
	createStructureType(PTMAlgorithm::HEX_DIAMOND, ParticleType::PredefinedStructureType::HEX_DIAMOND)->setEnabled(false);
	createStructureType(PTMAlgorithm::GRAPHENE, ParticleType::PredefinedStructureType::GRAPHENE)->setEnabled(false);

	// Create the visual element for the bonds (display is disabled by default).
	setBondsVis(new BondsVis(dataset));
	bondsVis()->setEnabled(false);
}

/******************************************************************************
* Creates a computation engine that finds the grains in a single frame.
******************************************************************************/
std::shared_ptr<GrainSegmentationEngine> GrainSegmentationModifier::createSegmentationEngine(TimePoint time, ModifierApplication* modApp, const PipelineFlowState& input)
{
	if(structureTypes().size() != PTMAlgorithm::NUM_STRUCTURE_TYPES)
		throwException(tr("The number of structure types has changed. Please remove this modifier from the modification pipeline and insert it again."));

	// Get modifier input.
	const ParticlesObject* particles = input.expectObject<ParticlesObject>();
	const PropertyObject* posProperty = particles->expectProperty(ParticlesObject::PositionProperty);
	const SimulationCellObject* simCell = input.expectObject<SimulationCellObject>();
	if(simCell->is2D())
		throwException(tr("The grain segmentation modifier does not support 2d simulation cells."));

	// Get particle selection.
	ConstPropertyPtr selectionProperty;
	if(onlySelectedParticles())
		selectionProperty = particles->expectProperty(ParticlesObject::SelectionProperty)->storage();

	// Initialize PTM library.
	ptm_initialize_global();

	// Create engine object. Pass all relevant modifier parameters to the engine as well as the input data.
	return std::make_shared<GrainSegmentationEngine>(particles, posProperty->storage(), simCell->data(),
			getTypesToIdentify(PTMAlgorithm::NUM_STRUCTURE_TYPES), std::move(selectionProperty),
			rmsdCutoff(), mergingThreshold(), minGrainAtomCount());
}

/******************************************************************************
* Creates a computation engine that will compute the modifier's results.
******************************************************************************/
Future<AsynchronousModifier::ComputeEnginePtr> GrainSegmentationModifier::createEngine(TimePoint time, ModifierApplication* modApp, const PipelineFlowState& input)
{
	return createSegmentationEngine(time, modApp, input);
}

/******************************************************************************
* Injects the computed results of the engine into the data pipeline.
******************************************************************************/
void GrainSegmentationEngine::emitResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	StructureIdentificationEngine::emitResults(time, modApp, state);

	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(modApp->modifier());
	OVITO_ASSERT(modifier);
	ParticlesObject* particles = state.expectMutableObject<ParticlesObject>();

	// Output per-particle properties.
	if(rmsd())
		particles->createProperty(rmsd());
	if(orientations())
		particles->createProperty(orientations());
	if(atomClusters())
		particles->createProperty(atomClusters());

#if 0
//TODO: put this back in

	// Output lattice neighbor bonds.
	if(latticeNeighborBonds()) {
		for(int i = output.objects().size()-1; i>=0; i--)
			if(dynamic_object_cast<BondProperty>(output.objects()[i]))
				output.removeObjectByIndex(i);
		OORef<BondProperty> topologyPropertyObj = BondProperty::createFromStorage(modifier->dataset(), latticeNeighborBonds());
		topologyPropertyObj->setVisElement(modifier->bondsVis());
		output.addObject(topologyPropertyObj);
		poh.setOutputBondCount(latticeNeighborBonds()->size());
		poh.outputProperty<BondProperty>(bondPBCShiftVectors());
		poh.outputProperty<BondProperty>(neighborDisorientationAngles());
	}
#endif

	// Store the RMSD histogram in the ModifierApplication.
	//static_object_cast<GrainSegmentationModifierApplication>(modApp)->setRmsdHistogram(rmsdHistogramData(), rmsdHistogramBinSize());

//TODO: put this back in
	//DataSeriesObject* seriesObj = state.createObject<DataSeriesObject>(QStringLiteral("grains-rmsd"), modApp, DataSeriesObject::Line, GrainSegmentationModifier::tr("RMSD distribution"), rmsdHistogram());
	//seriesObj->setAxisLabelX(GrainSegmentationModifier::tr("RMSD"));
	//seriesObj->setIntervalStart(0);
	//seriesObj->setIntervalEnd(rmsdHistogramRange());


	size_t numGrains = atomClusters()->size() == 0 ? 0 : (*std::max_element(atomClusters()->constDataInt64(), atomClusters()->constDataInt64() + atomClusters()->size()));
	//output.setStatus(PipelineStatus(PipelineStatus::Success, GrainSegmentationModifier::tr("Found %1 grains").arg(numGrains)));
	//return output;
}

#if 0
/******************************************************************************
* Performs grain tracking over the whole simulation trajectory.
******************************************************************************/
bool GrainSegmentationModifier::trackGrains(TaskManager& taskManager, ModifierApplication* modApp)
{
	Promise<> task = Promise<>::createSynchronous(&taskManager, true, true);

	// Determine the time interval over which grain tracking is performed.
	TimeInterval interval = dataset()->animationSettings()->animationInterval();
	if(interval.duration() <= 0)
		throwException(tr("Loaded simulation sequence consists only of a single frame. Grain tracking is not possible."));

	// Generate list of animation times which will be processed.
	std::vector<TimePoint> frames;
	for(TimePoint frameTime = interval.start(); frameTime <= interval.end(); frameTime += dataset()->animationSettings()->ticksPerFrame()) {
		frames.push_back(frameTime);
	}

	// Stage I: Perform grain segmentation on each frame:
	task.setProgressMaximum(frames.size());
	task.setProgressValue(0);

	// This is where we'll store the segmentation results for each grain.
	std::vector<std::shared_ptr<GrainSegmentationEngine>> grainSegmentationResults;

	// Step through the frames one by one.
	for(TimePoint frameTime : frames) {
		task.setProgressText(tr("Evaluating input pipeline (frame %1 of %2)").arg(task.progressValue()+1).arg(task.progressMaximum()));

		// Request frame data from input pipeline.
		SharedFuture<PipelineFlowState> stateFuture = modApp->evaluateInput(frameTime);
		if(!taskManager.waitForTask(stateFuture))
			return false;
		const PipelineFlowState& state = stateFuture.result();

		// Run regular grain segmentation algorithm for one frame.
		std::shared_ptr<GrainSegmentationEngine> grainSegmentationEngine = createSegmentationEngine(frameTime, modApp, state);
		if(!taskManager.waitForTask(taskManager.runTaskAsync(grainSegmentationEngine->task())))
			return false;
		// Release working data structures that are no longer needed.
		grainSegmentationEngine->cleanup();

		// Stash away the per-frame results for later use in the grain tracking stage below.
		grainSegmentationResults.push_back(std::move(grainSegmentationEngine));

		// Update progress display.
		task.setProgressValue(task.progressValue() + 1);
		if(task.isCanceled()) 
			return false;
	}

	task.setProgressText(tr("Tracking grains"));

	// Stage II: Spawn the grain tracking engine.
	// Pass the per-frame outputs of the grain segmentation step.
	auto grainTrackingEngine = std::make_shared<GrainTrackingEngine>(std::move(grainSegmentationResults));
	if(!taskManager.waitForTask(taskManager.runTaskAsync(grainTrackingEngine->task())))
		return false;

	// Release working data structures that are no longer needed.
	grainTrackingEngine->cleanup();
	
	return true;
}
#endif

}	// End of namespace
}	// End of namespace
}	// End of namespace
