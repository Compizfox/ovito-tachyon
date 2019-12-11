////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2019 Alexander Stukowski
//  Copyright 2019 Peter Mahler Larsen
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

#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <plugins/particles/objects/BondsVis.h>
#include <plugins/particles/objects/ParticlesObject.h>
#include <plugins/stdobj/simcell/SimulationCellObject.h>
#include <plugins/stdobj/properties/PropertyStorage.h>
#include <plugins/stdobj/series/DataSeriesObject.h>
#include <core/utilities/units/UnitsManager.h>
#include <core/dataset/pipeline/ModifierApplication.h>
#include <core/dataset/DataSet.h>
#include "GrainSegmentationModifier.h"
#include "GrainSegmentationEngine.h"

#include <ptm/ptm_functions.h>

namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

IMPLEMENT_OVITO_CLASS(GrainSegmentationModifier);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, rmsdCutoff);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, algorithmType);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, mergingThreshold);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, minGrainAtomCount);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, orphanAdoption);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, onlySelectedParticles);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, outputBonds);
DEFINE_REFERENCE_FIELD(GrainSegmentationModifier, bondsVis);
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, rmsdCutoff, "RMSD cutoff");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, algorithmType, "Linkage type");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, mergingThreshold, "Log merge threshold");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, minGrainAtomCount, "Minimum grain size (# of atoms)");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, orphanAdoption, "Adopt orphan atoms");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, onlySelectedParticles, "Use only selected particles");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, outputBonds, "Output bonds");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, bondsVis, "Bonds display");
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(GrainSegmentationModifier, rmsdCutoff, FloatParameterUnit, 0);
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(GrainSegmentationModifier, minGrainAtomCount, IntegerParameterUnit, 1);

/******************************************************************************
* Constructs the modifier object.
******************************************************************************/
GrainSegmentationModifier::GrainSegmentationModifier(DataSet* dataset) : StructureIdentificationModifier(dataset),
		_rmsdCutoff(0.1),
		_algorithmType(0),
		_minGrainAtomCount(100),
		_onlySelectedParticles(false),
		_mergingThreshold(0.0),
		_orphanAdoption(true),
		_outputBonds(false)
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

	// Create the visual element for the bonds.
	setBondsVis(new BondsVis(dataset));
}

/******************************************************************************
* Is called when the value of a property of this object has changed.
******************************************************************************/
void GrainSegmentationModifier::propertyChanged(const PropertyFieldDescriptor& field)
{
	if(field == PROPERTY_FIELD(mergingThreshold)) {
		// Immediately update viewports when threshold parameter is changed by the user.
		notifyDependents(ReferenceEvent::PreliminaryStateAvailable);
	}
	StructureIdentificationModifier::propertyChanged(field);
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
			rmsdCutoff(), algorithmType(), mergingThreshold(), minGrainAtomCount(), orphanAdoption(), outputBonds());
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

	// Output the edges of the neighbor graph.
	if(_outputBonds && modifier->outputBonds()) {
		particles->addBonds(latticeNeighborBonds(), modifier->bondsVis(), { neighborDisorientationAngles() });
	}

	// Output RMSD histogram.
	DataSeriesObject* seriesObj = state.createObject<DataSeriesObject>(QStringLiteral("grains-rmsd"), modApp, DataSeriesObject::Line, GrainSegmentationModifier::tr("RMSD distribution"), rmsdHistogram());
	seriesObj->setAxisLabelX(GrainSegmentationModifier::tr("RMSD"));
	seriesObj->setIntervalStart(0);
	seriesObj->setIntervalEnd(rmsdHistogramRange());

	if (_numClusters > 1) {
		// Output a data series object with the scatter points.
		DataSeriesObject* mergeSeriesObj = state.createObject<DataSeriesObject>(QStringLiteral("grains-merge"), modApp, DataSeriesObject::Scatter, GrainSegmentationModifier::tr("Merge size vs. Merge distance"));
		mergeSeriesObj->createProperty(std::move(mergeDistance()))->setName("Log merge distance");
		mergeSeriesObj->createProperty(std::move(mergeSize()))->setName("Merge size");
	}

	size_t numGrains = atomClusters()->size() == 0 ? 0 : (*std::max_element(atomClusters()->constDataInt64(), atomClusters()->constDataInt64() + atomClusters()->size()));
	state.setStatus(PipelineStatus(PipelineStatus::Success, GrainSegmentationModifier::tr("Found %1 grains").arg(numGrains)));
}

}	// End of namespace
}	// End of namespace
}	// End of namespace
