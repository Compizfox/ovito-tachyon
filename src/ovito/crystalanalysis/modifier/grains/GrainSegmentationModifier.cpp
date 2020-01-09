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

#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/particles/objects/BondsVis.h>
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include <ovito/stdobj/properties/PropertyStorage.h>
#include <ovito/stdobj/series/DataSeriesObject.h>
#include <ovito/core/utilities/units/UnitsManager.h>
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include <ovito/core/dataset/DataSet.h>
#include "GrainSegmentationModifier.h"
#include "GrainSegmentationEngine.h"

#include <ptm/ptm_functions.h>

namespace Ovito { namespace CrystalAnalysis {

IMPLEMENT_OVITO_CLASS(GrainSegmentationModifier);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, rmsdCutoff);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, algorithmType);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, mergingThreshold);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, minGrainAtomCount);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, orphanAdoption);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, onlySelectedParticles);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, outputBonds);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, colorParticlesByGrain);
DEFINE_REFERENCE_FIELD(GrainSegmentationModifier, bondsVis);
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, rmsdCutoff, "RMSD cutoff");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, algorithmType, "Linkage type");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, mergingThreshold, "Log merge threshold");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, minGrainAtomCount, "Minimum grain size (# of atoms)");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, orphanAdoption, "Adopt orphan atoms");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, onlySelectedParticles, "Use only selected particles");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, outputBonds, "Output bonds");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, colorParticlesByGrain, "Color particles by grain");
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
		_outputBonds(false),
		_colorParticlesByGrain(true)
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
	if(field == PROPERTY_FIELD(mergingThreshold) || field == PROPERTY_FIELD(minGrainAtomCount)) {
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
			rmsdCutoff(), algorithmType(), minGrainAtomCount(), orphanAdoption(), outputBonds());
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

	// Complete the segmentation by executing the merge steps up to the cutoff set by the user.
	executeMergeSequence(modifier->minGrainAtomCount(), modifier->mergingThreshold());

	ParticlesObject* particles = state.expectMutableObject<ParticlesObject>();

	// Output per-particle properties.
	if(rmsd())
		particles->createProperty(rmsd());
	if(orientations())
		particles->createProperty(orientations());
	if(atomClusters()) {
		particles->createProperty(atomClusters());

		if(modifier->colorParticlesByGrain()) {
			
			// Assign colors to particles according to the grains they belong to.
			ConstPropertyAccess<Color> grainColorsArray(_grainColors);
			PropertyAccess<Color> particleColorsArray = particles->createProperty(ParticlesObject::ColorProperty, false);
			boost::transform(ConstPropertyAccess<qlonglong>(atomClusters()), particleColorsArray.begin(), [&](qlonglong cluster) { 
				if(cluster != 0)
					return grainColorsArray[cluster - 1];
				else
					return Color(0.8, 0.8, 0.8); // Special color for non-crystalline particles not part of any grain.
			});
		}
	}


	// Output the edges of the neighbor graph.
	if(_outputBondsToPipeline && modifier->outputBonds()) {

		// Output disorientation angles as a bond property.
		PropertyAccessAndRef<FloatType> neighborDisorientationAngles = std::make_shared<PropertyStorage>(neighborBonds().size(), PropertyStorage::Float, 1, 0, QStringLiteral("Disorientation"), false);
		// Allocate the bonds array.
		std::vector<Bond> bonds(neighborBonds().size());

		ConstPropertyAccess<Point3> positionsArray(particles->expectProperty(ParticlesObject::PositionProperty));
		for(size_t i = 0; i < bonds.size(); i++) {
			Bond& bond = bonds[i];
			bond.index1 = neighborBonds()[i].a;
			bond.index2 = neighborBonds()[i].b;
			neighborDisorientationAngles[i] = neighborBonds()[i].disorientation;

			// Determine PBC bond shift using minimum image convention.
			Vector3 delta = positionsArray[bond.index1] - positionsArray[bond.index2];
			for(size_t dim = 0; dim < 3; dim++) {
				if(cell().pbcFlags()[dim])
					bond.pbcShift[dim] = (int)std::floor(cell().inverseMatrix().prodrow(delta, dim) + FloatType(0.5));
				else
					bond.pbcShift[dim] = 0;
			}
		}

		particles->addBonds(bonds, modifier->bondsVis(), { neighborDisorientationAngles.takeStorage() });
	}

	// Output RMSD histogram.
	DataSeriesObject* seriesObj = state.createObject<DataSeriesObject>(QStringLiteral("grains-rmsd"), modApp, DataSeriesObject::Line, GrainSegmentationModifier::tr("RMSD distribution"), rmsdHistogram());
	seriesObj->setAxisLabelX(GrainSegmentationModifier::tr("RMSD"));
	seriesObj->setIntervalStart(0);
	seriesObj->setIntervalEnd(rmsdHistogramRange());

	// Output a data series object with the dendrogram points.
	if(mergeSize() && mergeDistance())
		state.createObject<DataSeriesObject>(QStringLiteral("grains-merge"), modApp, DataSeriesObject::Scatter, GrainSegmentationModifier::tr("Merge size vs. Merge distance"), mergeSize(), mergeDistance());

	// Output a data series object with the list of grains.
	// The X-column consists of the grain IDs, the Y-column contains the grain sizes. 
	DataSeriesObject* grainListObj = state.createObject<DataSeriesObject>(QStringLiteral("grains"), modApp, DataSeriesObject::Scatter, GrainSegmentationModifier::tr("Grains"), _grainSizes, _grainIds);
	// Add extra columns to the table containing other per-grain data.
	grainListObj->createProperty(_grainColors);
	PropertyObject* grainStructureTypesProperty = grainListObj->createProperty(_grainStructureTypes);
	grainListObj->createProperty(_grainOrientations);

	// Transfer the set of crystal structure types to the structure column of the grain table. 
	for(const ElementType* type : modifier->structureTypes()) {
		if(type->enabled())
			grainStructureTypesProperty->addElementType(type);
	}

	size_t numGrains = 0;
	if(atomClusters()->size() != 0)
		numGrains = *boost::max_element(ConstPropertyAccess<qlonglong>(atomClusters()));

	state.addAttribute(QStringLiteral("GrainSegmentation.grain_count"), QVariant::fromValue(numGrains), modApp);
	state.setStatus(PipelineStatus(PipelineStatus::Success, GrainSegmentationModifier::tr("Found %1 grains").arg(numGrains)));
}

}	// End of namespace
}	// End of namespace
