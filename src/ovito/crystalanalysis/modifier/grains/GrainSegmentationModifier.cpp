////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2020 Alexander Stukowski
//  Copyright 2020 Peter Mahler Larsen
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
#include <ovito/stdobj/table/DataTable.h>
#include <ovito/core/utilities/units/UnitsManager.h>
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include <ovito/core/dataset/DataSet.h>
#include "GrainSegmentationModifier.h"
#include "GrainSegmentationEngine.h"

#include <ptm/ptm_functions.h>

namespace Ovito { namespace CrystalAnalysis {

IMPLEMENT_OVITO_CLASS(GrainSegmentationModifier);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, mergeAlgorithm);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, handleCoherentInterfaces);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, mergingThreshold);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, minGrainAtomCount);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, orphanAdoption);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, outputBonds);
DEFINE_PROPERTY_FIELD(GrainSegmentationModifier, colorParticlesByGrain);
DEFINE_REFERENCE_FIELD(GrainSegmentationModifier, bondsVis);
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, mergeAlgorithm, "Algorithm");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, handleCoherentInterfaces, "Handle coherent interfaces/stacking faults");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, mergingThreshold, "Merge threshold");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, minGrainAtomCount, "Minimum grain size (# of atoms)");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, orphanAdoption, "Adopt orphan atoms");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, outputBonds, "Output bonds");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, colorParticlesByGrain, "Color particles by grain");
SET_PROPERTY_FIELD_LABEL(GrainSegmentationModifier, bondsVis, "Bonds display");
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(GrainSegmentationModifier, minGrainAtomCount, IntegerParameterUnit, 1);

/******************************************************************************
* Constructs the modifier object.
******************************************************************************/
GrainSegmentationModifier::GrainSegmentationModifier(DataSet* dataset) : AsynchronousModifier(dataset),
		_mergeAlgorithm(GraphClusteringAutomatic),
		_minGrainAtomCount(100),
		_mergingThreshold(0.0),
		_orphanAdoption(true),
		_outputBonds(false),
		_colorParticlesByGrain(true),
		_handleCoherentInterfaces(true)
{
	// Create the visual element for the bonds.
	setBondsVis(new BondsVis(dataset));
}

/******************************************************************************
* Asks the modifier whether it can be applied to the given input data.
******************************************************************************/
bool GrainSegmentationModifier::OOMetaClass::isApplicableTo(const DataCollection& input) const
{
	return input.containsObject<ParticlesObject>();
}

/******************************************************************************
* Is called when the value of a property of this object has changed.
******************************************************************************/
void GrainSegmentationModifier::propertyChanged(const PropertyFieldDescriptor& field)
{
	if(field == PROPERTY_FIELD(mergingThreshold) || field == PROPERTY_FIELD(minGrainAtomCount) || field == PROPERTY_FIELD(colorParticlesByGrain) || field == PROPERTY_FIELD(orphanAdoption)) {
		// Immediately update viewports if parameters are changed by the user that don't require a full recalculation.
		notifyDependents(ReferenceEvent::PreliminaryStateAvailable);
	}
	AsynchronousModifier::propertyChanged(field);
}

/******************************************************************************
* Creates a computation engine that will compute the modifier's results.
******************************************************************************/
Future<AsynchronousModifier::EnginePtr> GrainSegmentationModifier::createEngine(const PipelineEvaluationRequest& request, ModifierApplication* modApp, const PipelineFlowState& input)
{
	// Get modifier input.
	const ParticlesObject* particles = input.expectObject<ParticlesObject>();
	particles->verifyIntegrity();
	const PropertyObject* posProperty = particles->expectProperty(ParticlesObject::PositionProperty);

	// Make sure the PTM modifier has been executed first and its output is available.
	const PropertyObject* structureProperty = particles->getProperty(ParticlesObject::StructureTypeProperty);
	if(!structureProperty)
		throwException(tr("Grain segmentation requires Polyhedral Template Matching (PTM) output. Please insert the PTM modifier into the pipeline first."));
	const PropertyObject* orientationProperty = particles->getProperty(ParticlesObject::OrientationProperty);
	if(!orientationProperty)
		throwException(tr("Grain segmentation requires lattice orientation information. Please activate the 'Lattice orientations' option of the PTM modifier."));
	const PropertyObject* correspondenceProperty = particles->expectProperty("Correspondences", PropertyStorage::Int64);

	const SimulationCellObject* simCell = input.expectObject<SimulationCellObject>();
	if(simCell->is2D())
		throwException(tr("The grain segmentation modifier does not support 2d simulation cells."));

	// Initialize PTM library.
	ptm_initialize_global();

	// Create engine object. Pass all relevant modifier parameters to the engine as well as the input data.
	return std::make_shared<GrainSegmentationEngine1>(
			particles,
			posProperty->storage(),
			structureProperty->storage(),
			orientationProperty->storage(),
			correspondenceProperty->storage(),
			simCell->data(),
			mergeAlgorithm(),
			handleCoherentInterfaces(),
			outputBonds());
}

/******************************************************************************
* Injects the computed results of the engine into the data pipeline.
******************************************************************************/
void GrainSegmentationEngine1::applyResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(modApp->modifier());
	OVITO_ASSERT(modifier);

	ParticlesObject* particles = state.expectMutableObject<ParticlesObject>();
	particles->verifyIntegrity();

	if(_inputFingerprint.hasChanged(particles))
		modApp->throwException(GrainSegmentationModifier::tr("Cached modifier results are obsolete, because the number or the storage order of input particles has changed."));

	// Output the edges of the neighbor graph.
	if(_outputBondsToPipeline && modifier->outputBonds()) {

		std::vector<Bond> bonds;
		std::vector<FloatType> disorientations;
		ConstPropertyAccess<Point3> positionsArray(particles->expectProperty(ParticlesObject::PositionProperty));

		for (auto edge: neighborBonds()) {
			if (isCrystallineBond(edge)) {
				Bond bond;
				bond.index1 = edge.a;
				bond.index2 = edge.b;
				disorientations.push_back(edge.disorientation);

				// Determine PBC bond shift using minimum image convention.
				Vector3 delta = positionsArray[bond.index1] - positionsArray[bond.index2];
				for(size_t dim = 0; dim < 3; dim++) {
					if(cell().pbcFlags()[dim])
						bond.pbcShift[dim] = (int)std::floor(cell().inverseMatrix().prodrow(delta, dim) + FloatType(0.5));
					else
						bond.pbcShift[dim] = 0;
				}

				bonds.push_back(bond);
			}
		}

		// Output disorientation angles as a bond property.
		PropertyAccessAndRef<FloatType> neighborDisorientationAngles = std::make_shared<PropertyStorage>(bonds.size(), PropertyStorage::Float, 1, 0, QStringLiteral("Disorientation"), false);
		for (size_t i=0;i<disorientations.size();i++) {
			neighborDisorientationAngles[i] = disorientations[i];
		}

		particles->addBonds(bonds, modifier->bondsVis(), { neighborDisorientationAngles.takeStorage() });
	}

	// Output a data plot with the dendrogram points.
	if(mergeSize() && mergeDistance())
		state.createObject<DataTable>(QStringLiteral("grains-merge"), modApp, DataTable::Scatter, GrainSegmentationModifier::tr("Merge size vs. distance"), mergeSize(), mergeDistance());

	// Output a data plot with the log-log dendrogram points.
	if(logMergeSize() && logMergeDistance())
		state.createObject<DataTable>(QStringLiteral("grains-log"), modApp, DataTable::Scatter, GrainSegmentationModifier::tr("Log distance vs. log merge size"), logMergeDistance(), logMergeSize());

	if(modifier->mergeAlgorithm() == GrainSegmentationModifier::GraphClusteringAutomatic)
		state.addAttribute(QStringLiteral("GrainSegmentation.auto_merge_threshold"), QVariant::fromValue(suggestedMergingThreshold()), modApp);
}

/******************************************************************************
* Injects the computed results of the engine into the data pipeline.
******************************************************************************/
void GrainSegmentationEngine2::applyResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	// Output the results from the 1st algorithm stage.
	_engine1->applyResults(time, modApp, state);

	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(modApp->modifier());
	OVITO_ASSERT(modifier);

	ParticlesObject* particles = state.expectMutableObject<ParticlesObject>();

	// Output per-particle properties.
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

	// Output a data table with the list of grains.
	// The X-column consists of the grain IDs, the Y-column contains the grain sizes. 
	DataTable* grainTable = state.createObject<DataTable>(QStringLiteral("grains"), modApp, DataTable::Scatter, GrainSegmentationModifier::tr("Grain list"), _grainSizes, _grainIds);
	// Add extra columns to the table containing other per-grain data.
	grainTable->createProperty(_grainColors);
	PropertyObject* grainStructureTypesProperty = grainTable->createProperty(_grainStructureTypes);
	grainTable->createProperty(_grainOrientations);

	// Transfer the set of PTM crystal structure types to the structure column of the grain table. 
	const PropertyObject* structureProperty = particles->expectProperty(ParticlesObject::StructureTypeProperty);
	for(const ElementType* type : structureProperty->elementTypes()) {
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
