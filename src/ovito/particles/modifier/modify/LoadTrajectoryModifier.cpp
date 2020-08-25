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
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/particles/objects/BondsObject.h>
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include <ovito/core/dataset/io/FileSource.h>
#include <ovito/core/dataset/animation/AnimationSettings.h>
#include <ovito/core/dataset/data/AttributeDataObject.h>
#include "LoadTrajectoryModifier.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(LoadTrajectoryModifier);
DEFINE_REFERENCE_FIELD(LoadTrajectoryModifier, trajectorySource);
SET_PROPERTY_FIELD_LABEL(LoadTrajectoryModifier, trajectorySource, "Trajectory source");

/******************************************************************************
* Constructs the modifier object.
******************************************************************************/
LoadTrajectoryModifier::LoadTrajectoryModifier(DataSet* dataset) : Modifier(dataset)
{
	// Create the file source object, which will be responsible for loading
	// and caching the trajectory data.
	OORef<FileSource> fileSource(new FileSource(dataset));

	setTrajectorySource(fileSource);
}

/******************************************************************************
* Asks the modifier whether it can be applied to the given input data.
******************************************************************************/
bool LoadTrajectoryModifier::OOMetaClass::isApplicableTo(const DataCollection& input) const
{
	return input.containsObject<ParticlesObject>();
}

/******************************************************************************
* Determines the time interval over which a computed pipeline state will remain valid.
******************************************************************************/
TimeInterval LoadTrajectoryModifier::validityInterval(const PipelineEvaluationRequest& request, const ModifierApplication* modApp) const
{
	TimeInterval iv = Modifier::validityInterval(request, modApp);

	if(trajectorySource())
		iv.intersect(trajectorySource()->validityInterval(request));

	return iv;
}

/******************************************************************************
* Modifies the input data.
******************************************************************************/
Future<PipelineFlowState> LoadTrajectoryModifier::evaluate(const PipelineEvaluationRequest& request, ModifierApplication* modApp, const PipelineFlowState& input)
{
	OVITO_ASSERT(input);

	// Get the trajectory data source.
	if(!trajectorySource())
		throwException(tr("No trajectory data source has been set."));

	// Obtain the trajectory frame from the secondary pipeline.
	SharedFuture<PipelineFlowState> trajStateFuture = trajectorySource()->evaluate(request);

	// Wait for the data to become available.
	return trajStateFuture.then(modApp->executor(), [state = input, modApp](const PipelineFlowState& trajState) mutable {

		if(LoadTrajectoryModifier* trajModifier = dynamic_object_cast<LoadTrajectoryModifier>(modApp->modifier())) {
			// Make sure the obtained configuration is valid and ready to use.
			if(trajState.status().type() == PipelineStatus::Error) {
				if(FileSource* fileSource = dynamic_object_cast<FileSource>(trajModifier->trajectorySource())) {
					if(fileSource->sourceUrls().empty())
						modApp->throwException(tr("Please pick a trajectory file."));
				}
				state.setStatus(trajState.status());
			}
			else {
				trajModifier->applyTrajectoryState(state, trajState);

				// Invalidate the synchronous state cache of the modifier application.
				// This is needed to force the pipeline system to call our evaluateSynchronous() method
				// again next time the system request a synchronous state from the pipeline.
				modApp->pipelineCache().invalidateSynchronousState();
			}
		}

		return std::move(state);
	});
}

/******************************************************************************
* Modifies the input data synchronously.
******************************************************************************/
void LoadTrajectoryModifier::evaluateSynchronous(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	if(trajectorySource()) {
		const PipelineFlowState& trajState = trajectorySource()->evaluateSynchronous(time);
		applyTrajectoryState(state, trajState);
	}
}

/******************************************************************************
* Transfers the particle positions from the trajectory frame to the current 
* pipeline input state.
******************************************************************************/
void LoadTrajectoryModifier::applyTrajectoryState(PipelineFlowState& state, const PipelineFlowState& trajState)
{
	if(!trajState)
		throwException(tr("Data source has not been specified yet or is empty. Please pick a trajectory file."));

	// Merge validity intervals of topology and trajectory datasets.
	state.intersectStateValidity(trajState.stateValidity());

	// Get the current particle positions.
	const ParticlesObject* trajectoryParticles = trajState.getObject<ParticlesObject>();
	if(!trajectoryParticles)
		throwException(tr("Trajectory dataset does not contain any particle dataset."));
	trajectoryParticles->verifyIntegrity();

	// Get the topology particle dataset.
	ParticlesObject* particles = state.expectMutableObject<ParticlesObject>();
	particles->verifyIntegrity();

	if(ConstPropertyAccess<Point3> trajectoryPosProperty = trajectoryParticles->getProperty(ParticlesObject::PositionProperty)) {

		// Build particle-to-particle index map.
		std::vector<size_t> indexToIndexMap(particles->elementCount());
		ConstPropertyAccess<qlonglong> identifierProperty = particles->getProperty(ParticlesObject::IdentifierProperty);
		ConstPropertyAccess<qlonglong> trajIdentifierProperty = trajectoryParticles->getProperty(ParticlesObject::IdentifierProperty);
		if(identifierProperty && trajIdentifierProperty) {

			// Build map of particle identifiers in trajectory dataset.
			std::map<qlonglong, size_t> refMap;
			size_t index = 0;
			for(qlonglong id : trajIdentifierProperty) {
				if(refMap.insert(std::make_pair(id, index++)).second == false)
					throwException(tr("Particles with duplicate identifiers detected in trajectory dataset."));
			}

			// Check for duplicate identifiers in topology dataset.
			std::vector<size_t> idSet(identifierProperty.cbegin(), identifierProperty.cend());
			boost::sort(idSet);
			if(boost::adjacent_find(idSet) != idSet.cend())
				throwException(tr("Particles with duplicate identifiers detected in topology dataset."));

			// Build mapping of particle indices from the topology dataset to the corresponding indices in the trajectory dataset.
			const qlonglong* id = identifierProperty.cbegin();
			for(auto& mappedIndex : indexToIndexMap) {
				auto iter = refMap.find(*id);
				if(iter == refMap.end())
					throwException(tr("Particle id %1 from topology dataset not found in trajectory dataset.").arg(*id));
				mappedIndex = iter->second;
				refMap.erase(iter);
				++id;
			}

			// Check if the trajectory dataset contains excess particles that are not present in the topology dataset yet.
			if(!refMap.empty()) {
				// Insert the new particles after the existing particles in the topology dataset.
				particles->setElementCount(particles->elementCount() + refMap.size());
				indexToIndexMap.reserve(indexToIndexMap.size() + refMap.size());

				// Extend index mapping and particle identifier property.
				PropertyAccess<qlonglong> identifierProperty = particles->expectMutableProperty(ParticlesObject::IdentifierProperty);
				qlonglong* id = identifierProperty.begin() + indexToIndexMap.size();
				for(const auto& entry : refMap) {
					*id++ = entry.first;
					indexToIndexMap.push_back(entry.second);
				}
				OVITO_ASSERT(id == identifierProperty.end());
				OVITO_ASSERT(indexToIndexMap.size() == particles->elementCount());
			}
		}
		else {
			// Topology dataset and trajectory data must contain the same number of particles.
			if(trajectoryPosProperty.size() != particles->elementCount()) {
				throwException(tr("Cannot apply trajectories to current particle dataset. Numbers of particles in the trajectory file and in the topology file do not match."));
			}

			// When particle identifiers are not available, use trivial 1-to-1 mapping.
			std::iota(indexToIndexMap.begin(), indexToIndexMap.end(), size_t(0));
		}

		// Transfer particle properties from the trajectory file.
		for(const PropertyObject* property : trajectoryParticles->properties()) {
			if(property->type() == ParticlesObject::IdentifierProperty)
				continue;

			// Get or create the output particle property.
			PropertyObject* outputProperty;
			bool replacingProperty;
			if(property->type() != ParticlesObject::UserProperty) {
				replacingProperty = (particles->getProperty(property->type()) != nullptr);
				outputProperty = particles->createProperty(property->type(), true);
				if(outputProperty->dataType() != property->dataType()
					|| outputProperty->componentCount() != property->componentCount())
					continue; // Types of source property and output property are not compatible.
			}
			else {
				replacingProperty = (particles->getProperty(property->name()) != nullptr);
				outputProperty = particles->createProperty(property->name(),
					property->dataType(), property->componentCount(),
					0, true);
			}
			OVITO_ASSERT(outputProperty->stride() == property->stride());

			// Copy and reorder property data.
			property->mappedCopyTo(outputProperty, indexToIndexMap);

			// Transfer the visual element(s) unless the property already existed in the topology dataset.
			if(!replacingProperty) {
				outputProperty->setVisElements(property->visElements());
			}
		}

		// Transfer box geometry.
		const SimulationCellObject* topologyCell = state.getObject<SimulationCellObject>();
		const SimulationCellObject* trajectoryCell = trajState.getObject<SimulationCellObject>();
		if(topologyCell && trajectoryCell) {
			SimulationCellObject* outputCell = state.makeMutable(topologyCell);
			outputCell->setCellMatrix(trajectoryCell->cellMatrix());
			const AffineTransformation& simCell = trajectoryCell->cellMatrix();

			// Trajectories of atoms may cross periodic boundaries and if atomic positions are
			// stored in wrapped coordinates, then it becomes necessary to fix bonds using the minimum image convention.
			std::array<bool, 3> pbc = topologyCell->pbcFlags();
			if((pbc[0] || pbc[1] || pbc[2]) && particles->bonds() && std::abs(simCell.determinant()) > FLOATTYPE_EPSILON) {
				ConstPropertyAccess<Point3> outputPosProperty = particles->expectProperty(ParticlesObject::PositionProperty);
				AffineTransformation inverseCellMatrix = simCell.inverse();

				BondsObject* bonds = particles->makeBondsMutable();
				if(ConstPropertyAccess<ParticleIndexPair> topologyProperty = bonds->getProperty(BondsObject::TopologyProperty)) {
					PropertyAccess<Vector3I> periodicImageProperty = bonds->createProperty(BondsObject::PeriodicImageProperty, true);

					// Wrap bonds crossing a periodic boundary by resetting their PBC shift vectors.
					SimulationCell cell = trajectoryCell->data();
					for(size_t bondIndex = 0; bondIndex < topologyProperty.size(); bondIndex++) {
						size_t particleIndex1 = topologyProperty[bondIndex][0];
						size_t particleIndex2 = topologyProperty[bondIndex][1];
						if(particleIndex1 >= outputPosProperty.size() || particleIndex2 >= outputPosProperty.size())
							continue;
						const Point3& p1 = outputPosProperty[particleIndex1];
						const Point3& p2 = outputPosProperty[particleIndex2];
						Vector3 delta = p1 - p2;
						for(int dim = 0; dim < 3; dim++) {
							if(pbc[dim])
								periodicImageProperty[bondIndex][dim] = std::lround(inverseCellMatrix.prodrow(delta, dim));
						}
					}
				}
			}
		}
	}
	else if(const BondsObject* trajectoryBonds = trajectoryParticles->bonds()) {
		trajectoryBonds->verifyIntegrity();

		// Create a mutable copy of the particles object.
		ParticlesObject* particles = state.expectMutableObject<ParticlesObject>();
		particles->verifyIntegrity();

		// If the trajectory file contains bond topology, completely replace all existing bonds 
		// from the topology dataset with the new set of bonds. 
		if(trajectoryBonds->getProperty(BondsObject::TopologyProperty)) {
			if(OORef<BondsObject> oldBonds = particles->bonds()) {
				// Replace the property arrays, but make sure BondType instances  
				// as well as the visual elements from the topology dataset are preserved.
				particles->makeBondsMutable()->setContent(trajectoryBonds->elementCount(), trajectoryBonds->properties());
			}
			else {
				// We can simply adopt the bonds object from the trajectory dataset as a whole
				// if the topology dataset didn't contain any bonds yet.
				particles->setBonds(trajectoryBonds);
			}

			// Compute the PBC shift vectors of the bonds based on current particle positions.
			if(const SimulationCellObject* simCellObj = state.getObject<SimulationCellObject>()) {
				if(simCellObj->pbcX() || simCellObj->pbcY() || simCellObj->pbcZ()) {
					particles->makeBondsMutable()->generatePeriodicImageProperty(particles, simCellObj);
				}
			}
		}
		else if(particles->bonds()) {
			// If the trajectory dataset doesn't contain the "Topology" bond property, 
			// then add the bond properties to the existing bonds from the topology dataset.
			// This requires that the number of bonds remains constant.
			if(trajectoryBonds->elementCount() != particles->bonds()->elementCount()) {
				throwException(tr("Cannot merge bond properties of trajectory dataset with topology dataset, because numbers of bonds in the two datasets do not match."));
			}

			if(!trajectoryBonds->properties().empty()) {
				BondsObject* bonds = particles->makeBondsMutable();

				// Add the properties to the existing bonds, overwriting existing values if necessary.
				for(const PropertyObject* newProperty : trajectoryBonds->properties()) {
					const PropertyObject* existingPropertyObj = (newProperty->type() != 0) ? bonds->getProperty(newProperty->type()) : bonds->getProperty(newProperty->name());
					if(existingPropertyObj) {
						bonds->makeMutable(existingPropertyObj)->setStorage(newProperty->storage());
					}
					else {
						bonds->addProperty(newProperty);
					}
				}
			}
		}
		else {
			throwException(tr("Neither the trajectory nor the topology dataset contain bond connectivity information."));
		}
	}
	else {
		throwException(tr("Trajectory dataset does not contain any particle positions."));
	}

	// Merge global attributes of topology and trajectory datasets.
	// If there is a naming collision, attributes from the trajectory dataset override those from the topology dataset.
	for(const DataObject* obj : trajState.data()->objects()) {
		if(const AttributeDataObject* attribute = dynamic_object_cast<AttributeDataObject>(obj)) {
			const AttributeDataObject* existingAttribute = nullptr;
			for(const DataObject* obj2 : state.data()->objects()) {
				if(const AttributeDataObject* attribute2 = dynamic_object_cast<AttributeDataObject>(obj2)) {
					if(attribute2->identifier() == attribute->identifier()) {
						existingAttribute = attribute2;
						break;
					}
				}
			}
			if(existingAttribute)
				state.mutableData()->replaceObject(existingAttribute, attribute);
			else
				state.addObject(attribute);
		}
	}
}

/******************************************************************************
* Is called when a RefTarget referenced by this object has generated an event.
******************************************************************************/
bool LoadTrajectoryModifier::referenceEvent(RefTarget* source, const ReferenceEvent& event)
{
	if(event.type() == ReferenceEvent::AnimationFramesChanged && source == trajectorySource()) {
		// Propagate animation interval events from the trajectory source.
		return true;
	}
	return Modifier::referenceEvent(source, event);
}

/******************************************************************************
* Gets called when the data object of the node has been replaced.
******************************************************************************/
void LoadTrajectoryModifier::referenceReplaced(const PropertyFieldDescriptor& field, RefTarget* oldTarget, RefTarget* newTarget)
{
	if(field == PROPERTY_FIELD(trajectorySource) && !isBeingLoaded()) {
		// The animation length might have changed when the trajectory source has been replaced.
		notifyDependents(ReferenceEvent::AnimationFramesChanged);
	}
	Modifier::referenceReplaced(field, oldTarget, newTarget);
}

}	// End of namespace
}	// End of namespace
