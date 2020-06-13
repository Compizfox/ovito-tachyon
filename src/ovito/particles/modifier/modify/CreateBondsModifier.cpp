////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2019 Alexander Stukowski
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
#include <ovito/particles/util/CutoffNeighborFinder.h>
#include <ovito/particles/objects/BondsVis.h>
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include <ovito/core/dataset/DataSet.h>
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include <ovito/core/utilities/concurrent/ParallelFor.h>
#include <ovito/core/utilities/units/UnitsManager.h>
#include "CreateBondsModifier.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(CreateBondsModifier);
DEFINE_PROPERTY_FIELD(CreateBondsModifier, cutoffMode);
DEFINE_PROPERTY_FIELD(CreateBondsModifier, uniformCutoff);
DEFINE_PROPERTY_FIELD(CreateBondsModifier, pairwiseCutoffs);
DEFINE_PROPERTY_FIELD(CreateBondsModifier, minimumCutoff);
DEFINE_PROPERTY_FIELD(CreateBondsModifier, onlyIntraMoleculeBonds);
DEFINE_PROPERTY_FIELD(CreateBondsModifier, autoDisableBondDisplay);
DEFINE_REFERENCE_FIELD(CreateBondsModifier, bondType);
DEFINE_REFERENCE_FIELD(CreateBondsModifier, bondsVis);
SET_PROPERTY_FIELD_LABEL(CreateBondsModifier, cutoffMode, "Cutoff mode");
SET_PROPERTY_FIELD_LABEL(CreateBondsModifier, uniformCutoff, "Cutoff radius");
SET_PROPERTY_FIELD_LABEL(CreateBondsModifier, pairwiseCutoffs, "Pair-wise cutoffs");
SET_PROPERTY_FIELD_LABEL(CreateBondsModifier, minimumCutoff, "Lower cutoff");
SET_PROPERTY_FIELD_LABEL(CreateBondsModifier, onlyIntraMoleculeBonds, "Suppress inter-molecular bonds");
SET_PROPERTY_FIELD_LABEL(CreateBondsModifier, bondType, "Bond type");
SET_PROPERTY_FIELD_LABEL(CreateBondsModifier, bondsVis, "Visual element");
SET_PROPERTY_FIELD_LABEL(CreateBondsModifier, autoDisableBondDisplay, "Auto-disable bond display");
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(CreateBondsModifier, uniformCutoff, WorldParameterUnit, 0);
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(CreateBondsModifier, minimumCutoff, WorldParameterUnit, 0);

/******************************************************************************
* Constructs the modifier object.
******************************************************************************/
CreateBondsModifier::CreateBondsModifier(DataSet* dataset) : AsynchronousModifier(dataset),
	_cutoffMode(UniformCutoff),
	_uniformCutoff(3.2),
	_onlyIntraMoleculeBonds(false),
	_minimumCutoff(0),
	_autoDisableBondDisplay(true)
{
	// Create the bond type that will be assigned to the newly created bonds.
	setBondType(new BondType(dataset));

	// Create the vis element for rendering the bonds generated by the modifier.
	setBondsVis(new BondsVis(dataset));
}

/******************************************************************************
* Is called when a RefTarget referenced by this object has generated an event.
******************************************************************************/
bool CreateBondsModifier::referenceEvent(RefTarget* source, const ReferenceEvent& event)
{
	if(source == bondsVis() && event.type() == ReferenceEvent::TargetEnabledOrDisabled && bondsVis()->isEnabled()) {
		// If the user explicitly re-enables the display of bonds, then the modifier should stop turning it off
		// again in the future.
		setAutoDisableBondDisplay(false);
	}
	return AsynchronousModifier::referenceEvent(source, event);
}

/******************************************************************************
* Asks the modifier whether it can be applied to the given input data.
******************************************************************************/
bool CreateBondsModifier::OOMetaClass::isApplicableTo(const DataCollection& input) const
{
	return input.containsObject<ParticlesObject>();
}

/******************************************************************************
* Sets the cutoff radius for a pair of particle types.
******************************************************************************/
void CreateBondsModifier::setPairwiseCutoff(const QVariant& typeA, const QVariant& typeB, FloatType cutoff)
{
	PairwiseCutoffsList newList = pairwiseCutoffs();
	if(cutoff > 0) {
		newList[qMakePair(typeA, typeB)] = cutoff;
		newList[qMakePair(typeB, typeA)] = cutoff;
	}
	else {
		newList.remove(qMakePair(typeA, typeB));
		newList.remove(qMakePair(typeB, typeA));
	}
	setPairwiseCutoffs(newList);
}

/******************************************************************************
* Returns the pair-wise cutoff radius for a pair of particle types.
******************************************************************************/
FloatType CreateBondsModifier::getPairwiseCutoff(const QVariant& typeA, const QVariant& typeB) const
{
	auto iter = pairwiseCutoffs().find(qMakePair(typeA, typeB));
	if(iter != pairwiseCutoffs().end()) return iter.value();
	iter = pairwiseCutoffs().find(qMakePair(typeB, typeA));
	if(iter != pairwiseCutoffs().end()) return iter.value();
	return 0;
}

/******************************************************************************
* This method is called by the system when the modifier has been inserted
* into a pipeline.
******************************************************************************/
void CreateBondsModifier::initializeModifier(ModifierApplication* modApp)
{
	AsynchronousModifier::initializeModifier(modApp);

	// Adopt the upstream BondsVis object if there is already is one.
	// Also initialize the numeric ID of the type ID to not conflict with any existing bond types.
	int bondTypeId = 1;
	const PipelineFlowState& input = modApp->evaluateInputSynchronous(dataset()->animationSettings()->time());
	if(const ParticlesObject* particles = input.getObject<ParticlesObject>()) {
		if(particles->bonds()) {
			if(BondsVis* bondsVis = particles->bonds()->visElement<BondsVis>()) {
				setBondsVis(bondsVis);
			}
			if(const PropertyObject* bondTypeProperty = particles->bonds()->getProperty(BondsObject::TypeProperty)) {
				bondTypeId = bondTypeProperty->generateUniqueElementTypeId();
			}
		}
	}
	if(bondType() && bondType()->numericId() == 0) {
		bondType()->setNumericId(bondTypeId);
	}
}

/******************************************************************************
* Looks up a particle type in the type list based on the name or the numeric ID.
******************************************************************************/
const ElementType* CreateBondsModifier::lookupParticleType(const PropertyObject* typeProperty, const QVariant& typeSpecification)
{
	if(typeSpecification.type() == QVariant::Int) {
		return typeProperty->elementType(typeSpecification.toInt());
	}
	else {
		const QString& name = typeSpecification.toString();
		for(const ElementType* type : typeProperty->elementTypes())
			if(type->nameOrNumericId() == name)
				return type;
		return nullptr;
	}
}

/******************************************************************************
* Creates and initializes a computation engine that will compute the
* modifier's results.
******************************************************************************/
Future<AsynchronousModifier::EnginePtr> CreateBondsModifier::createEngine(const PipelineEvaluationRequest& request, ModifierApplication* modApp, const PipelineFlowState& input)
{
	// Get modifier input.
	const ParticlesObject* particles = input.expectObject<ParticlesObject>();
	particles->verifyIntegrity();
	const SimulationCellObject* simCell = input.expectObject<SimulationCellObject>();
	const PropertyObject* posProperty = particles->expectProperty(ParticlesObject::PositionProperty);

	// The neighbor list cutoff.
	FloatType maxCutoff = uniformCutoff();

	// Build table of pair-wise cutoff radii.
	const PropertyObject* typeProperty = nullptr;
	std::vector<std::vector<FloatType>> pairCutoffSquaredTable;
	if(cutoffMode() == PairCutoff) {
		typeProperty = particles->expectProperty(ParticlesObject::TypeProperty);
		if(typeProperty) {
			maxCutoff = 0;
			for(auto entry = pairwiseCutoffs().begin(); entry != pairwiseCutoffs().end(); ++entry) {
				FloatType cutoff = entry.value();
				if(cutoff > 0) {
					const ElementType* ptype1 = lookupParticleType(typeProperty, entry.key().first);
					const ElementType* ptype2 = lookupParticleType(typeProperty, entry.key().second);
					if(ptype1 && ptype2 && ptype1->numericId() >= 0 && ptype2->numericId() >= 0) {
						int id1 = ptype1->numericId();
						int id2 = ptype2->numericId();
						if((int)pairCutoffSquaredTable.size() <= std::max(id1, id2)) pairCutoffSquaredTable.resize(std::max(id1, id2) + 1);
						if((int)pairCutoffSquaredTable[id1].size() <= id2) pairCutoffSquaredTable[id1].resize(id2 + 1, FloatType(0));
						if((int)pairCutoffSquaredTable[id2].size() <= id1) pairCutoffSquaredTable[id2].resize(id1 + 1, FloatType(0));
						pairCutoffSquaredTable[id1][id2] = cutoff * cutoff;
						pairCutoffSquaredTable[id2][id1] = cutoff * cutoff;
						if(cutoff > maxCutoff) maxCutoff = cutoff;
					}
				}
			}
			if(maxCutoff <= 0)
				throwException(tr("At least one positive bond cutoff must be set for a valid pair of particle types."));
		}
	}

	// Get molecule IDs.
	ConstPropertyPtr moleculeProperty = onlyIntraMoleculeBonds() ? particles->getPropertyStorage(ParticlesObject::MoleculeProperty) : nullptr;

	// Create engine object. Pass all relevant modifier parameters to the engine as well as the input data.
	return std::make_shared<BondsEngine>(particles, posProperty->storage(),
			typeProperty ? typeProperty->storage() : nullptr, simCell->data(), cutoffMode(),
			maxCutoff, minimumCutoff(), std::move(pairCutoffSquaredTable), std::move(moleculeProperty));
}

/******************************************************************************
* Performs the actual analysis. This method is executed in a worker thread.
******************************************************************************/
void CreateBondsModifier::BondsEngine::perform()
{
	setProgressText(tr("Generating bonds"));

	// Prepare the neighbor list.
	CutoffNeighborFinder neighborFinder;
	if(!neighborFinder.prepare(_maxCutoff, _positions, _simCell, {}, this))
		return;

	FloatType minCutoffSquared = _minCutoff * _minCutoff;

	ConstPropertyAccess<qlonglong> moleculeIDsArray(_moleculeIDs);
	ConstPropertyAccess<int> particleTypesArray(_particleTypes);

	// Generate bonds.
	size_t particleCount = _positions->size();
	setProgressMaximum(particleCount);
	if(!particleTypesArray) {
		for(size_t particleIndex = 0; particleIndex < particleCount; particleIndex++) {
			for(CutoffNeighborFinder::Query neighborQuery(neighborFinder, particleIndex); !neighborQuery.atEnd(); neighborQuery.next()) {
				if(neighborQuery.distanceSquared() < minCutoffSquared)
					continue;
				if(moleculeIDsArray && moleculeIDsArray[particleIndex] != moleculeIDsArray[neighborQuery.current()])
					continue;

				Bond bond = { particleIndex, neighborQuery.current(), neighborQuery.unwrappedPbcShift() };

				// Skip every other bond to create only one bond per particle pair.
				if(!bond.isOdd())
					bonds().push_back(bond);
			}
			// Update progress indicator.
			if(!setProgressValueIntermittent(particleIndex))
				return;
		}
	}
	else {
		for(size_t particleIndex = 0; particleIndex < particleCount; particleIndex++) {
			for(CutoffNeighborFinder::Query neighborQuery(neighborFinder, particleIndex); !neighborQuery.atEnd(); neighborQuery.next()) {
				if(neighborQuery.distanceSquared() < minCutoffSquared)
					continue;
				if(moleculeIDsArray && moleculeIDsArray[particleIndex] != moleculeIDsArray[neighborQuery.current()])
					continue;
				int type1 = particleTypesArray[particleIndex];
				int type2 = particleTypesArray[neighborQuery.current()];
				if(type1 >= 0 && type1 < (int)_pairCutoffsSquared.size() && type2 >= 0 && type2 < (int)_pairCutoffsSquared[type1].size()) {
					if(neighborQuery.distanceSquared() <= _pairCutoffsSquared[type1][type2]) {
						Bond bond = { particleIndex, neighborQuery.current(), neighborQuery.unwrappedPbcShift() };
						// Skip every other bond to create only one bond per particle pair.
						if(!bond.isOdd())
							bonds().push_back(bond);
					}
				}
			}
			// Update progress indicator.
			if(!setProgressValueIntermittent(particleIndex))
				return;
		}
	}
	setProgressValue(particleCount);

	// Release data that is no longer needed.
	_positions.reset();
	_particleTypes.reset();
	_moleculeIDs.reset();
}

/******************************************************************************
* Injects the computed results of the engine into the data pipeline.
******************************************************************************/
void CreateBondsModifier::BondsEngine::applyResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	CreateBondsModifier* modifier = static_object_cast<CreateBondsModifier>(modApp->modifier());
	OVITO_ASSERT(modifier);

	// Add our bonds to the system.
	ParticlesObject* particles = state.expectMutableObject<ParticlesObject>();

	// Bonds have been created for a specific particles ordering. Make sure it's still the same.
	if(_inputFingerprint.hasChanged(particles))
		modApp->throwException(tr("Cached modifier results are obsolete, because the number or the storage order of input particles has changed."));

	particles->addBonds(bonds(), modifier->bondsVis(), {}, modifier->bondType());

	size_t bondsCount = bonds().size();
	state.addAttribute(QStringLiteral("CreateBonds.num_bonds"), QVariant::fromValue(bondsCount), modApp);

	// If the number of bonds is unusually high, we better turn off bonds display to prevent the program from freezing.
	if(bondsCount > 1000000 && modifier->autoDisableBondDisplay() && modifier->bondsVis() && Application::instance()->executionContext() == Application::ExecutionContext::Interactive) {
		modifier->bondsVis()->setEnabled(false);
		state.setStatus(PipelineStatus(PipelineStatus::Warning, tr("Created %1 bonds, which is a lot. As a precaution, the display of bonds has been disabled. You can manually enable it again if needed.").arg(bondsCount)));
	}
	else {
		state.setStatus(PipelineStatus(PipelineStatus::Success, tr("Created %1 bonds.").arg(bondsCount)));
	}
}

}	// End of namespace
}	// End of namespace
