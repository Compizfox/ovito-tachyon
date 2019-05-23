///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2019) Alexander Stukowski
//  Copyright (2019) Henrik Andersen Sveinsson
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

#include <plugins/particles/Particles.h>
#include <core/app/Application.h>
#include <core/dataset/DataSet.h>
#include <plugins/particles/util/CutoffNeighborFinder.h>
#include <plugins/stdobj/simcell/SimulationCellObject.h>
#include <plugins/stdobj/series/DataSeriesObject.h>
#include <core/dataset/pipeline/ModifierApplication.h>
#include <core/utilities/units/UnitsManager.h>
#include <core/utilities/concurrent/ParallelFor.h>
#include "ChillPlusModifier.h"

namespace Ovito { namespace Particles { OVITO_BEGIN_INLINE_NAMESPACE(Modifiers) OVITO_BEGIN_INLINE_NAMESPACE(Analysis)

IMPLEMENT_OVITO_CLASS(ChillPlusModifier);
DEFINE_PROPERTY_FIELD(ChillPlusModifier, cutoff);
SET_PROPERTY_FIELD_LABEL(ChillPlusModifier, cutoff, "Cutoff radius");
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(ChillPlusModifier, cutoff, WorldParameterUnit, 0);

/******************************************************************************
* Constructs the modifier object.
******************************************************************************/
ChillPlusModifier::ChillPlusModifier(DataSet* dataset) : StructureIdentificationModifier(dataset),
    _cutoff(3.5)
{
    createStructureType(CUBIC_ICE, ParticleType::PredefinedStructureType::CUBIC_ICE);
	createStructureType(HEXAGONAL_ICE, ParticleType::PredefinedStructureType::HEXAGONAL_ICE);
	createStructureType(INTERFACIAL_ICE, ParticleType::PredefinedStructureType::INTERFACIAL_ICE);
	createStructureType(HYDRATE, ParticleType::PredefinedStructureType::HYDRATE);
    createStructureType(INTERFACIAL_HYDRATE, ParticleType::PredefinedStructureType::INTERFACIAL_HYDRATE);
    createStructureType(OTHER, ParticleType::PredefinedStructureType::OTHER);
}



/******************************************************************************
* Creates and initializes a computation engine that will compute the 
* modifier's results.
******************************************************************************/
Future<AsynchronousModifier::ComputeEnginePtr> ChillPlusModifier::createEngine(TimePoint time, ModifierApplication* modApp, const PipelineFlowState& input)
{
	if(structureTypes().size() != NUM_STRUCTURE_TYPES)
		throwException(tr("The number of structure types has changed. Please remove this modifier from the pipeline and insert it again."));

	// Get modifier input.
	const ParticlesObject* particles = input.expectObject<ParticlesObject>();
	const PropertyObject* posProperty = particles->expectProperty(ParticlesObject::PositionProperty);
	const SimulationCellObject* simCell = input.expectObject<SimulationCellObject>();
	if(simCell->is2D())
		throwException(tr("Chill+ modifier does not support 2d simulation cells."));

	// Get particle selection.
	ConstPropertyPtr selectionProperty;
	if(onlySelectedParticles())
		selectionProperty = particles->expectProperty(ParticlesObject::SelectionProperty)->storage();

	// Create engine object. Pass all relevant modifier parameters to the engine as well as the input data.
	return std::make_shared<ChillPlusEngine>(particles, posProperty->storage(), simCell->data(), getTypesToIdentify(NUM_STRUCTURE_TYPES), std::move(selectionProperty), cutoff());
}

/******************************************************************************
* Performs the actual analysis. This method is executed in a worker thread.
******************************************************************************/
void ChillPlusModifier::ChillPlusEngine::perform()
{
	task()->setProgressText(tr("Computing q_lm values in Chill+ analysis"));

	// Prepare the neighbor list.
	CutoffNeighborFinder neighborListBuilder;
	if(!neighborListBuilder.prepare(cutoff(), *positions(), cell(), selection().get(), task().get()))
		return;

	// Create output storage.
	PropertyStorage& output = *structures();

    // Find all relevant q_lm
    // create matrix of q_lm 
    size_t particleCount = positions()->size();
	task()->setProgressValue(0);
	task()->setProgressMaximum(particleCount);

    q_values = boost::numeric::ublas::matrix<std::complex<float>>(particleCount, 7);

    task()->setProgressText(tr("Computing c_ij values i Chill+"));
    // Parallel calculation loop:
	std::mutex mutex;
	parallelFor(particleCount, *task(), [this, &neighborListBuilder, &output](size_t index) {
            int coordination = 0;
            for (int m = -3; m<=3; m++) 
            {
                q_values(index, m+3) = compute_q_lm(neighborListBuilder, index, 3, m);
            }
    });

    // For each particle, count the bonds and determine structure
    parallelFor(particleCount, *task(), [this, &neighborListBuilder, &output](size_t index) {
        output.setInt(index, determineStructure(neighborListBuilder, index, typesToIdentify()));
    });
}


std::complex<float> ChillPlusModifier::ChillPlusEngine::compute_q_lm(CutoffNeighborFinder& neighFinder, size_t particleIndex, int l, int m) 
{
    
    std::complex<float> q = 0;
    for(CutoffNeighborFinder::Query neighQuery(neighFinder, particleIndex); !neighQuery.atEnd(); neighQuery.next()) {
        Vector3 delta = neighQuery.delta();
        std::pair<float, float> angles = polar_asimuthal(delta);
        q += boost::math::spherical_harmonic(l, m, angles.first, angles.second);
    }
    return q;
}

/******************************************************************************
* Determines the structure of each atom based on the number of eclipsed and 
* staggered bonds.
******************************************************************************/
ChillPlusModifier::StructureType ChillPlusModifier::ChillPlusEngine::determineStructure(CutoffNeighborFinder& neighFinder, size_t particleIndex, const QVector<bool>& typesToIdentify)
{
    int num_eclipsed = 0;
    int num_staggered = 0;
    int coordination = 0;
    for(CutoffNeighborFinder::Query neighQuery(neighFinder, particleIndex); !neighQuery.atEnd(); neighQuery.next()) {
        // Compute c(i,j)
        std::complex<float> c1 = 0;
        std::complex<float> c2 = 0;
        std::complex<float> c3 = 0;
        std::complex<float> q_i = 0;
        std::complex<float> q_j = 0;
        for (int m = -3; m<=3; m++)
        {
            q_i = q_values(particleIndex, m+3);
            q_j = q_values(neighQuery.current(), m+3);
            c1 += q_i*std::conj(q_j);
            c2 += q_i*std::conj(q_i);
            c3 += q_j*std::conj(q_j);
        }
        std::complex<float> c_ij = c1/(std::sqrt(c2)*std::sqrt(c3));
        if (std::real(c_ij) > -0.35 && std::real(c_ij) < 0.25)
        {
            num_eclipsed ++;
        }
        if (std::real(c_ij)<-0.8)
        {
            num_staggered ++;
        }
        coordination ++;
        
    }           

    if (coordination == 4)
    {
        if (num_eclipsed == 4) 
        {
            return HYDRATE;
        }
        else if (num_eclipsed == 3)
        {
            return INTERFACIAL_HYDRATE;
        }
        else if (num_staggered == 4)
        {
            return CUBIC_ICE;
        }
        else if (num_staggered == 3 && num_eclipsed == 1)
        {
            return HEXAGONAL_ICE;
        }
        else if (num_staggered == 3 && num_eclipsed == 0)
        {
            return INTERFACIAL_ICE;
        }
        else if (num_staggered == 2)
        {
            return INTERFACIAL_ICE;
        }
    }
    else {return OTHER;}
}



/******************************************************************************
* Injects the computed results of the engine into the data pipeline.
******************************************************************************/
void ChillPlusModifier::ChillPlusEngine::emitResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	StructureIdentificationEngine::emitResults(time, modApp, state);

	// Also output structure type counts, which have been computed by the base class.
	state.addAttribute(QStringLiteral("ChillPlus.counts.OTHER"), QVariant::fromValue(getTypeCount(OTHER)), modApp);
	state.addAttribute(QStringLiteral("ChillPlus.counts.CUBIC_ICE"), QVariant::fromValue(getTypeCount(CUBIC_ICE)), modApp);
	state.addAttribute(QStringLiteral("ChillPlus.counts.HEXAGONAL_ICE"), QVariant::fromValue(getTypeCount(HEXAGONAL_ICE)), modApp);
	state.addAttribute(QStringLiteral("ChillPlus.counts.INTERFACIAL_ICE"), QVariant::fromValue(getTypeCount(INTERFACIAL_ICE)), modApp);
	state.addAttribute(QStringLiteral("ChillPlus.counts.HYDRATE"), QVariant::fromValue(getTypeCount(HYDRATE)), modApp);
    state.addAttribute(QStringLiteral("ChillPlus.counts.INTERFACIAL_HYDRATE"), QVariant::fromValue(getTypeCount(INTERFACIAL_HYDRATE)), modApp);
}

std::pair<float, float> ChillPlusModifier::ChillPlusEngine::polar_asimuthal(Vector3 delta)
{
  float asimuthal = std::atan2(delta.y(), delta.x());
  float xy_distance = std::sqrt(delta.x()*delta.x()+delta.y()*delta.y());
  float polar = std::atan2(xy_distance, delta.z());
  return std::pair<float, float>(polar, asimuthal);
}

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace
}	// End of namespace
