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

#pragma once


#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <plugins/particles/modifier/analysis/StructureIdentificationModifier.h>
#include <plugins/particles/objects/BondsVis.h>
#include <plugins/particles/objects/ParticlesObject.h>
#include <plugins/stdobj/series/DataSeriesObject.h>

namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

class GrainSegmentationEngine;  // defined in GrainSegmentationEngine.h

/*
 * Decomposes a polycrystalline microstructure into individual grains.
 */
class OVITO_CRYSTALANALYSIS_EXPORT GrainSegmentationModifier : public StructureIdentificationModifier
{
	Q_OBJECT
	OVITO_CLASS(GrainSegmentationModifier)

	Q_CLASSINFO("DisplayName", "Grain segmentation");
	Q_CLASSINFO("ModifierCategory", "Analysis");

public:

	/// Constructor.
	Q_INVOKABLE GrainSegmentationModifier(DataSet* dataset);

	/// Performs grain tracking over the whole simulation trajectory.
	bool trackGrains(TaskManager& taskManager, ModifierApplication* modApp);

protected:

	/// Creates a computation engine that finds the grains in a single frame.
	std::shared_ptr<GrainSegmentationEngine> createSegmentationEngine(TimePoint time, ModifierApplication* modApp, const PipelineFlowState& input);

	/// Creates a computation engine that will compute the modifier's results.
	virtual Future<ComputeEnginePtr> createEngine(TimePoint time, ModifierApplication* modApp, const PipelineFlowState& input) override;

private:

	/// The RMSD cutoff for the PTM algorithm.
	DECLARE_MODIFIABLE_PROPERTY_FIELD_FLAGS(FloatType, rmsdCutoff, setRmsdCutoff, PROPERTY_FIELD_MEMORIZE);

	/// Controls the amount of noise allowed inside a grain.
	DECLARE_MODIFIABLE_PROPERTY_FIELD(FloatType, mergingThreshold, setMergingThreshold);

	/// The minimum number of crystalline atoms per grain.
	DECLARE_MODIFIABLE_PROPERTY_FIELD(int, minGrainAtomCount, setMinGrainAtomCount);

	/// Controls whether to adopt orphan atoms
	DECLARE_MODIFIABLE_PROPERTY_FIELD_FLAGS(bool, orphanAdoption, setOrphanAdoption, PROPERTY_FIELD_MEMORIZE);

	/// Controls whether only selected particles should be taken into account.
	DECLARE_MODIFIABLE_PROPERTY_FIELD(bool, onlySelectedParticles, setOnlySelectedParticles);

	/// The visual element for rendering the bonds created by the modifier.
	DECLARE_MODIFIABLE_REFERENCE_FIELD_FLAGS(BondsVis, bondsVis, setBondsVis, PROPERTY_FIELD_DONT_PROPAGATE_MESSAGES | PROPERTY_FIELD_MEMORIZE);

	/// Controls the output of bonds by the modifier.
	DECLARE_MODIFIABLE_PROPERTY_FIELD(bool, outputBonds, setOutputBonds);
};

}	// End of namespace
}	// End of namespace
}	// End of namespace
