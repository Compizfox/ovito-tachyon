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

#pragma once


#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/particles/modifier/analysis/StructureIdentificationModifier.h>
#include <ovito/particles/objects/BondsVis.h>

namespace Ovito { namespace CrystalAnalysis {

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

	/// This method indicates whether cached computation results of the modifier should be discarded whenever
	/// a parameter of the modifier changes.
	virtual bool discardResultsOnModifierChange(const PropertyFieldEvent& event) const override {
		// Avoid a recomputation from scratch if just the threshold value is changed.
		if(event.field() == &PROPERTY_FIELD(mergingThreshold)) return false;
		return StructureIdentificationModifier::discardResultsOnModifierChange(event);
	}

protected:

	/// Is called when the value of a property of this object has changed.
	virtual void propertyChanged(const PropertyFieldDescriptor& field) override;

	/// Creates a computation engine that finds the grains in a single frame.
	std::shared_ptr<GrainSegmentationEngine> createSegmentationEngine(TimePoint time, ModifierApplication* modApp, const PipelineFlowState& input);

	/// Creates a computation engine that will compute the modifier's results.
	virtual Future<ComputeEnginePtr> createEngine(TimePoint time, ModifierApplication* modApp, const PipelineFlowState& input) override;

private:

	/// The RMSD cutoff for the PTM algorithm.
	DECLARE_MODIFIABLE_PROPERTY_FIELD_FLAGS(FloatType, rmsdCutoff, setRmsdCutoff, PROPERTY_FIELD_MEMORIZE);

	/// The merging algorithm to use.
	DECLARE_MODIFIABLE_PROPERTY_FIELD(bool, algorithmType, setAlgorithmType);

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
