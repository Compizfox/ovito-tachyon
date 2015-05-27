///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2015) Alexander Stukowski
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

#ifndef __OVITO_DISLOCATION_ANALYSIS_ENGINE_H
#define __OVITO_DISLOCATION_ANALYSIS_ENGINE_H

#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <plugins/particles/modifier/AsynchronousParticleModifier.h>
#include <plugins/crystalanalysis/util/DelaunayTessellation.h>
#include "StructureAnalysis.h"
#include "ElasticMapping.h"
#include "InterfaceMesh.h"

namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

/*
 * Computation engine of the DislocationAnalysisModifier, which performs the actual dislocation analysis.
 */
class DislocationAnalysisEngine : public AsynchronousParticleModifier::ComputeEngine
{
public:

	/// Constructor.
	DislocationAnalysisEngine(const TimeInterval& validityInterval, ParticleProperty* positions, const SimulationCell& simCell);

	/// Computes the modifier's results and stores them in this object for later retrieval.
	virtual void perform() override;

	/// Returns the generated defect mesh.
	HalfEdgeMesh* defectMesh() { return _defectMesh.data(); }

	/// Indicates whether the entire simulation cell is part of the 'bad' crystal region.
	bool isDefectRegionEverywhere() const { return _isDefectRegionEverywhere; }

private:

	QExplicitlySharedDataPointer<HalfEdgeMesh> _defectMesh;
	StructureAnalysis _structureAnalysis;
	DelaunayTessellation _tessellation;
	ElasticMapping _elasticMapping;
	InterfaceMesh _interfaceMesh;

	bool _isDefectRegionEverywhere;
};

}	// End of namespace
}	// End of namespace
}	// End of namespace

#endif // __OVITO_DISLOCATION_ANALYSIS_ENGINE_H
