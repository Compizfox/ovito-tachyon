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

#pragma once


#include <ovito/particles/Particles.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include <ovito/particles/util/NearestNeighborFinder.h>
#include "PTMAlgorithm.h"
//#include <3rdparty/ptm/ptm_functions.h>


namespace Ovito { namespace Particles {

/**
 * \brief This class is... TODO: complete description
 */
class OVITO_PARTICLES_EXPORT PTMNeighborFinder
{
public:
	PTMNeighborFinder(NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS>& neighQuery,
						ConstPropertyAccess<qlonglong> correspondenceArray,
						ConstPropertyAccess<PTMAlgorithm::StructureType> structuresArray);

	PTMAlgorithm::StructureType structureType;

	/// Contains information about a single neighbor of the central particle.
	struct Neighbor
	{
		Vector3 delta;
		FloatType distanceSq;
		//NeighborListAtom* atom;
		size_t index;
	};

private:
	bool fill_neighbors(NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS>& neighQuery,
						   size_t particleIndex,
						   size_t offset,
						   size_t num,
						   double* delta)
	{
		neighQuery.findNeighbors(particleIndex);
		int numNeighbors = neighQuery.results().size();

		if (numNeighbors < num) {
			return false;
		}

		if (offset == 0) {
			_env.atom_indices[0] = particleIndex;
			_env.points[0][0] = 0;
			_env.points[0][1] = 0;
			_env.points[0][2] = 0;
		}

		for(int i = 0; i < num; i++) {
			int p = _env.correspondences[i + 1 + offset] - 1;
			_env.atom_indices[i + 1 + offset] = neighQuery.results()[p].index;
			_env.points[i + 1 + offset][0] = neighQuery.results()[p].delta.x() + delta[0];
			_env.points[i + 1 + offset][1] = neighQuery.results()[p].delta.y() + delta[1];
			_env.points[i + 1 + offset][2] = neighQuery.results()[p].delta.z() + delta[2];
		}

		return true;
	}

	NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> _neighQuery;
	ConstPropertyAccess<qlonglong> _correspondenceArray;
	ConstPropertyAccess<PTMAlgorithm::StructureType> _structuresArray;
	ptm_atomicenv_t _env;
	std::vector<Neighbor> _results;

public:
	void findNeighbors(size_t particleIndex)
	{
		structureType = (PTMAlgorithm::StructureType)_structuresArray[particleIndex];
		int ptm_type = PTMAlgorithm::ovito_to_ptm_structure_type(structureType);

		_neighQuery.findNeighbors(particleIndex);
		int numNeighbors = _neighQuery.results().size();
		int num_inner = ptm_num_nbrs[ptm_type], num_outer = 0;

		if (ptm_type == PTM_MATCH_NONE) {
			for (int i=0;i<PTM_MAX_INPUT_POINTS;i++) {
				_env.correspondences[i] = i;
			}

			num_inner = numNeighbors;
		}
		else {
			numNeighbors = ptm_num_nbrs[ptm_type];
			ptm_decode_correspondences(ptm_type, _correspondenceArray[particleIndex], _env.correspondences);
		}

		_env.num = numNeighbors + 1;

		if (ptm_type == PTM_MATCH_DCUB || ptm_type == PTM_MATCH_DHEX) {
			num_inner = 4;
			num_outer = 3;
		}
		else if (ptm_type == PTM_MATCH_GRAPHENE) {
			num_inner = 3;
			num_outer = 2;
		}

		fill_neighbors(_neighQuery, particleIndex, 0, num_inner, _env.points[0]);
		if (num_outer) {
			for (int i=0;i<num_inner;i++) {
				fill_neighbors(_neighQuery, _env.atom_indices[1 + i], num_inner + i * num_outer, num_outer, _env.points[i + 1]);
			}
		}

		_results.clear();
		for (int i=0;i<numNeighbors;i++) {
			Neighbor n;
			double* point = _env.points[i + 1];
			n.delta = Vector3(point[0], point[1], point[2]);
			n.distanceSq = n.delta.squaredLength();
			n.index = _env.atom_indices[1 + i];
			//n.atom = atom;
			_results.push_back(n);
		}
	}

	std::vector<Neighbor> results() {return _results;}
};

}	// End of namespace
}	// End of namespace
