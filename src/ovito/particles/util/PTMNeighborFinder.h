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
#include <ovito/particles/modifier/analysis/ptm/PTMAlgorithm.h>

namespace Ovito { namespace Particles {

/**
 * \brief This utility class finds the *k* nearest neighbors of a particle or around some point in space.
 *        *k* is a positive integer.
 *
 * OVITO provides two facilities for finding the neighbors of particles: The CutoffNeighborFinder class, which
 * finds all neighbors within a certain cutoff radius, and the PTMNeighborFinder class, which finds
 * the *k* nearest neighbor of a particle, where *k* is some positive integer. Note that the cutoff-based neighbor finder
 * can return an unknown number of neighbor particles, while the nearest neighbor finder will return exactly
 * the requested number of nearest neighbors (ordered by increasing distance from the central particle).
 * Whether CutoffNeighborFinder or PTMNeighborFinder is the right choice depends on the application.
 *
 * The PTMNeighborFinder class must be initialized by a call to prepare(). This function sorts all input particles
 * in a binary search for fast nearest neighbor queries.
 *
 * After the PTMNeighborFinder has been initialized, one can find the nearest neighbors of some central
 * particle by constructing an instance of the PTMNeighborFinder::Query class. This is a light-weight class generates
 * the sorted list of nearest neighbors of a particle.
 *
 * The PTMNeighborFinder class takes into account periodic boundary conditions. With periodic boundary conditions,
 * a particle can be appear multiple times in the neighbor list of another particle. Note, however, that a different neighbor *vector* is
 * reported for each periodic image of a neighbor.
 */
class OVITO_PARTICLES_EXPORT PTMNeighborFinder
{
public:
	//// Constructor
	PTMNeighborFinder() {}

	/// \brief Prepares the tree data structure.
	/// \param posProperty The positions of the particles.
	/// \param cellData The simulation cell data.
	/// \param selectionProperty Determines which particles are included in the neighbor search (optional).
	/// \param promis A callback object that will be used to the report progress.
	/// \return \c false when the operation has been canceled by the user;
	///         \c true on success.
	/// \throw Exception on error.
	bool prepare(NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> *neighQuery,
				 ConstPropertyAccess<qlonglong> correspondenceArray,
				 ConstPropertyAccess<PTMAlgorithm::StructureType> structuresArray,
				 Task* promise);

#if 0
	/// Returns the number of input particles in the system for which the PTMNeighborFinder was created.
	size_t particleCount() const {
		return atoms.size();
	}

	/// Returns the coordinates of the i-th input particle.
	const Point3& particlePos(size_t index) const {
		OVITO_ASSERT(index >= 0 && index < atoms.size());
		return atoms[index].pos;
	}
#endif

public:
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
	bool fill_neighbors(NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> *neighQuery,
						   size_t particleIndex,
						   size_t offset,
						   size_t num,
						   double* delta)
	{
		neighQuery->findNeighbors(particleIndex);
		int numNeighbors = neighQuery->results().size();

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
			_env.atom_indices[i + 1 + offset] = neighQuery->results()[p].index;
			_env.points[i + 1 + offset][0] = neighQuery->results()[p].delta.x() + delta[0];
			_env.points[i + 1 + offset][1] = neighQuery->results()[p].delta.y() + delta[1];
			_env.points[i + 1 + offset][2] = neighQuery->results()[p].delta.z() + delta[2];
		}

		return true;
	}

public:
	void findNeighbors(size_t particleIndex)
	{
		structureType = (PTMAlgorithm::StructureType)_structuresArray[particleIndex];
		int ptm_type = PTMAlgorithm::ovito_to_ptm_structure_type(structureType);

		_neighQuery->findNeighbors(particleIndex);
		int numNeighbors = _neighQuery->results().size();
		int num_inner = ptm_num_nbrs[ptm_type], num_outer = 0;

		int best_template_index = 0;
		if (ptm_type == PTM_MATCH_NONE) {
			for (int i=0;i<PTM_MAX_INPUT_POINTS;i++) {
				_env.correspondences[i] = i;
			}

			num_inner = numNeighbors;
		}
		else {
			numNeighbors = ptm_num_nbrs[ptm_type];
			ptm_decode_correspondences(ptm_type,
									   _correspondenceArray[particleIndex],
									   _env.correspondences,
									   &best_template_index);
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

private:
	NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> *_neighQuery;
	ConstPropertyAccess<qlonglong> _correspondenceArray;
	ConstPropertyAccess<PTMAlgorithm::StructureType> _structuresArray;
	ptm_atomicenv_t _env;
	std::vector<Neighbor> _results;

	/// The internal list of atoms.
	//std::vector<NeighborListAtom> atoms;
};

}	// End of namespace
}	// End of namespace



#if 0
	if (_structureType != OTHER && qtarget != nullptr) {

		//arrange orientation in PTM format
		double qtarget_ptm[4] = { qtarget->w(), qtarget->x(), qtarget->y(), qtarget->z() };

		double disorientation = 0;	//TODO: this is probably not needed
		int template_index = ptm_remap_template(type, true, _bestTemplateIndex, qtarget_ptm, _q, &disorientation, _env.correspondences, &_bestTemplate);
		if (template_index < 0)
			return _structureType;

		_bestTemplateIndex = template_index;
	}
#endif

#if 0
////--------replace this with saved PTM data--------------------

//todo: when getting stored PTM data (when PTM is not called here), assert that output_conventional_orientations=true was used.

	//TODO: don't hardcode input flags
	int32_t flags = PTM_CHECK_FCC | PTM_CHECK_DCUB | PTM_CHECK_GRAPHENE;// | PTM_CHECK_HCP | PTM_CHECK_BCC;

	// Call PTM library to identify local structure.
	int32_t type = PTM_MATCH_NONE, alloy_type = PTM_ALLOY_NONE;
	double scale, interatomic_distance;
	double rmsd;
	double q[4];
	int8_t correspondences[PolyhedralTemplateMatchingModifier::MAX_NEIGHBORS+1];
	ptm_index(handle, numNeighbors + 1, points, nullptr, flags, true,
			&type, &alloy_type, &scale, &rmsd, q, nullptr, nullptr, nullptr, nullptr,
			&interatomic_distance, nullptr, nullptr, nullptr, correspondences);
	if (rmsd > 0.1) {	//TODO: don't hardcode threshold
		type = PTM_MATCH_NONE;
		rmsd = INFINITY;
	}

	if (type == PTM_MATCH_NONE)
		return;
////------------------------------------------------------------

	double qmapped[4];
	double qtarget[4] = {qw, qx, qy, qz};	//PTM quaterion ordering
	const double (*best_template)[3] = NULL;
	int _template_index = ptm_remap_template(type, true, input_template_index, qtarget, q, qmapped, &disorientation, correspondences, &best_template);
	if (_template_index < 0)
		return;

	//arrange orientation in OVITO format
	qres[0] = qmapped[1];	//qx
	qres[1] = qmapped[2];	//qy
	qres[2] = qmapped[3];	//qz
	qres[3] = qmapped[0];	//qw

	//structure_type = type;
	template_index = _template_index;
	for (int i=0;i<ptm_num_nbrs[type];i++) {

		int index = correspondences[i+1] - 1;
		auto r = neighQuery.results()[index];

		Vector3 tcoords(best_template[i+1][0], best_template[i+1][1], best_template[i+1][2]);
		result.emplace_back(r.delta, tcoords, r.distanceSq, r.index);
	}
#endif

#if 0
		/// Returns the number of neighbors for the PTM structure found for the current particle.
		int numTemplateNeighbors() const;

		/// Returns the number of nearest neighbors found for the current particle.
		int numNearestNeighbors() const { return results().size(); }

		/// Returns the neighbor information for the i-th nearest neighbor of the current particle.
		const NearestNeighborFinder::Neighbor& getNearestNeighbor(int index) const {
			OVITO_ASSERT(index >= 0 && index < results().size());
			return results()[index];
		}

		/// Returns the neighbor information corresponding to the i-th neighbor in the PTM template
		/// identified for the current particle.
		const NearestNeighborFinder::Neighbor& getTemplateNeighbor(int index) const;

		/// Returns the ideal vector corresponding to the i-th neighbor in the PTM template
		/// identified for the current particle.
		const Vector_3<double>& getIdealNeighborVector(int index) const;
#endif
