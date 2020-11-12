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
	PTMNeighborFinder(bool _all_properties)
	{
		all_properties = _all_properties;
	}

	/// \brief Prepares the tree data structure.
	/// \param posProperty The positions of the particles.
	/// \param cellData The simulation cell data.
	/// \param selectionProperty Determines which particles are included in the neighbor search (optional).
	/// \param promise A callback object that will be used to the report progress.
	/// \return \c false when the operation has been canceled by the user;
	///         \c true on success.
	/// \throw Exception on error.
	bool prepare(NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> *neighQuery,
				 ConstPropertyAccess<PTMAlgorithm::StructureType> structuresArray,
				 ConstPropertyAccess<Quaternion> orientationsArray,
				 ConstPropertyAccess<qlonglong> correspondencesArray,
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
	Quaternion orientation;
	FloatType rmsd;

	/// Contains information about a single neighbor of the central particle.
	struct Neighbor
	{
		Vector3 delta;
		FloatType distanceSq;
		Vector3 idealVector;
		//NeighborListAtom* atom;
		size_t index;
		FloatType disorientation;
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

	bool getNeighbors(size_t particleIndex, int ptm_type)
	{
		_neighQuery->findNeighbors(particleIndex);
		_numNeighbors = _neighQuery->results().size();
		_templateIndex = 0;

		int num_inner = ptm_num_nbrs[ptm_type], num_outer = 0;
		if (ptm_type == PTM_MATCH_NONE) {
			for (int i=0;i<PTM_MAX_INPUT_POINTS;i++) {
				_env.correspondences[i] = i;
			}

			num_inner = _numNeighbors;
		}
		else {
			_numNeighbors = ptm_num_nbrs[ptm_type];
			ptm_decode_correspondences(ptm_type,
									   _correspondencesArray[particleIndex],
									   _env.correspondences,
									   &_templateIndex);
		}

		_env.num = _numNeighbors + 1;

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
				fill_neighbors(_neighQuery,
							   _env.atom_indices[1 + i],
							   num_inner + i * num_outer,
							   num_outer,
							   _env.points[i + 1]);
			}
		}

		return true;
	}

	void calculate_rmsd_scale(const double (*ptmTemplate)[3])
	{
		// get neighbor points
		std::vector<Vector3> centered;
		for (int i=0;i<_numNeighbors+1;i++) {
			Vector3 v = *reinterpret_cast<const Vector_3<double>*>(_env.points[i]);
			centered.push_back(v);
		}

		Vector3 barycenter = Vector_3<double>(0, 0, 0);
		for (auto v: centered) {
			barycenter += v;
		}

		barycenter /= centered.size();
		for (int i=0;i<centered.size();i++) {
			centered[i] -= barycenter;
		}

		// get template points
		std::vector<Vector3> rotatedTemplate;
		for (int i=0;i<_numNeighbors+1;i++) {
			Vector3 v = *reinterpret_cast<const Vector_3<double>*>(ptmTemplate[i]);
			rotatedTemplate.push_back(orientation * v);
		}

		// calculate scale
		// (s.a - b)^2 = s^2.a^2 - 2.s.a.b + b^2
		// d/ds (s^2.a^2 - 2.s.a.b + b^2) = 2.s.a^2 - 2.a.b
		// s.a^2 = a.b
		// s = a.b / (a.a)
		FloatType numerator = 0, denominator = 0;
		for (int i=0;i<centered.size();i++) {
			numerator += centered[i].dot(rotatedTemplate[i]);
			denominator += centered[i].squaredLength();
		}
		FloatType scale = numerator / denominator;

		// calculate RMSD
		rmsd = 0;
		for (int i=0;i<centered.size();i++) {
			auto delta = scale * centered[i] - rotatedTemplate[i];
			rmsd += delta.squaredLength();
		}
		rmsd = sqrt(rmsd / centered.size());
	}

public:
	void findNeighbors(size_t particleIndex)
	{
		structureType = _structuresArray[particleIndex];
		orientation = _orientationsArray[particleIndex];
		rmsd = std::numeric_limits<FloatType>::infinity();

		int ptm_type = PTMAlgorithm::ovito_to_ptm_structure_type(structureType);
		getNeighbors(particleIndex, ptm_type);
		const double (*ptmTemplate)[3] = PTMAlgorithm::get_template(structureType, _templateIndex);;

		if (structureType != PTMAlgorithm::OTHER) {
			calculate_rmsd_scale(ptmTemplate);
		}

		_results.clear();
		for (int i=0;i<_numNeighbors;i++) {
			Neighbor n;
			double* point = _env.points[i + 1];
			n.delta = Vector3(point[0], point[1], point[2]);
			n.distanceSq = n.delta.squaredLength();
			n.index = _env.atom_indices[1 + i];

			if (structureType == PTMAlgorithm::OTHER || !all_properties) {
				n.disorientation = std::numeric_limits<FloatType>::max();
				n.idealVector = Vector_3<double>(0, 0, 0);
			}
			else {
				n.disorientation = PTMAlgorithm::calculate_disorientation(structureType,
																		  _structuresArray[n.index],
																		  orientation,
																		  _orientationsArray[n.index]);
				int index = _env.correspondences[i + 1] - 1;
				n.idealVector = *reinterpret_cast<const Vector_3<double>*>(ptmTemplate[index + 1]);
			}

			//n.atom = atom;
			_results.push_back(n);
		}
	}

	std::vector<Neighbor> results() {return _results;}

private:
	bool all_properties;
	NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> *_neighQuery;
	ConstPropertyAccess<PTMAlgorithm::StructureType> _structuresArray;
	ConstPropertyAccess<Quaternion> _orientationsArray;
	ConstPropertyAccess<qlonglong> _correspondencesArray;
	int _numNeighbors;
	int _templateIndex;
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
		double qmapped[4];
		const double (*best_template)[3] = NULL;
		double disorientation = 0;	//TODO: this is probably not needed
		int _template_index = ptm_remap_template(type, true, input_template_index, qtarget_ptm, q, qmapped, &disorientation, correspondences, &best_template);
		if (_template_index < 0)
			return;

		//arrange orientation in OVITO format
		qres[0] = qmapped[1];	//qx
		qres[1] = qmapped[2];	//qy
		qres[2] = qmapped[3];	//qz
		qres[3] = qmapped[0];	//qw
	}
#endif
