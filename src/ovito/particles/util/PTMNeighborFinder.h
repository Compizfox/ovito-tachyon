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
 * \brief This utility class finds the neighbors of a particle whose local crystalline order has been determined
 *        with the PolyhedralTemplateMatching modifier. In order to use this class, the "output_orientation"
 *        parameter must be set (in the scripting interface) or "Lattice orientations" (in the GUI).
 *
 * The PTMNeighborFinder class must be initialized by a call to prepare().
 *
 * After the PTMNeighborFinder has been initialized, one can find the nearest neighbors of some central
 * particle by constructing an instance of the PTMNeighborFinder::Query class. This class also contains some
 * properties computed by PTM.
 */
class OVITO_PARTICLES_EXPORT PTMNeighborFinder
{
public:
	//// Constructor
	/// \param _all_properties Determines whether to compute disorientations between neighbors.
	PTMNeighborFinder(bool _all_properties)
	{
		// scripting interface should always set this to true;
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
		std::vector<Vector3> rotatedTemplate;
		centered.push_back(Vector_3<double>(0, 0, 0));
		rotatedTemplate.push_back(Vector_3<double>(0, 0, 0));
		Vector3 barycenter = Vector_3<double>(0, 0, 0);

		for (auto nbr: _results) {
			rotatedTemplate.push_back(orientation * nbr.idealVector);
			centered.push_back(nbr.delta);
			barycenter += nbr.delta;
		}

		barycenter /= centered.size();
		for (int i=0;i<centered.size();i++) {
			centered[i] -= barycenter;
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
	void findNeighbors(size_t particleIndex, Quaternion* targetOrientation)
	{
		structureType = _structuresArray[particleIndex];
		orientation = _orientationsArray[particleIndex];
		rmsd = std::numeric_limits<FloatType>::infinity();

		int ptm_type = PTMAlgorithm::ovito_to_ptm_structure_type(structureType);
		getNeighbors(particleIndex, ptm_type);

		int8_t remap_permutation[PTM_MAX_INPUT_POINTS];
		for (int i=0;i<PTM_MAX_INPUT_POINTS;i++) {
			remap_permutation[i] = i;
		}

		if (structureType != PTMAlgorithm::OTHER && targetOrientation != nullptr) {
			//arrange orientation in PTM format
			double qtarget[4] = {targetOrientation->w(),
								 targetOrientation->x(),
								 targetOrientation->y(),
								 targetOrientation->z()};
			double qptm[4] = {orientation.w(),
							  orientation.x(),
							  orientation.y(),
							  orientation.z()};
			double dummy = 0;
			_templateIndex = ptm_remap_template(ptm_type, true, _templateIndex,
												qtarget, qptm, &dummy,
												remap_permutation, nullptr);
			//arrange orientation in OVITO format
			orientation.w() = qptm[0];
			orientation.x() = qptm[1];
			orientation.y() = qptm[2];
			orientation.z() = qptm[3];
		}

		const double (*ptmTemplate)[3] = PTMAlgorithm::get_template(structureType, _templateIndex);
		_results.clear();
		for (int i=0;i<_numNeighbors;i++) {
			Neighbor n;
			int index = remap_permutation[i + 1];
			n.index = _env.atom_indices[index];
			double* p = _env.points[index];
			n.delta = Vector3(p[0], p[1], p[2]);
			n.distanceSq = n.delta.squaredLength();

			if (structureType == PTMAlgorithm::OTHER) {
				n.idealVector = Vector_3<double>(0, 0, 0);
			}
			else {
				const double* q = ptmTemplate[i + 1];
				n.idealVector = Vector3(q[0], q[1], q[2]);
			}

			if (all_properties && structureType != PTMAlgorithm::OTHER) {
				n.disorientation = PTMAlgorithm::calculate_disorientation(structureType,
																		  _structuresArray[n.index],
																		  orientation,
																		  _orientationsArray[n.index]);
			}
			else {
				n.disorientation = std::numeric_limits<FloatType>::max();
			}

			//n.atom = atom;
			_results.push_back(n);
		}

		if (structureType != PTMAlgorithm::OTHER) {
			calculate_rmsd_scale(ptmTemplate);
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

