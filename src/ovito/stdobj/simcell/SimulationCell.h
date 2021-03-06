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

#pragma once


#include <ovito/stdobj/StdObj.h>

namespace Ovito { namespace StdObj {

/**
* \brief Stores the geometry and boundary conditions of a simulation box.
 *
 * The simulation box geometry is a parallelepiped defined by three edge vectors.
 * A fourth vector specifies the origin of the simulation box in space.
 */
class OVITO_STDOBJ_EXPORT SimulationCell
{
public:

	/// Default constructor.
	SimulationCell() {
		_simulationCell = AffineTransformation::Zero();
		_reciprocalSimulationCell = AffineTransformation::Zero();
		_pbcFlags.fill(false);
		_is2D = false;
		_isValid = false;
	}

	/// Initialization constructor.
	SimulationCell(const AffineTransformation& cellMatrix, const std::array<bool,3>& pbcFlags, bool is2D = false) : _simulationCell(cellMatrix), _pbcFlags(pbcFlags), _is2D(is2D) {
		computeInverseMatrix();
	}

	/// Indicates that the simulation cell geometry has been specified and is non-degenerate. 
	bool isValid() const { return _isValid; }

	/// Returns whether this is a 2D system.
	bool is2D() const { return _is2D; }

	/// Sets whether this is a 2D system.
	void set2D(bool is2D) {
		_is2D = is2D;
		if(is2D) _pbcFlags[2] = false;
		computeInverseMatrix();
	}

	/// Returns the current simulation cell matrix.
	const AffineTransformation& matrix() const { return _simulationCell; }

	/// Returns the current reciprocal simulation cell matrix.
	const AffineTransformation& inverseMatrix() const { return _reciprocalSimulationCell; }

	/// Sets the simulation cell matrix.
	void setMatrix(const AffineTransformation& cellMatrix) {
		_simulationCell = cellMatrix;
		computeInverseMatrix();
	}

	/// Returns the PBC flags.
	const std::array<bool,3>& pbcFlags() const { return _pbcFlags; }

	/// Returns whether the simulation cell has periodic boundary conditions applied in the given direction.
	bool hasPbc(size_t dim) const { OVITO_ASSERT(dim < 3); return _pbcFlags[dim]; }

	/// Returns whether the simulation cell has periodic boundary conditions applied in at least one direction.
	bool hasPbc() const { return _pbcFlags[0] || _pbcFlags[1] || _pbcFlags[2]; }

	/// Sets the PBC flags.
	void setPbcFlags(const std::array<bool,3>& flags) { _pbcFlags = flags; }

	/// Sets the PBC flags.
	void setPbcFlags(bool pbcX, bool pbcY, bool pbcZ) { _pbcFlags[0] = pbcX; _pbcFlags[1] = pbcY; _pbcFlags[2] = pbcZ; }

	/// Computes the (positive) volume of the three-dimensional cell.
	FloatType volume3D() const {
		return std::abs(_simulationCell.determinant());
	}

	/// Computes the (positive) volume of the two-dimensional cell.
	FloatType volume2D() const {
		return _simulationCell.column(0).cross(_simulationCell.column(1)).length();
	}

	/// Returns true if the three edges of the cell are parallel to the three
	/// coordinates axes.
	bool isAxisAligned() const {
		if(matrix()(1,0) != 0 || matrix()(2,0) != 0) return false;
		if(matrix()(0,1) != 0 || matrix()(2,1) != 0) return false;
		if(matrix()(0,2) != 0 || matrix()(1,2) != 0) return false;
		return true;
	}

	/// Checks if two simulation cells are identical.
	bool operator==(const SimulationCell& other) const {
		return (_simulationCell == other._simulationCell && _pbcFlags == other._pbcFlags && _is2D == other._is2D);
	}

	/// Converts a point given in reduced cell coordinates to a point in absolute coordinates.
	Point3 reducedToAbsolute(const Point3& reducedPoint) const { return _simulationCell * reducedPoint; }

	/// Converts a point given in absolute coordinates to a point in reduced cell coordinates.
	Point3 absoluteToReduced(const Point3& absPoint) const { return _reciprocalSimulationCell * absPoint; }

	/// Converts a vector given in reduced cell coordinates to a vector in absolute coordinates.
	Vector3 reducedToAbsolute(const Vector3& reducedVec) const { return _simulationCell * reducedVec; }

	/// Converts a vector given in absolute coordinates to a point in vector cell coordinates.
	Vector3 absoluteToReduced(const Vector3& absVec) const { return _reciprocalSimulationCell * absVec; }

	/// Wraps a point at the periodic boundaries of the cell.
	Point3 wrapPoint(const Point3& p) const {
		Point3 pout = p;
		for(size_t dim = 0; dim < 3; dim++) {
			if(_pbcFlags[dim]) {
				if(FloatType s = std::floor(_reciprocalSimulationCell.prodrow(p, dim)))
					pout -= s * _simulationCell.column(dim);
			}
		}
		return pout;
	}

	/// Wraps a vector at the periodic boundaries of the cell using minimum image convention.
	Vector3 wrapVector(const Vector3& v) const {
		Vector3 vout = v;
		for(size_t dim = 0; dim < 3; dim++) {
			if(_pbcFlags[dim]) {
				if(FloatType s = std::floor(_reciprocalSimulationCell.prodrow(v, dim) + FloatType(0.5)))
					vout -= s * _simulationCell.column(dim);
			}
		}
		return vout;
	}

	/// Calculates the normal vector of the given simulation cell side.
	Vector3 cellNormalVector(size_t dim) const {
		Vector3 normal = _simulationCell.column((dim+1)%3).cross(_simulationCell.column((dim+2)%3));
		// Flip normal if necessary.
		if(normal.dot(_simulationCell.column(dim)) < 0)
			return normal / (-normal.length());
		else
			return normal.safelyNormalized();
	}

	/// Tests if a vector so long that it would be wrapped at a periodic boundary when using the minimum image convention.
	bool isWrappedVector(const Vector3& v) const {
		for(size_t dim = 0; dim < 3; dim++) {
			if(_pbcFlags[dim]) {
				if(std::abs(_reciprocalSimulationCell.prodrow(v, dim)) >= FloatType(0.5))
					return true;
			}
		}
		return false;
	}

	/// \brief Helper function that computes the modulo operation for two integer numbers k and n.
	///
	/// This function can handle negative numbers k. This allows mapping any number k that is
	/// outside the interval [0,n) back into the interval. Use this to implement periodic boundary conditions.
	static inline int modulo(int k, int n) {
		return ((k %= n) < 0) ? k+n : k;
	}

	/// \brief Helper function that computes the modulo operation for two floating-point numbers k and n.
	///
	/// This function can handle negative numbers k. This allows mapping any number k that is
	/// outside the interval [0,n) back into the interval. Use this to implement periodic boundary conditions.
	static inline FloatType modulo(FloatType k, FloatType n) {
		k = std::fmod(k, n);
		return (k < 0) ? k+n : k;
	}

private:

	/// Computes the inverse of the cell matrix.
	void computeInverseMatrix() {
		if(!is2D()) {
			_isValid = _simulationCell.inverse(_reciprocalSimulationCell);
			if(!_isValid) {
				_reciprocalSimulationCell.setIdentity();
				_pbcFlags.fill(false);
			}
		}
		else {
			_reciprocalSimulationCell.setIdentity();
			FloatType det = _simulationCell(0,0)*_simulationCell(1,1) - _simulationCell(0,1)*_simulationCell(1,0);
			_isValid = (std::abs(det) > FLOATTYPE_EPSILON);
			if(_isValid) {
				_reciprocalSimulationCell(0,0) = _simulationCell(1,1) / det;
				_reciprocalSimulationCell(1,0) = -_simulationCell(1,0) / det;
				_reciprocalSimulationCell(0,1) = -_simulationCell(0,1) / det;
				_reciprocalSimulationCell(1,1) = _simulationCell(0,0) / det;
				_reciprocalSimulationCell.translation().x() = -(_reciprocalSimulationCell(0,0) * _simulationCell.translation().x() + _reciprocalSimulationCell(0,1) * _simulationCell.translation().y());
				_reciprocalSimulationCell.translation().y() = -(_reciprocalSimulationCell(1,0) * _simulationCell.translation().x() + _reciprocalSimulationCell(1,1) * _simulationCell.translation().y());
			}
			else _pbcFlags.fill(false);
		}
	}

	/// The geometry of the cell.
	AffineTransformation _simulationCell;

	/// The reciprocal cell matrix.
	AffineTransformation _reciprocalSimulationCell;

	/// PBC flags.
	std::array<bool,3> _pbcFlags;

	/// Indicates that it is a 2D system.
	bool _is2D;

	/// Indicates whether the simulation cell geometry has been specified.
	bool _isValid;
};

}	// End of namespace
}	// End of namespace
