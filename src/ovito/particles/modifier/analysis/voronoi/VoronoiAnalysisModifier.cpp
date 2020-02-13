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
#include <ovito/particles/util/NearestNeighborFinder.h>
#include <ovito/particles/objects/BondsObject.h>
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include <ovito/core/utilities/concurrent/ParallelFor.h>
#include <ovito/core/utilities/units/UnitsManager.h>
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include "VoronoiAnalysisModifier.h"

#include <voro++.hh>

namespace Ovito { namespace Particles {

constexpr int VoronoiAnalysisModifier::VoronoiAnalysisEngine::FaceOrderStorageLimit;

IMPLEMENT_OVITO_CLASS(VoronoiAnalysisModifier);
DEFINE_PROPERTY_FIELD(VoronoiAnalysisModifier, onlySelected);
DEFINE_PROPERTY_FIELD(VoronoiAnalysisModifier, useRadii);
DEFINE_PROPERTY_FIELD(VoronoiAnalysisModifier, computeIndices);
DEFINE_PROPERTY_FIELD(VoronoiAnalysisModifier, computeBonds);
DEFINE_PROPERTY_FIELD(VoronoiAnalysisModifier, computePolyhedra);
DEFINE_PROPERTY_FIELD(VoronoiAnalysisModifier, edgeThreshold);
DEFINE_PROPERTY_FIELD(VoronoiAnalysisModifier, faceThreshold);
DEFINE_PROPERTY_FIELD(VoronoiAnalysisModifier, relativeFaceThreshold);
DEFINE_REFERENCE_FIELD(VoronoiAnalysisModifier, bondsVis);
DEFINE_REFERENCE_FIELD(VoronoiAnalysisModifier, polyhedraVis);
SET_PROPERTY_FIELD_LABEL(VoronoiAnalysisModifier, onlySelected, "Use only selected particles");
SET_PROPERTY_FIELD_LABEL(VoronoiAnalysisModifier, useRadii, "Use particle radii");
SET_PROPERTY_FIELD_LABEL(VoronoiAnalysisModifier, computeIndices, "Compute Voronoi indices");
SET_PROPERTY_FIELD_LABEL(VoronoiAnalysisModifier, computeBonds, "Generate neighbor bonds");
SET_PROPERTY_FIELD_LABEL(VoronoiAnalysisModifier, computePolyhedra, "Generate Voronoi polyhedra");
SET_PROPERTY_FIELD_LABEL(VoronoiAnalysisModifier, edgeThreshold, "Edge length threshold");
SET_PROPERTY_FIELD_LABEL(VoronoiAnalysisModifier, faceThreshold, "Absolute face area threshold");
SET_PROPERTY_FIELD_LABEL(VoronoiAnalysisModifier, relativeFaceThreshold, "Relative face area threshold");
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(VoronoiAnalysisModifier, edgeThreshold, WorldParameterUnit, 0);
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(VoronoiAnalysisModifier, faceThreshold, FloatParameterUnit, 0);
SET_PROPERTY_FIELD_UNITS_AND_RANGE(VoronoiAnalysisModifier, relativeFaceThreshold, PercentParameterUnit, 0, 1);

/******************************************************************************
* Constructs the modifier object.
******************************************************************************/
VoronoiAnalysisModifier::VoronoiAnalysisModifier(DataSet* dataset) : AsynchronousModifier(dataset),
	_onlySelected(false),
	_useRadii(false),
	_edgeThreshold(0),
	_faceThreshold(0),
	_computeIndices(false),
	_computeBonds(false),
	_computePolyhedra(false),
	_relativeFaceThreshold(0)
{
	// Create the vis element for rendering the bonds generated by the modifier.
	setBondsVis(new BondsVis(dataset));

	// Create the vis element for rendering the Coronoi polyhedra generated by the modifier.
	setPolyhedraVis(new SurfaceMeshVis(dataset));
	polyhedraVis()->setShowCap(false);
	polyhedraVis()->setSmoothShading(false);
	polyhedraVis()->setSurfaceTransparency(FloatType(0.25));
	polyhedraVis()->setHighlightEdges(true);
	polyhedraVis()->setObjectTitle(tr("Voronoi polyhedra"));
}

/******************************************************************************
* Asks the modifier whether it can be applied to the given input data.
******************************************************************************/
bool VoronoiAnalysisModifier::OOMetaClass::isApplicableTo(const DataCollection& input) const
{
	return input.containsObject<ParticlesObject>();
}

/******************************************************************************
* Creates and initializes a computation engine that will compute the modifier's results.
******************************************************************************/
Future<AsynchronousModifier::ComputeEnginePtr> VoronoiAnalysisModifier::createEngine(const PipelineEvaluationRequest& request, ModifierApplication* modApp, const PipelineFlowState& input)
{
	// Get the input particles.
	const ParticlesObject* particles = input.expectObject<ParticlesObject>();
	particles->verifyIntegrity();
	const PropertyObject* posProperty = particles->expectProperty(ParticlesObject::PositionProperty);

	// Get simulation cell.
	const SimulationCellObject* inputCell = input.expectObject<SimulationCellObject>();
	if(inputCell->is2D())
		throwException(tr("The Voronoi modifier does not support 2d simulation cells."));

	// Get selection particle property.
	ConstPropertyPtr selectionProperty;
	if(onlySelected())
		selectionProperty = particles->expectProperty(ParticlesObject::SelectionProperty)->storage();

	// Get particle radii.
	std::vector<FloatType> radii;
	if(useRadii())
		radii = particles->inputParticleRadii();

	// The Voro++ library uses 32-bit integers. It cannot handle more than 2^31 input points.
	if(posProperty->size() > (size_t)std::numeric_limits<int>::max())
		throwException(tr("Voronoi analysis modifier is limited to a maximum of %1 particles in the current program version.").arg(std::numeric_limits<int>::max()));

	// Create engine object. Pass all relevant modifier parameters to the engine as well as the input data.
	return std::make_shared<VoronoiAnalysisEngine>(
			input.stateValidity(),
			particles,
			posProperty->storage(),
			selectionProperty,
			particles->getPropertyStorage(ParticlesObject::IdentifierProperty),
			std::move(radii),
			inputCell->data(),
			computeIndices(),
			computeBonds(),
			computePolyhedra(),
			edgeThreshold(),
			faceThreshold(),
			relativeFaceThreshold());
}

/******************************************************************************
* Performs the actual computation. This method is executed in a worker thread.
******************************************************************************/
void VoronoiAnalysisModifier::VoronoiAnalysisEngine::perform()
{
	task()->setProgressText(tr("Performing Voronoi analysis"));
	task()->beginProgressSubSteps(_computePolyhedra ? 2 : 1);

	// Stores the starting vertex index and the vertex count for each Voronoi polyhedron. 
	std::vector<std::pair<SurfaceMeshData::vertex_index, SurfaceMeshData::size_type>> polyhedraVertices;

	if(_computePolyhedra) {

		// Create the "Region" mesh face property.
		_polyhedraMesh.createFaceProperty(SurfaceMeshFaces::RegionProperty, false);

		// Create the "Adjacent Cell" face property, which stores the index of the neighboring Voronoi cell.
		_adjacentCellProperty = _polyhedraMesh.createFaceProperty(PropertyStorage::Int, 1, 0, QStringLiteral("Adjacent Cell"), false).get();

		// Create as many mesh regions as there are input particles.
		_polyhedraMesh.createRegions(_positions->size());
		polyhedraVertices.resize(_polyhedraMesh.regionCount());

		// Create the "Particle Identifier" region property, which indicates the ID of the particles that are at the center of each Voronoi polyhedron.
		_centerParticleProperty = _polyhedraMesh.createRegionProperty(PropertyStorage::Int64, 1, 0, QStringLiteral("Particle Identifier"), true).get();
		if(_particleIdentifiers) {
			OVITO_ASSERT(_centerParticleProperty->size() == _particleIdentifiers->size());
			boost::copy(ConstPropertyAccess<qlonglong>(_particleIdentifiers), PropertyAccess<qlonglong>(_centerParticleProperty).begin());
		}
		else {
			boost::algorithm::iota_n(PropertyAccess<qlonglong>(_centerParticleProperty).begin(), (qlonglong)1, _centerParticleProperty->size());
		}

		// Create the "Volume" region property, which stores the volume of each Voronoi cell.
		_cellVolumeProperty = _polyhedraMesh.createRegionProperty(SurfaceMeshRegions::VolumeProperty, true).get();

		// Create the "Coordination" region property, which stores the number of faces of each Voronoi cell.
		_cellCoordinationProperty = _polyhedraMesh.createRegionProperty(PropertyStorage::Int, 1, 0, QStringLiteral("Coordination"), true).get();

		// Create the "Surface Area" region property, which stores the face area of each Voronoi cell.
		_surfaceAreaProperty = _polyhedraMesh.createRegionProperty(SurfaceMeshRegions::SurfaceAreaProperty, true).get();
	}

	if(_positions->size() == 0 || _simCell.volume3D() == 0) {
		if(maxFaceOrders()) {
			_voronoiIndices = std::make_shared<PropertyStorage>(_positions->size(), PropertyStorage::Int, 3, 0, QStringLiteral("Voronoi Index"), true);
			// Re-use the output particle property as an output mesh region property.
			if(_computePolyhedra) {
				_polyhedraMesh.addRegionProperty(voronoiIndices());
				_polyhedraMesh.addRegionProperty(maxFaceOrders());
			}
		}
		// Nothing else to do if there are no particles.
		return;
	}

	// The squared edge length threshold.
	// Apply additional prefactor of 4, because Voronoi cell vertex coordinates are all scaled by factor of 2.
	const FloatType sqEdgeThreshold = _edgeThreshold * _edgeThreshold * 4;

	// Prepare output data arrays.
	PropertyAccess<FloatType> atomicVolumesArray(atomicVolumes());
	PropertyAccess<int> coordinationNumbersArray(coordinationNumbers());
	PropertyAccess<int> maxFaceOrdersArray(maxFaceOrders());
	PropertyAccess<int> adjacentCellArray(_adjacentCellProperty);
	PropertyAccess<FloatType> cellVolumeArray(_cellVolumeProperty);
	PropertyAccess<FloatType> surfaceAreaArray(_surfaceAreaProperty);
	PropertyAccess<int> cellCoordinationArray(_cellCoordinationProperty);

	// Prepare input data array.
	ConstPropertyAccess<int> selectionArray(_selection);
	ConstPropertyAccess<Point3> positionsArray(_positions);

	auto processCell = [&](voro::voronoicell_neighbor& v, size_t index,
		std::vector<int>& voronoiBuffer, std::vector<size_t>& voronoiBufferIndex, QMutex* bondMutex)
	{
		// Compute cell volume.
		double vol = v.volume();
		atomicVolumesArray[index] = (FloatType)vol;

		// Accumulate total volume of Voronoi cells.
		// Loop is for lock-free write access to shared max counter.
		double prevVolumeSum = voronoiVolumeSum();
		while(!voronoiVolumeSum().compare_exchange_weak(prevVolumeSum, prevVolumeSum + vol));

		// Compute total surface area of Voronoi cell when relative area threshold is used to
		// filter out small faces.
		double faceAreaThreshold = _faceThreshold;
		if(_relativeFaceThreshold > 0)
			faceAreaThreshold = std::max(v.surface_area() * _relativeFaceThreshold, faceAreaThreshold);

		int localMaxFaceOrder = 0;
		int localVoronoiIndex[FaceOrderStorageLimit] = {0};
		int coordNumber = 0;
		FloatType cellFaceArea = 0;

		// Create Voronoi cell mesh vertices.
		SurfaceMeshData::vertex_index meshVertexBaseIndex;
		SurfaceMeshData::region_index meshRegionIndex = index;
		if(_computePolyhedra) {
			const Point3& center = positionsArray[index];
			QMutexLocker locker(bondMutex);
			cellVolumeArray[meshRegionIndex] = vol;
			meshVertexBaseIndex = _polyhedraMesh.vertexCount();
			const double* ptsp = v.pts;
			for(int i = 0; i < v.p; i++, ptsp += 3) {
				_polyhedraMesh.createVertex(Point3(center.x() + 0.5*ptsp[0], center.y() + 0.5*ptsp[1], center.z() + 0.5*ptsp[2]));
			}
			// Store the base vertex index and the vertex count in the look-up map.
			polyhedraVertices[meshRegionIndex].first = meshVertexBaseIndex;
			polyhedraVertices[meshRegionIndex].second = v.p;
		}

		// Iterate over the Voronoi faces and their edges.
		for(int i = 1; i < v.p; i++) {
			for(int j = 0; j < v.nu[i]; j++) {
				int k = v.ed[i][j];
				if(k >= 0) {
					int neighbor_id = v.ne[i][j];
					int faceOrder = 0;
					FloatType area = 0;

					// Create Voronoi cell mesh face.
					SurfaceMeshData::face_index meshFace;
					if(_computePolyhedra) {
						QMutexLocker locker(bondMutex);
						meshFace = _polyhedraMesh.createFace({}, meshRegionIndex);
						adjacentCellArray[meshFace] = neighbor_id;
						_polyhedraMesh.createEdge(meshVertexBaseIndex + i, meshVertexBaseIndex + k, meshFace);
					}

					// Compute length of first face edge.
					Vector3 d(v.pts[3*k] - v.pts[3*i], v.pts[3*k+1] - v.pts[3*i+1], v.pts[3*k+2] - v.pts[3*i+2]);
					if(d.squaredLength() > sqEdgeThreshold)
						faceOrder++;
					v.ed[i][j] = -1 - k;
					int l = v.cycle_up(v.ed[i][v.nu[i]+j], k);
					do {
						int m = v.ed[k][l];
						if(_computePolyhedra) {
							QMutexLocker locker(bondMutex);
							_polyhedraMesh.createEdge(meshVertexBaseIndex + k, meshVertexBaseIndex + m, meshFace);
						}
						// Compute length of current edge.
						if(sqEdgeThreshold != 0) {
							Vector3 u(v.pts[3*m] - v.pts[3*k], v.pts[3*m+1] - v.pts[3*k+1], v.pts[3*m+2] - v.pts[3*k+2]);
							if(u.squaredLength() > sqEdgeThreshold)
								faceOrder++;
						}
						else faceOrder++;
						if(faceAreaThreshold != 0 || _computePolyhedra) {
							Vector3 w(v.pts[3*m] - v.pts[3*i], v.pts[3*m+1] - v.pts[3*i+1], v.pts[3*m+2] - v.pts[3*i+2]);
							area += d.cross(w).length() / 8;
							d = w;
						}
						v.ed[k][l] = -1 - m;
						l = v.cycle_up(v.ed[k][v.nu[k]+l], m);
						k = m;
					}
					while(k != i);
					cellFaceArea += area;

					if((faceAreaThreshold == 0 || area > faceAreaThreshold) && faceOrder >= 3) {
						coordNumber++;
						if(faceOrder > localMaxFaceOrder)
							localMaxFaceOrder = faceOrder;
						faceOrder--;
						if(maxFaceOrders() && faceOrder < FaceOrderStorageLimit)
							localVoronoiIndex[faceOrder]++;
						if(_computeBonds && neighbor_id >= 0 && neighbor_id != index) {
							OVITO_ASSERT(neighbor_id < _positions->size());
							Vector3 delta = positionsArray[index] - positionsArray[neighbor_id];
							Vector3I pbcShift = Vector3I::Zero();
							for(size_t dim = 0; dim < 3; dim++) {
								if(_simCell.pbcFlags()[dim])
									pbcShift[dim] = (int)floor(_simCell.inverseMatrix().prodrow(delta, dim) + FloatType(0.5));
							}
							Bond bond = { index, (size_t)neighbor_id, pbcShift };
							if(!bond.isOdd()) {
								QMutexLocker locker(bondMutex);
								bonds().push_back(bond);
							}
						}
					}
				}
			}
		}

		// Store computed result.
		coordinationNumbersArray[index] = coordNumber;
		if(maxFaceOrdersArray) {
			maxFaceOrdersArray[index] = localMaxFaceOrder;
			voronoiBufferIndex.push_back(index);
			voronoiBuffer.insert(voronoiBuffer.end(), localVoronoiIndex, localVoronoiIndex + std::min(localMaxFaceOrder, FaceOrderStorageLimit));
		}
		if(_computePolyhedra) {
			surfaceAreaArray[meshRegionIndex] = cellFaceArea;
			cellCoordinationArray[meshRegionIndex] = coordNumber;
		}

		// Keep track of the maximum number of edges per face.
		// Loop is for lock-free write access to shared max counter.
		int prevMaxFaceOrder = maxFaceOrder();
		while(localMaxFaceOrder > prevMaxFaceOrder && !maxFaceOrder().compare_exchange_weak(prevMaxFaceOrder, localMaxFaceOrder));
	};

	std::vector<int> voronoiBuffer;
	std::vector<size_t> voronoiBufferIndex;

	// Decide whether to use Voro++ container class or our own implementation.
	if(_simCell.isAxisAligned()) {
		// Use Voro++ container.
		double ax = _simCell.matrix()(0,3);
		double ay = _simCell.matrix()(1,3);
		double az = _simCell.matrix()(2,3);
		double bx = ax + _simCell.matrix()(0,0);
		double by = ay + _simCell.matrix()(1,1);
		double bz = az + _simCell.matrix()(2,2);
		if(ax > bx) std::swap(ax,bx);
		if(ay > by) std::swap(ay,by);
		if(az > bz) std::swap(az,bz);
		double volumePerCell = (bx - ax) * (by - ay) * (bz - az) * voro::optimal_particles / _positions->size();
		double cellSize = pow(volumePerCell, 1.0/3.0);
		int nx = (int)std::ceil((bx - ax) / cellSize);
		int ny = (int)std::ceil((by - ay) / cellSize);
		int nz = (int)std::ceil((bz - az) / cellSize);

		size_t count = 0;
		if(_radii.empty()) {
			// All particles have a uniform size.
			voro::container voroContainer(ax, bx, ay, by, az, bz, nx, ny, nz,
					_simCell.pbcFlags()[0], _simCell.pbcFlags()[1], _simCell.pbcFlags()[2], (int)std::ceil(voro::optimal_particles));

			// Insert particles into Voro++ container.
			for(size_t index = 0; index < positionsArray.size(); index++) {
				// Skip unselected particles (if requested).
				if(selectionArray && selectionArray[index] == 0)
					continue;
				const Point3& p = positionsArray[index];
				voroContainer.put(index, p.x(), p.y(), p.z());
				count++;
			}
			if(!count) return;

			task()->setProgressValue(0);
			task()->setProgressMaximum(count);

			voro::c_loop_all cl(voroContainer);
			voro::voronoicell_neighbor v;
			if(cl.start()) {
				do {
					if(!task()->incrementProgressValue())
						return;
					if(!voroContainer.compute_cell(v,cl))
						continue;
					processCell(v, cl.pid(), voronoiBuffer, voronoiBufferIndex, nullptr);
					count--;
				}
				while(cl.inc());
			}
			if(count)
				throw Exception(tr("Could not compute Voronoi cell for some particles."));
		}
		else {
			// Particles have non-uniform sizes -> Compute polydisperse Voronoi tessellation.
			voro::container_poly voroContainer(ax, bx, ay, by, az, bz, nx, ny, nz,
					_simCell.pbcFlags()[0], _simCell.pbcFlags()[1], _simCell.pbcFlags()[2],
					(int)std::ceil(voro::optimal_particles));

			// Insert particles into Voro++ container.
			for(size_t index = 0; index < positionsArray.size(); index++) {
				// Skip unselected particles (if requested).
				if(selectionArray && selectionArray[index] == 0)
					continue;
				const Point3& p = positionsArray[index];
				voroContainer.put(index, p.x(), p.y(), p.z(), _radii[index]);
				count++;
			}

			if(!count) return;
			task()->setProgressValue(0);
			task()->setProgressMaximum(count);

			voro::c_loop_all cl(voroContainer);
			voro::voronoicell_neighbor v;
			if(cl.start()) {
				do {
					if(!task()->incrementProgressValue())
						return;
					if(!voroContainer.compute_cell(v,cl))
						continue;
					processCell(v, cl.pid(), voronoiBuffer, voronoiBufferIndex, nullptr);
					count--;
				}
				while(cl.inc());
			}
		}
		if(count)
			throw Exception(tr("Could not compute Voronoi cell for some particles."));
	}
	else {
		// Special code path for non-orthogonal simulation cells:

		// Prepare the nearest neighbor list generator.
		NearestNeighborFinder nearestNeighborFinder;
		if(!nearestNeighborFinder.prepare(positions(), _simCell, selection(), task().get()))
			return;

		// Squared particle radii (input was just radii).
		for(auto& r : _radii)
			r = r*r;

		// This is the size we use to initialize Voronoi cells. Must be larger than the simulation box.
		double boxDiameter = sqrt(
				  _simCell.matrix().column(0).squaredLength()
				+ _simCell.matrix().column(1).squaredLength()
				+ _simCell.matrix().column(2).squaredLength());

		// The normal vectors of the three cell planes.
		std::array<Vector3,3> planeNormals;
		planeNormals[0] = _simCell.cellNormalVector(0);
		planeNormals[1] = _simCell.cellNormalVector(1);
		planeNormals[2] = _simCell.cellNormalVector(2);

		Point3 corner1 = Point3::Origin() + _simCell.matrix().column(3);
		Point3 corner2 = corner1 + _simCell.matrix().column(0) + _simCell.matrix().column(1) + _simCell.matrix().column(2);

		QMutex bondMutex;
		QMutex indexMutex;

		// Perform analysis, particle-wise parallel.
		task()->setProgressMaximum(_positions->size());
		parallelForChunks(_positions->size(), *task(), [&](size_t startIndex, size_t chunkSize, Task& promise) {
			std::vector<int> localVoronoiBuffer;
			std::vector<size_t> localVoronoiBufferIndex;
			for(size_t index = startIndex; chunkSize--; index++) {
				if(promise.isCanceled()) return;
				if((index % 256) == 0) promise.incrementProgressValue(256);

				// Skip unselected particles (if requested).
				if(selectionArray && selectionArray[index] == 0)
					continue;

				// Build Voronoi cell.
				voro::voronoicell_neighbor v;

				// Initialize the Voronoi cell to be a cube larger than the simulation cell, centered at the origin.
				v.init(-boxDiameter, boxDiameter, -boxDiameter, boxDiameter, -boxDiameter, boxDiameter);

				// Cut Voronoi cell at simulation cell boundaries in non-periodic directions.
				bool skipParticle = false;
				for(size_t dim = 0; dim < 3; dim++) {
					if(!_simCell.pbcFlags()[dim]) {
						double r;
						r = 2 * planeNormals[dim].dot(corner2 - positionsArray[index]);
						if(r <= 0) skipParticle = true;
						v.nplane(planeNormals[dim].x() * r, planeNormals[dim].y() * r, planeNormals[dim].z() * r, r*r, -1);
						r = 2 * planeNormals[dim].dot(positionsArray[index] - corner1);
						if(r <= 0) skipParticle = true;
						v.nplane(-planeNormals[dim].x() * r, -planeNormals[dim].y() * r, -planeNormals[dim].z() * r, r*r, -1);
					}
				}
				// Skip particles that are located outside of non-periodic box boundaries.
				if(skipParticle)
					continue;

				// This function will be called for every neighbor particle.
				int nvisits = 0;
				auto visitFunc = [&](const NearestNeighborFinder::Neighbor& n, FloatType& mrs) {
					// Skip unselected particles (if requested).
					OVITO_ASSERT(!selectionArray || selectionArray[n.index]);
					FloatType rs = n.distanceSq;
					if(!_radii.empty())
						rs += _radii[index] - _radii[n.index];
					v.nplane(n.delta.x(), n.delta.y(), n.delta.z(), rs, n.index);
					if(nvisits == 0) {
						mrs = v.max_radius_squared();
						nvisits = 100;
					}
					nvisits--;
				};

				// Visit all neighbors of the current particles.
				nearestNeighborFinder.visitNeighbors(nearestNeighborFinder.particlePos(index), visitFunc);

				processCell(v, index, localVoronoiBuffer, localVoronoiBufferIndex, &bondMutex);
			}
			if(!localVoronoiBufferIndex.empty()) {
				QMutexLocker locker(&indexMutex);
				voronoiBufferIndex.insert(voronoiBufferIndex.end(), localVoronoiBufferIndex.cbegin(), localVoronoiBufferIndex.cend());
				voronoiBuffer.insert(voronoiBuffer.end(), localVoronoiBuffer.cbegin(), localVoronoiBuffer.cend());
			}
		});
	}

	if(maxFaceOrders()) {
		size_t componentCount = std::min(_maxFaceOrder.load(), FaceOrderStorageLimit);
		_voronoiIndices = std::make_shared<PropertyStorage>(_positions->size(), PropertyStorage::Int, componentCount, 0, QStringLiteral("Voronoi Index"), true);
		PropertyAccess<int,true> voronoiIndicesArray(_voronoiIndices);
		ConstPropertyAccess<int> maxFaceOrdersArray(maxFaceOrders());
		auto indexData = voronoiBuffer.cbegin();
		for(size_t particleIndex : voronoiBufferIndex) {
			size_t c = std::min(maxFaceOrdersArray[particleIndex], FaceOrderStorageLimit);
			for(size_t i = 0; i < c; i++) {
				voronoiIndicesArray.set(particleIndex, i, *indexData++);
			}
		}
		OVITO_ASSERT(indexData == voronoiBuffer.cend());

		// Re-use the output particle property as an output mesh region property.
		if(_computePolyhedra) {
			_polyhedraMesh.addRegionProperty(voronoiIndices());
			_polyhedraMesh.addRegionProperty(maxFaceOrders());
		}
	}

	// Finalize the polyhedral mesh.
	if(_computePolyhedra) {
		task()->nextProgressSubStep();
		task()->beginProgressSubStepsWithWeights({1,12,1,1,1});

		// First, connect adjacent faces from the same Voronoi cell.
		_polyhedraMesh.connectOppositeHalfedges();

		// The polyhedral cells should now be closed manifolds.
		OVITO_ASSERT(_polyhedraMesh.topology()->isClosed());
		task()->nextProgressSubStep();
		task()->setProgressMaximum(_polyhedraMesh.faceCount());

		// Merge mesh vertices that are shared by adjacent Voronoi polyhedra.

		// Initialize disjoint set data structure to keep track which vertices have been merged with which.
		std::vector<SurfaceMeshData::vertex_index> parents(_polyhedraMesh.vertexCount());
		std::vector<SurfaceMeshData::vertex_index> ranks(_polyhedraMesh.vertexCount(), 0);
		std::iota(parents.begin(), parents.end(), (SurfaceMeshData::vertex_index)0);

		// Iterate over all Voronoi faces.
		for(SurfaceMeshData::face_index face = 0; face < _polyhedraMesh.faceCount(); face++) {
			if(!task()->setProgressValueIntermittent(face)) return;
			SurfaceMeshData::region_index region = _polyhedraMesh.faceRegion(face);

			// We know for each Voronoi face which Voronoi polyhedron is on the other side. 
			SurfaceMeshData::region_index adjacentRegion = adjacentCellArray[face];
			// Skip faces that are at the outer surface.
			if(adjacentRegion < 0) continue;
			// Skip faces that belong to a periodic polyhedron.	
			if(adjacentRegion == region) continue;

			// Iterate over all vertices of the current Voronoi face.
			SurfaceMeshData::edge_index ffe = _polyhedraMesh.firstFaceEdge(face);
			SurfaceMeshData::edge_index edge = ffe;
			do {
				// Get the coordinates of the current vertex.
				SurfaceMeshData::vertex_index vertex = _polyhedraMesh.vertex2(edge);
				const Point3& vertex_pos = _polyhedraMesh.vertexPosition(vertex);

				// Iterate over all vertices of the adjacent Voronoi cell.
				FloatType longest_dist = 0;
				FloatType shortest_dist = std::numeric_limits<FloatType>::max();
				SurfaceMeshData::vertex_index closest_vertex = HalfEdgeMesh::InvalidIndex;
				for(SurfaceMeshData::vertex_index other_vertex = polyhedraVertices[adjacentRegion].first, end_vertex = other_vertex + polyhedraVertices[adjacentRegion].second; other_vertex != end_vertex; ++other_vertex) {

					// Check if vertex has an adjacent face leading back to the current Voronoi cell.
					bool isCandidateVertex = false;
					for(SurfaceMeshData::edge_index adj_edge = _polyhedraMesh.firstVertexEdge(other_vertex); adj_edge != HalfEdgeMesh::InvalidIndex; adj_edge = _polyhedraMesh.nextVertexEdge(adj_edge)) {
						if(adjacentCellArray[_polyhedraMesh.adjacentFace(adj_edge)] == region) {
							isCandidateVertex = true;
							break;
						}
					}
					if(!isCandidateVertex) continue;

					// Compute distance of other vertex to current vertex.
					FloatType squared_dist = _polyhedraMesh.cell().wrapVector(_polyhedraMesh.vertexPosition(other_vertex) - vertex_pos).squaredLength();

					// Determine the closest vertex and longest distance (as a measure of the cell size).
					if(squared_dist > longest_dist) longest_dist = squared_dist;
					if(squared_dist < shortest_dist) {
						shortest_dist = squared_dist;
						closest_vertex = other_vertex;
					}
				}
				OVITO_ASSERT(closest_vertex != HalfEdgeMesh::InvalidIndex);

				// Determine a threshold distance for testing whether the two vertices should be merged.
				FloatType distance_threshold = sqrt(longest_dist) * FloatType(1e-9);
				if(shortest_dist <= distance_threshold) {
					// Merge the two vertices.

					// Find root and make root as parent (path compression)
					SurfaceMeshData::vertex_index parentA = parents[vertex];
					while(parentA != parents[parentA]) {
						parentA = parents[parentA];
					}
					parents[vertex] = parentA;
					SurfaceMeshData::vertex_index parentB = parents[closest_vertex];
					while(parentB != parents[parentB]) {
						parentB = parents[parentB];
					}
					parents[closest_vertex] = parentB;
					if(parentA != parentB) {

						// Attach smaller rank tree under root of high rank tree (Union by Rank)
						if(ranks[parentA] < ranks[parentB]) {
							parents[parentA] = parentB;
						}
						else {
							parents[parentB] = parentA;

							// If ranks are same, then make one as root and increment its rank by one.
							if(ranks[parentA] == ranks[parentB])
								ranks[parentA]++;
						}
					}
				}

				edge = _polyhedraMesh.nextFaceEdge(edge);
			}
			while(edge != ffe);
		}
		task()->nextProgressSubStep();

		// Transfer edges from vertices that are going to be deleted to ramining vertices.
		for(SurfaceMeshData::edge_index edge = 0; edge < _polyhedraMesh.edgeCount(); edge++) {
			SurfaceMeshData::vertex_index new_vertex = parents[_polyhedraMesh.vertex2(edge)];
			_polyhedraMesh.topology()->transferFaceBoundaryToVertex(edge, new_vertex);
			if(task()->isCanceled()) return;
		}
		task()->nextProgressSubStep();

		// Delete unused vertices.
		for(SurfaceMeshData::vertex_index vertex = _polyhedraMesh.vertexCount() - 1; vertex >= 0; vertex--) {
			if(parents[vertex] != vertex) {
				_polyhedraMesh.deleteVertex(vertex);
				if(task()->isCanceled()) return;
			}
		}
		task()->nextProgressSubStep();
		task()->setProgressMaximum(_polyhedraMesh.faceCount());

		// Connect pairs of internal Voronoi faces.
		for(SurfaceMeshData::face_index face = 0; face < _polyhedraMesh.faceCount(); face++) {
			if(_polyhedraMesh.hasOppositeFace(face)) continue;
			if(!task()->setProgressValueIntermittent(face)) return;

			// We know for each Voronoi face which Voronoi polyhedron is on the other side.
			SurfaceMeshData::region_index adjacentRegion = adjacentCellArray[face];
			// Skip faces that belong to the outer surface.
			if(adjacentRegion < 0) continue;
			// Periodic polyhedra pose a problem.	
			if(adjacentRegion == _polyhedraMesh.faceRegion(face)) {
				throw Exception(tr("Cannot generate polyhedron mesh for this input, because at least one Voronoi cell is touching a periodic image of itself. To avoid this error you can try to use the Replicate modifier or turn off periodic boundary conditions for the simulation cell."));
			}

			SurfaceMeshData::edge_index first_edge = _polyhedraMesh.firstFaceEdge(face);
			SurfaceMeshData::vertex_index vertex1 = _polyhedraMesh.vertex1(first_edge);
			SurfaceMeshData::vertex_index vertex2 = _polyhedraMesh.vertex2(first_edge);

			// Iterate over all edges/faces adjacent to one of the vertices.
			for(SurfaceMeshData::edge_index edge = _polyhedraMesh.firstVertexEdge(vertex1); edge != HalfEdgeMesh::InvalidIndex; edge = _polyhedraMesh.nextVertexEdge(edge)) {
				SurfaceMeshData::face_index adjacentFace = _polyhedraMesh.adjacentFace(edge);
				if(_polyhedraMesh.faceRegion(adjacentFace) != adjacentRegion) continue;
				SurfaceMeshData::edge_index oppositeEdge = _polyhedraMesh.findEdge(adjacentFace, vertex2, vertex1);
				if(oppositeEdge != HalfEdgeMesh::InvalidIndex) {
					OVITO_ASSERT(!_polyhedraMesh.hasOppositeFace(adjacentFace));
					_polyhedraMesh.topology()->linkOppositeFaces(face, adjacentFace);
					break;
				}
			}
			OVITO_ASSERT(_polyhedraMesh.hasOppositeFace(face));
		}

		// Remove the "Adjacent Cell" property from the mesh faces, because the user is typically not interested in it.
		_polyhedraMesh.removeFaceProperty(_adjacentCellProperty);

		task()->endProgressSubSteps();
	}

	task()->endProgressSubSteps();
}

/******************************************************************************
* Injects the computed results of the engine into the data pipeline.
******************************************************************************/
void VoronoiAnalysisModifier::VoronoiAnalysisEngine::emitResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	VoronoiAnalysisModifier* modifier = static_object_cast<VoronoiAnalysisModifier>(modApp->modifier());
	ParticlesObject* particles = state.expectMutableObject<ParticlesObject>();

	if(_inputFingerprint.hasChanged(particles))
		modApp->throwException(tr("Cached modifier results are obsolete, because the number or the storage order of input particles has changed."));

	particles->createProperty(coordinationNumbers());
	particles->createProperty(atomicVolumes());

	if(modifier->computeIndices()) {
		if(voronoiIndices())
			particles->createProperty(voronoiIndices());
		if(maxFaceOrders())
			particles->createProperty(maxFaceOrders());

		state.setStatus(PipelineStatus(PipelineStatus::Success,
			tr("Maximum face order: %1").arg(maxFaceOrder().load())));
	}

	// Check computed Voronoi cell volume sum.
	FloatType simulationBoxVolume = _simCell.volume3D();
	if(particles->elementCount() != 0 && std::abs(voronoiVolumeSum() - simulationBoxVolume) > 1e-8 * particles->elementCount() * simulationBoxVolume) {
		state.setStatus(PipelineStatus(PipelineStatus::Warning,
				tr("The volume sum of all Voronoi cells does not match the simulation box volume. "
						"This may be a result of particles being located outside of the simulation box boundaries. "
						"See user manual for more information.\n"
						"Simulation box volume: %1\n"
						"Voronoi cell volume sum: %2").arg(simulationBoxVolume).arg(voronoiVolumeSum())));
	}

	if(modifier->computeBonds()) {
		// Insert output object into the pipeline.
		particles->addBonds(bonds(), modifier->bondsVis());
	}
	
	if(modifier->computePolyhedra()) {
		// Output the surface mesh representing the computed Voronoi polyhedra.
		SurfaceMesh* meshObj = state.createObject<SurfaceMesh>(QStringLiteral("voronoi-polyhedra"), modApp, tr("Voronoi polyhedra"));
		polyhedraMesh().transferTo(meshObj);
		meshObj->setDomain(state.getObject<SimulationCellObject>());
		meshObj->setVisElement(modifier->polyhedraVis());
	}

	state.addAttribute(QStringLiteral("Voronoi.max_face_order"), QVariant::fromValue(maxFaceOrder().load()), modApp);
}

}	// End of namespace
}	// End of namespace
