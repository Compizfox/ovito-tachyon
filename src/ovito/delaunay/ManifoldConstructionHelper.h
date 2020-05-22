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
#include <ovito/stdobj/simcell/SimulationCell.h>
#include <ovito/stdobj/properties/PropertyStorage.h>
#include <ovito/mesh/surface/SurfaceMeshData.h>
#include <ovito/core/utilities/concurrent/Task.h>
#include <ovito/delaunay/DelaunayTessellation.h>

#if 0
#include <boost/functional/hash.hpp>
#endif

namespace Ovito { namespace Delaunay {

using namespace Ovito::Mesh;

/**
 * Constructs a SurfaceMesh structure from a DelaunayTessellation representing the separating surface manifold
 * between different spatial regions of the tessellation.
 */
template<bool CreateTwoSidedMesh, bool CreateRegions>
class ManifoldConstructionHelper
{
public:

	// A no-op face-preparation functor.
	struct DefaultPrepareMeshFaceFunc {
		void operator()(HalfEdgeMesh::face_index face,
				const std::array<size_t,3>& vertexIndices,
				const std::array<DelaunayTessellation::VertexHandle,3>& vertexHandles,
				DelaunayTessellation::CellHandle cell) {}
	};

	// A no-op vertex-preparation functor.
	struct DefaultPrepareMeshVertexFunc {
		void operator()(HalfEdgeMesh::vertex_index vertex, size_t particleIndex) {}
	};

public:

	/// Constructor.
	ManifoldConstructionHelper(DelaunayTessellation& tessellation, SurfaceMeshData& outputMesh, FloatType alpha,
			const PropertyStorage& positions) : _tessellation(tessellation), _mesh(outputMesh), _alpha(alpha), _positions(positions) {}

	/// Returns the number of filled regions that have been identified.
	SurfaceMeshData::size_type filledRegionCount() const { return _filledRegionCount; }
	
	/// Returns the number of empty regions that have been identified.
	SurfaceMeshData::size_type emptyRegionCount() const { return _emptyRegionCount; }

	/// This is the main function, which constructs the manifold triangle mesh.
	template<typename CellRegionFunc, typename PrepareMeshFaceFunc = DefaultPrepareMeshFaceFunc, typename PrepareMeshVertexFunc = DefaultPrepareMeshVertexFunc>
	bool construct(CellRegionFunc&& determineCellRegion, Task& task,
			PrepareMeshFaceFunc&& prepareMeshFaceFunc = PrepareMeshFaceFunc(), PrepareMeshVertexFunc&& prepareMeshVertexFunc = PrepareMeshVertexFunc())
	{
		// Algorithm is divided into several sub-steps.
		if(CreateRegions)
			task.beginProgressSubStepsWithWeights({1,8,2,1});
		else
			task.beginProgressSubStepsWithWeights({1,1,2});

		// Assign tetrahedra to spatial regions.
		if(!classifyTetrahedra(std::move(determineCellRegion), task))
			return false;

		task.nextProgressSubStep();

		// Aggregate connected tetrahedra into regions.
		if(CreateRegions) {

			// Create the "Region" face property in the output mesh.
			_mesh.createFaceProperty(SurfaceMeshFaces::RegionProperty);

			if(!formFilledRegions(task))
				return false;
			task.nextProgressSubStep();
		}

		// Create triangle facets at interfaces between two different regions.
		if(!createInterfaceFacets(std::move(prepareMeshFaceFunc), std::move(prepareMeshVertexFunc), task))
			return false;

		task.nextProgressSubStep();

		// Connect triangles with one another to form a closed manifold.
		if(!linkHalfedges(task))
			return false;

#ifdef OVITO_DEBUG
		// Verify that generated manifold connectivity is correct.
		if(CreateTwoSidedMesh) {
			for(SurfaceMeshData::edge_index edge = 0; edge < _mesh.edgeCount(); edge++) {
				// Each edge must be part of at least two opposite manifolds.
				OVITO_ASSERT(_mesh.countManifolds(edge) >= 2);
				// Manifold count of the edge must be consistent with its opposite edge.
				OVITO_ASSERT(_mesh.countManifolds(edge) == _mesh.countManifolds(_mesh.oppositeEdge(edge)));
			}
		}
#endif

		task.endProgressSubSteps();

		return !task.isCanceled();
	}

	/// Finds disconnected exterior regions and computes their volumes.
	bool formEmptyRegions(Task& task) {
		OVITO_ASSERT(CreateTwoSidedMesh && CreateRegions);
		task.beginProgressSubStepsWithWeights({1,1});

		// Count number of identified regions.
		_emptyRegionCount = 0;

		// Flags indicating which periodic cell directions are connected by a surface through the cell boundary.
		const SimulationCell simCell = _tessellation.simCell();
		bool surfaceCrossesBoundaries[3] = { false, false, false };
		bool detectBoundaryCrossings = simCell.hasPbc();

		// Stack of faces to visit. Used for implementing recursive algorithm.
		std::deque<SurfaceMeshData::face_index> facesToProcess;

		// Identify disconnected components of the surface mesh bordering to an empty region.
		task.setProgressMaximum(_mesh.faceCount() / 2); // Note: Dividing by two, because only every other face is oriented towards the empty region.
		for(SurfaceMeshData::face_index face = 0; face < _mesh.faceCount(); face++) {
			// Look for mesh faces that are not adjacent to a filled region and which have not been visited yet.
			SurfaceMeshData::region_index faceRegion = _mesh.faceRegion(face);
			if(faceRegion != HalfEdgeMesh::InvalidIndex)
				continue;

			// Found a first seed face. Start a new mesh component.
			faceRegion = _filledRegionCount + (_emptyRegionCount++);
			_mesh.setFaceRegion(face, faceRegion);
			size_t faceCount = 1;

			facesToProcess.push_back(face);
			do {
				// Take next face from the stack.
				SurfaceMeshData::face_index currentFace = facesToProcess.front();
				facesToProcess.pop_front();
				if(!task.incrementProgressValue())
					return false;

				// Visit neighbors of current face.
				SurfaceMeshData::edge_index firstEdge = _mesh.firstFaceEdge(currentFace);
				SurfaceMeshData::edge_index edge = firstEdge;
				OVITO_ASSERT(firstEdge != HalfEdgeMesh::InvalidIndex);
				do {
					// Determine whether this edge crosses a periodic simulation cell boundary.
					if(detectBoundaryCrossings) {
						const Point3& p1 = _mesh.vertexPosition(_mesh.vertex1(edge));
						const Point3& p2 = _mesh.vertexPosition(_mesh.vertex2(edge));
						Vector3 delta = p2 - p1;
						for(size_t dim = 0; dim < 3; dim++) {
							if(!surfaceCrossesBoundaries[dim] && simCell.hasPbc(dim)) {
								// The edge is crossing the periodic boundary if it spans more than half of the simulation cell in that direction.
								if(std::abs(simCell.inverseMatrix().prodrow(delta, dim)) >= FloatType(0.5)) {
									surfaceCrossesBoundaries[dim] = true;
									detectBoundaryCrossings = (simCell.hasPbc(0) && !surfaceCrossesBoundaries[0]) || (simCell.hasPbc(1) && !surfaceCrossesBoundaries[1]) || (simCell.hasPbc(2) && !surfaceCrossesBoundaries[2]);
								}
							}
						}
					}

					// Get mesh face adjacent to the opposite half-edge of the current half-edge. 
					SurfaceMeshData::edge_index oppositeEdge = _mesh.oppositeEdge(edge);
					OVITO_ASSERT(oppositeEdge != HalfEdgeMesh::InvalidIndex);
					SurfaceMeshData::face_index neighborFace = _mesh.adjacentFace(oppositeEdge);
					OVITO_ASSERT(neighborFace != HalfEdgeMesh::InvalidIndex);

					if(_mesh.faceRegion(neighborFace) == HalfEdgeMesh::InvalidIndex) {
						// Assign neighbor face to the current empty region.
						_mesh.setFaceRegion(neighborFace, faceRegion);
						// Put on recursive stack.
						facesToProcess.push_back(neighborFace);
						faceCount++;
					}

					edge = _mesh.nextFaceEdge(edge);
				}
				while(edge != firstEdge);
			}
			while(!facesToProcess.empty());
		}

		task.nextProgressSubStep();
		std::deque<DelaunayTessellation::CellHandle> cellsToProcess;

		// Disjoint set data structure for merging empty regions:
		std::vector<size_t> regionParents(_emptyRegionCount);
		std::vector<size_t> regionSizes(_emptyRegionCount);
		std::vector<double> regionVolumes(_emptyRegionCount, 0.0);
		std::iota(regionParents.begin(), regionParents.end(), (size_t)0);
		std::fill(regionSizes.begin(), regionSizes.end(), 1);

		// "Find" part of Union-Find:
		auto findRegion = [&](size_t index) {
			OVITO_ASSERT(index < regionParents.size());
			// Find root and make root as parent of i (path compression).
			size_t parent = regionParents[index];
			while(parent != regionParents[parent]) {
				parent = regionParents[parent];
			}
			regionParents[index] = parent;
			return parent;
		};
		// "Union" part of Union-Find.
		auto mergeRegions = [&](size_t regionA, size_t regionB) {
			if(regionA == regionB) return;
			size_t parentA = findRegion(regionA - _filledRegionCount);
			size_t parentB = findRegion(regionB - _filledRegionCount);
			if(parentA != parentB) {
				// Attach smaller tree under root of larger tree.
				if(regionSizes[parentA] < regionSizes[parentB]) {
					regionParents[parentA] = parentB;
					regionSizes[parentB] += regionSizes[parentA];
					regionVolumes[parentB] += regionVolumes[parentA];
				}
				else {
					regionParents[parentB] = parentA;
					regionSizes[parentA] += regionSizes[parentB];
					regionVolumes[parentA] += regionVolumes[parentB];
				}
			}
		};

		// The ID of the empty region that is split by a periodic simulation box boundary.
		SurfaceMeshData::region_index splitPeriodicRegion = HalfEdgeMesh::InvalidIndex;

		// Loop over all cells to cluster them.
		task.setProgressMaximum(_tessellation.numberOfTetrahedra() - _numFilledCells);
		for(DelaunayTessellation::CellIterator cellIter = _tessellation.begin_cells(); cellIter != _tessellation.end_cells() && !task.isCanceled(); ++cellIter) {
			DelaunayTessellation::CellHandle cell = *cellIter;

			// Only consider finite cells.
			if(!_tessellation.isFiniteCell(cell))
				continue;

			// Skip filled cells, or cells that have been visited before.
			if(_tessellation.getUserField(cell) != HalfEdgeMesh::InvalidIndex)
				continue;

			// Determine whether the cell is adjacent to a face of the generated surface mesh.
			SurfaceMeshData::region_index emptyRegion = HalfEdgeMesh::InvalidIndex;
			for(int f = 0; f < 4; f++) {

				SurfaceMeshData::face_index adjacentMeshFace = findCellFace(_tessellation.mirrorFacet(cell, f));
				if(adjacentMeshFace == HalfEdgeMesh::InvalidIndex)
					continue;
				SurfaceMeshData::face_index oppositeFace = _mesh.oppositeFace(adjacentMeshFace);

				OVITO_ASSERT(oppositeFace != HalfEdgeMesh::InvalidIndex);
				OVITO_ASSERT(emptyRegion == HalfEdgeMesh::InvalidIndex || emptyRegion == _mesh.faceRegion(oppositeFace));
				// The surface mesh face should be bordering an empty and not a filled region.
				if(_mesh.faceRegion(oppositeFace) >= _filledRegionCount)
					emptyRegion = _mesh.faceRegion(oppositeFace);
			}

			// Skip the cell if it is not adjacent to a face of the surface mesh.
			if(emptyRegion == HalfEdgeMesh::InvalidIndex)
				continue;

			// Start recursive algorithm to gather all connected empty tetrahedra.
			cellsToProcess.push_back(cell);
			_tessellation.setUserField(cell, emptyRegion);
			double regionVolume = 0;
			do {
				// Take next tetrahedron from the stack.
				DelaunayTessellation::CellHandle currentCell = cellsToProcess.front();
				cellsToProcess.pop_front();
				if(!task.incrementProgressValue())
					return false;

				// Flags indicating whether the current tetrahedron crosses the cell boundaries.
				std::array<bool,3> cellCrossesBoundaries = {false, false, false};

				// Compute the overlap of the tetrahedron with the simulation box.
				FloatType overlapVolume = calculateVolumeOverlap(currentCell, cellCrossesBoundaries);
				// Stop at cells that are completely outside of the simulation box.
				if(overlapVolume == 0)
					continue;
				regionVolume += overlapVolume;

				// Determine whether the current spatial region is split by a periodic box boundary
				// that is not crossed by the surface mesh. In this case, we need to merge the two empty
				// regions that were created separately on either side of the simulation box boundary. 
				for(size_t dim = 0; dim < 3; dim++) {
					if(cellCrossesBoundaries[dim] && !surfaceCrossesBoundaries[dim] && simCell.hasPbc(dim)) {
						if(splitPeriodicRegion == HalfEdgeMesh::InvalidIndex)
							splitPeriodicRegion = emptyRegion;
						else
							mergeRegions(emptyRegion, splitPeriodicRegion);
					}
				}

				// Continue with neighboring cells.
				for(int f = 0; f < 4; f++) {

					// Check if this is a region border or not.
					DelaunayTessellation::Facet mirrorFacet = _tessellation.mirrorFacet(currentCell, f);
					SurfaceMeshData::face_index adjacentMeshFace = findCellFace(mirrorFacet);
					if(adjacentMeshFace != HalfEdgeMesh::InvalidIndex) {
						// We have hit a region border. Merge the two empty region IDs.
						SurfaceMeshData::face_index oppositeFace = _mesh.oppositeFace(adjacentMeshFace);
						SurfaceMeshData::region_index secondEmptyRegion = _mesh.faceRegion(oppositeFace);
						if(secondEmptyRegion >= _filledRegionCount) {
							mergeRegions(emptyRegion, secondEmptyRegion);
						}
					}
					else {
						DelaunayTessellation::CellHandle neighborCell = mirrorFacet.first;

						// Only consider finite cells.
						if(!_tessellation.isFiniteCell(neighborCell))
							continue;

						// Skip filled cells, or cells that have been visited before.
						SurfaceMeshData::region_index neighborRegion = _tessellation.getUserField(neighborCell);
						if(neighborRegion != HalfEdgeMesh::InvalidIndex) {
							if(neighborRegion >= _filledRegionCount)
								mergeRegions(emptyRegion, neighborRegion);
							continue;
						}

						// Put neighbor cell on processing stack.
						_tessellation.setUserField(neighborCell, emptyRegion);
						cellsToProcess.push_back(neighborCell);
					}
				}
			}
			while(!cellsToProcess.empty());

			// Add the total volume of the visited tetrahedra to the region's volume.
			regionVolumes[findRegion(emptyRegion - _filledRegionCount)] += regionVolume;
		}
		if(task.isCanceled())
			return false;

		// Remap merged regions to contiguous range.
		std::vector<SurfaceMeshData::region_index> regionMapping(_emptyRegionCount);
		_emptyRegionCount = 0;
		for(size_t i = 0; i < regionParents.size(); i++) {
			if(findRegion(i) == i) {
				regionMapping[i] = _emptyRegionCount + _filledRegionCount;
				_mesh.createRegion(0, regionVolumes[i]);
				_emptyRegionCount++;
			}
		}

		// Map each mesh face's preliminary region ID to the final ID.
		for(SurfaceMeshData::region_index& region : _mesh.faceRegionsRange()) {
			if(region >= _filledRegionCount) {
				region = regionMapping[findRegion(region - _filledRegionCount)];
				OVITO_ASSERT(region < _mesh.regionCount());
			}
		}

		// Create one space-filling empty region if there is no filled region.
		if(_emptyRegionCount == 0 && _filledRegionCount == 0) {
			_mesh.createRegion(0, simCell.volume3D());
			_emptyRegionCount = 1;
		}

		task.endProgressSubSteps();
		return !task.isCanceled();
	}

private:

	/// Assigns each tetrahedron to a region.
	template<typename CellRegionFunc>
	bool classifyTetrahedra(CellRegionFunc&& determineCellRegion, Task& task)
	{
		task.setProgressValue(0);
		task.setProgressMaximum(_tessellation.numberOfTetrahedra());

		_numFilledCells = 0;
		size_t progressCounter = 0;
		_mesh.setSpaceFillingRegion(HalfEdgeMesh::InvalidIndex);
		bool spaceFillingRegionUndetermined = true;
		bool isSpaceFilling = true;
		for(DelaunayTessellation::CellIterator cellIter = _tessellation.begin_cells(); cellIter != _tessellation.end_cells(); ++cellIter) {
			DelaunayTessellation::CellHandle cell = *cellIter;

			// Update progress indicator.
			if(!task.setProgressValueIntermittent(progressCounter++))
				return false;

			// Alpha-shape criterion: This determines whether the Delaunay tetrahedron is part of the filled region.
			bool isFilledTetrehedron = _tessellation.isFiniteCell(cell) && _tessellation.alphaTest(cell, _alpha);

			SurfaceMeshData::region_index region = HalfEdgeMesh::InvalidIndex;
			if(isFilledTetrehedron) {
				region = determineCellRegion(cell);
				OVITO_ASSERT(region >= 0 || region == HalfEdgeMesh::InvalidIndex);
				OVITO_ASSERT(!CreateRegions || (region == 0 || region == HalfEdgeMesh::InvalidIndex));
				OVITO_ASSERT(CreateRegions || region < _mesh.regionCount());
			}
			_tessellation.setUserField(cell, region);

			if(!_tessellation.isGhostCell(cell)) {
				if(spaceFillingRegionUndetermined) {
					_mesh.setSpaceFillingRegion(region);
					spaceFillingRegionUndetermined = false;
				}
				else {
					if(isSpaceFilling && _mesh.spaceFillingRegion() != region) {
						_mesh.setSpaceFillingRegion(HalfEdgeMesh::InvalidIndex);
						isSpaceFilling = false;
					}
				}
			}

			if(region != HalfEdgeMesh::InvalidIndex && !_tessellation.isGhostCell(cell)) {
				_tessellation.setCellIndex(cell, _numFilledCells++);
			}
			else {
				_tessellation.setCellIndex(cell, -1);
			}
		}

		return !task.isCanceled();
	}

	/// Computes the volume of a Delaunay tetrahedron.
	FloatType cellVolume(DelaunayTessellation::CellHandle cell) const
	{
		Point3 p0 = _tessellation.vertexPosition(_tessellation.cellVertex(cell, 0));
		Vector3 ad = _tessellation.vertexPosition(_tessellation.cellVertex(cell, 1)) - p0;
		Vector3 bd = _tessellation.vertexPosition(_tessellation.cellVertex(cell, 2)) - p0;
		Vector3 cd = _tessellation.vertexPosition(_tessellation.cellVertex(cell, 3)) - p0;
		return std::abs(ad.dot(cd.cross(bd))) / FloatType(6);
	}

	/// Aggregates Delaunay tetrahedra into connected regions.
	bool formFilledRegions(Task& task)
	{
		task.beginProgressSubStepsWithWeights({2,3,1});

		// Create a lookup map that allows retrieving the primary image of a Delaunay cell for a triangular face formed by three particles.
		if(!createCellMap(task))
			return false;

		// Make sure no regions have been defined so far.
		OVITO_ASSERT(_mesh.regionCount() == 0);

		// Create the output property arrays for the identified regions.
 		_mesh.createRegionProperty(SurfaceMeshRegions::VolumeProperty, true);

		std::deque<size_t> toProcess;

		// Loop over all cells to cluster them.
		task.nextProgressSubStep();
		task.setProgressMaximum(_tessellation.numberOfTetrahedra());
		for(DelaunayTessellation::CellIterator cellIter = _tessellation.begin_cells(); cellIter != _tessellation.end_cells() && !task.isCanceled(); ++cellIter) {
			DelaunayTessellation::CellHandle cell = *cellIter;

			// Skip empty cells and cells that have already been assigned to a region.
			if(_tessellation.getUserField(cell) != 0)
				continue;

			// Skip ghost cells.
			if(_tessellation.isGhostCell(cell))
				continue;

			// Start a new region.
			int currentRegionId = _mesh.regionCount() + 1;
			OVITO_ASSERT(currentRegionId >= 1);
			_tessellation.setUserField(cell, currentRegionId);
			// Compute total region volume as we go while visiting all tetrahedra.
			double regionVolume = 0;

			// Now recursively iterate over all neighbors of the seed cell and add them to the current region.
			toProcess.push_back(cell);
			do {
				DelaunayTessellation::CellHandle currentCell = toProcess.front();
				toProcess.pop_front();
				if(!task.incrementProgressValue())
					return false;

				// Add the volume of the current cell to the total region volume.
				regionVolume += cellVolume(currentCell);

				// Loop over the 4 facets of the cell
				for(int f = 0; f < 4; f++) {

					// Get the 3 vertices of the facet.
					// Note that we reverse their order to find the opposite face.
					std::array<size_t, 3> vertices;
					for(int v = 0; v < 3; v++)
						vertices[v] = _tessellation.vertexIndex(_tessellation.cellVertex(currentCell, DelaunayTessellation::cellFacetVertexIndex(f, 2-v)));

					// Bring vertices into a well-defined order, which can be used as lookup key to find the adjacent tetrahedron.
					reorderFaceVertices(vertices);

					// Look up the adjacent Delaunay cell.
					auto neighborCell = _cellLookupMap.find(vertices);
					if(neighborCell != _cellLookupMap.end()) {
						// Add adjacent cell to the deque if it has not been visited yet.
						if(_tessellation.getUserField(neighborCell->second) == 0) {
							toProcess.push_back(neighborCell->second);
							_tessellation.setUserField(neighborCell->second, currentRegionId);
						}
					}
				}
			}
			while(!toProcess.empty());

			// Create a spatial region in the output mesh.
			SurfaceMeshData::region_index regionId = _mesh.createRegion(0, regionVolume);
			OVITO_ASSERT(regionId + 1 == currentRegionId);
		}
		task.nextProgressSubStep();

		_filledRegionCount = _mesh.regionCount();
		if(_mesh.regionCount() > 0) {
			// Shift filled region IDs to start at index 0.
			for(DelaunayTessellation::CellIterator cellIter = _tessellation.begin_cells(); cellIter != _tessellation.end_cells(); ++cellIter) {
				int region = _tessellation.getUserField(*cellIter);
				if(region > 0)
					_tessellation.setUserField(*cellIter, region - 1);
			}

			// Copy assigned region IDs from primary tetrahedra to ghost tetrahedra.
			task.setProgressMaximum(_tessellation.numberOfTetrahedra());
			for(DelaunayTessellation::CellIterator cellIter = _tessellation.begin_cells(); cellIter != _tessellation.end_cells(); ++cellIter) {
				DelaunayTessellation::CellHandle cell = *cellIter;
				if(_tessellation.isGhostCell(cell) && _tessellation.getUserField(cell) == 0) {
					if(!task.setProgressValueIntermittent(cell))
						return false;

					// Get the 3 vertices of the first face of the tet.
					std::array<size_t, 3> vertices;
					for(int v = 0; v < 3; v++)
						vertices[v] = _tessellation.vertexIndex(_tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(0, v)));

					// Bring vertices into a well-defined order, which can be used as lookup key.
					reorderFaceVertices(vertices);

					// Find the primary tet whose face connects the same three particles.
					auto neighborCell = _cellLookupMap.find(vertices);
					if(neighborCell != _cellLookupMap.end()) {
						_tessellation.setUserField(cell, _tessellation.getUserField(neighborCell->second));
					}
				}
			}
		}
		task.endProgressSubSteps();

		return true;
	}

	/// Creates a lookup map that allows to retrieve the primary Delaunay cell image that belongs to a
	/// triangular face formed by three particles.
	bool createCellMap(Task& task)
	{
		task.setProgressMaximum(_tessellation.numberOfTetrahedra());
		for(DelaunayTessellation::CellIterator cellIter = _tessellation.begin_cells(); cellIter != _tessellation.end_cells(); ++cellIter) {
			DelaunayTessellation::CellHandle cell = *cellIter;

			// Skip cells that belong to the exterior region.
			if(_tessellation.getUserField(cell) == HalfEdgeMesh::InvalidIndex)
				continue;

			// Skip ghost cells.
			if(_tessellation.isGhostCell(cell))
				continue;

			// Update progress indicator.
			if(!task.setProgressValueIntermittent(cell))
				return false;

			// Loop over the 4 facets of the cell.
			for(int f = 0; f < 4; f++) {

				// Get the 3 vertices of the facet.
				std::array<size_t, 3> vertices;
				for(int v = 0; v < 3; v++)
					vertices[v] = _tessellation.vertexIndex(_tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(f, v)));

				// Bring vertices into a well-defined order, which can be used as lookup key.
				reorderFaceVertices(vertices);

				OVITO_ASSERT(_cellLookupMap.find(vertices) == _cellLookupMap.end());

				// Add facet and its adjacent cell to the loopup map.
				_cellLookupMap.emplace(vertices, cell);
			}
		}
		return !task.isCanceled();
	}

	/// Constructs the triangle facets that separate different regions in the tetrahedral mesh.
	template<typename PrepareMeshFaceFunc, typename PrepareMeshVertexFunc>
	bool createInterfaceFacets(PrepareMeshFaceFunc&& prepareMeshFaceFunc, PrepareMeshVertexFunc&& prepareMeshVertexFunc, Task& task)
	{
		// Stores the triangle mesh vertices created for the vertices of the tetrahedral mesh.
		std::vector<SurfaceMeshData::vertex_index> vertexMap(_positions.size(), HalfEdgeMesh::InvalidIndex);
		_tetrahedraFaceList.clear();
		_faceLookupMap.clear();

		task.setProgressValue(0);
		task.setProgressMaximum(_numFilledCells);

		for(DelaunayTessellation::CellIterator cellIter = _tessellation.begin_cells(); cellIter != _tessellation.end_cells(); ++cellIter) {
			DelaunayTessellation::CellHandle cell = *cellIter;

			// Consider only filled local tetrahedra.
			if(_tessellation.getCellIndex(cell) == -1) continue;
			SurfaceMeshData::region_index filledRegion = _tessellation.getUserField(cell);
			OVITO_ASSERT(filledRegion != HalfEdgeMesh::InvalidIndex);

			// Update progress indicator.
			if(!task.setProgressValueIntermittent(_tessellation.getCellIndex(cell)))
				return false;

			Point3 unwrappedVerts[4];
			for(int i = 0; i < 4; i++)
				unwrappedVerts[i] = _tessellation.vertexPosition(_tessellation.cellVertex(cell, i));

			// Check validity of tessellation. 
			// Delaunay edges (of filled tetrahedro) should never span more than half of the simulation box in periodic directions.
			Vector3 ad = unwrappedVerts[0] - unwrappedVerts[3];
			Vector3 bd = unwrappedVerts[1] - unwrappedVerts[3];
			Vector3 cd = unwrappedVerts[2] - unwrappedVerts[3];
			if(_tessellation.simCell().isWrappedVector(ad) || _tessellation.simCell().isWrappedVector(bd) || _tessellation.simCell().isWrappedVector(cd))
				throw Exception("Cannot construct manifold. Simulation cell length is too small for the given probe sphere radius parameter.");

			// Iterate over the four faces of the tetrahedron cell.
			_tessellation.setCellIndex(cell, -1);
			for(int f = 0; f < 4; f++) {

				// Check if the adjacent tetrahedron belongs to a different region.
				std::pair<DelaunayTessellation::CellHandle,int> mirrorFacet = _tessellation.mirrorFacet(cell, f);
				DelaunayTessellation::CellHandle adjacentCell = mirrorFacet.first;
				if(_tessellation.getUserField(adjacentCell) == filledRegion)
					continue;

				// Create the three vertices of the face or use existing output vertices.
				std::array<SurfaceMeshData::vertex_index,3> facetVertices;
				std::array<DelaunayTessellation::VertexHandle,3> vertexHandles;
				std::array<size_t,3> vertexIndices;
				for(int v = 0; v < 3; v++) {
					vertexHandles[v] = _tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(f, _flipOrientation ? v : (2-v)));
					size_t vertexIndex = vertexIndices[v] = _tessellation.vertexIndex(vertexHandles[v]);
					OVITO_ASSERT(vertexIndex < vertexMap.size());
					if(vertexMap[vertexIndex] == HalfEdgeMesh::InvalidIndex) {
						vertexMap[vertexIndex] = _mesh.createVertex(_positions[vertexIndex]);
						prepareMeshVertexFunc(vertexMap[vertexIndex], vertexIndex);
					}
					facetVertices[v] = vertexMap[vertexIndex];
				}

				// Create a new triangle facet.
				SurfaceMeshData::face_index face = _mesh.createFace(facetVertices.begin(), facetVertices.end(), filledRegion);

				// Tell client code about the new facet.
				prepareMeshFaceFunc(face, vertexIndices, vertexHandles, cell);

				// Create additional face for exterior region if requested.
				if(CreateTwoSidedMesh && _tessellation.getUserField(adjacentCell) == HalfEdgeMesh::InvalidIndex) {

					// Build face vertex list.
					std::reverse(std::begin(vertexHandles), std::end(vertexHandles));
					std::array<size_t,3> reverseVertexIndices;
					for(int v = 0; v < 3; v++) {
						vertexHandles[v] = _tessellation.cellVertex(adjacentCell, DelaunayTessellation::cellFacetVertexIndex(mirrorFacet.second, _flipOrientation ? v : (2-v)));
						size_t vertexIndex = reverseVertexIndices[v] = _tessellation.vertexIndex(vertexHandles[v]);
						OVITO_ASSERT(vertexIndex < vertexMap.size());
						OVITO_ASSERT(vertexMap[vertexIndex] != HalfEdgeMesh::InvalidIndex);
						facetVertices[v] = vertexMap[vertexIndex];
					}

					// Create a new triangle facet.
					HalfEdgeMesh::face_index oppositeFace = _mesh.createFace(facetVertices.begin(), facetVertices.end(), HalfEdgeMesh::InvalidIndex);

					// Tell client code about the new facet.
					prepareMeshFaceFunc(oppositeFace, reverseVertexIndices, vertexHandles, adjacentCell);

					// Insert new facet into lookup map.
					reorderFaceVertices(reverseVertexIndices);
					_faceLookupMap.emplace(reverseVertexIndices, oppositeFace);
				}

				// Insert new facet into lookup map.
				reorderFaceVertices(vertexIndices);
				_faceLookupMap.emplace(vertexIndices, face);

				// Insert into contiguous list of tetrahedron faces.
				if(_tessellation.getCellIndex(cell) == -1) {
					_tessellation.setCellIndex(cell, _tetrahedraFaceList.size());
					_tetrahedraFaceList.push_back(std::array<HalfEdgeMesh::face_index, 4>{{ HalfEdgeMesh::InvalidIndex, HalfEdgeMesh::InvalidIndex, HalfEdgeMesh::InvalidIndex, HalfEdgeMesh::InvalidIndex }});
				}
				_tetrahedraFaceList[_tessellation.getCellIndex(cell)][f] = face;
			}
		}

		return !task.isCanceled();
	}

	SurfaceMeshData::face_index findAdjacentFace(DelaunayTessellation::CellHandle cell, int f, int e)
	{
		int vertexIndex1, vertexIndex2;
		if(!_flipOrientation) {
			vertexIndex1 = DelaunayTessellation::cellFacetVertexIndex(f, 2-e);
			vertexIndex2 = DelaunayTessellation::cellFacetVertexIndex(f, (4-e)%3);
		}
		else {
			vertexIndex1 = DelaunayTessellation::cellFacetVertexIndex(f, (e+1)%3);
			vertexIndex2 = DelaunayTessellation::cellFacetVertexIndex(f, e);
		}
		DelaunayTessellation::FacetCirculator circulator_start = _tessellation.incident_facets(cell, vertexIndex1, vertexIndex2, cell, f);
		DelaunayTessellation::FacetCirculator circulator = circulator_start;
		OVITO_ASSERT((*circulator).first == cell);
		OVITO_ASSERT((*circulator).second == f);
		--circulator;
		OVITO_ASSERT(circulator != circulator_start);
		int region = _tessellation.getUserField(cell);
		do {
			// Look for the first cell while going around the edge that belongs to a different region.
			if(_tessellation.getUserField((*circulator).first) != region)
				break;
			--circulator;
		}
		while(circulator != circulator_start);
		OVITO_ASSERT(circulator != circulator_start);

		// Get the current adjacent cell, which is part of the same region as the first tet.
		std::pair<DelaunayTessellation::CellHandle,int> mirrorFacet = _tessellation.mirrorFacet(*circulator);
		OVITO_ASSERT(_tessellation.getUserField(mirrorFacet.first) == region);

		SurfaceMeshData::face_index adjacentFace = findCellFace(mirrorFacet);
		if(adjacentFace == HalfEdgeMesh::InvalidIndex)
			throw Exception("Cannot construct mesh for this input dataset. Adjacent cell face not found.");
		return adjacentFace;
	}

	bool linkHalfedges(Task& task)
	{
		task.setProgressValue(0);
		task.setProgressMaximum(_tetrahedraFaceList.size());

		auto tet = _tetrahedraFaceList.cbegin();
		for(DelaunayTessellation::CellIterator cellIter = _tessellation.begin_cells(); cellIter != _tessellation.end_cells(); ++cellIter) {
			DelaunayTessellation::CellHandle cell = *cellIter;

			// Look for filled cells being adjacent to at least one mesh face.
			if(_tessellation.getCellIndex(cell) == -1) continue;
			OVITO_ASSERT(_tetrahedraFaceList.cbegin() + _tessellation.getCellIndex(cell) == tet);

			// Update progress indicator.
			if(!task.setProgressValueIntermittent(_tessellation.getCellIndex(cell)))
				return false;

			// Visit the mesh faces adjacent to the current cell.
			for(int f = 0; f < 4; f++) {
				SurfaceMeshData::face_index facet = (*tet)[f];
				if(facet == HalfEdgeMesh::InvalidIndex) continue;

				// Link within manifold.
				SurfaceMeshData::edge_index edge = _mesh.firstFaceEdge(facet);
				for(int e = 0; e < 3; e++, edge = _mesh.nextFaceEdge(edge)) {
					if(_mesh.hasOppositeEdge(edge)) continue;
					SurfaceMeshData::face_index adjacentFace = findAdjacentFace(cell, f, e);
					SurfaceMeshData::edge_index oppositeEdge = _mesh.findEdge(adjacentFace, _mesh.vertex2(edge), _mesh.vertex1(edge));
					if(oppositeEdge == HalfEdgeMesh::InvalidIndex)
						throw Exception("Cannot construct mesh for this input dataset. Opposite half-edge not found.");
					_mesh.linkOppositeEdges(edge, oppositeEdge);
				}

				if(CreateTwoSidedMesh) {
					std::pair<DelaunayTessellation::CellHandle,int> oppositeFacet = _tessellation.mirrorFacet(cell, f);
					OVITO_ASSERT(_tessellation.getUserField(oppositeFacet.first) != _tessellation.getUserField(cell));
					SurfaceMeshData::face_index outerFacet = findCellFace(oppositeFacet);
					OVITO_ASSERT(outerFacet != HalfEdgeMesh::InvalidIndex);

					// Link opposite mesh faces (Note: they may have already been linked).
					_mesh.topology()->linkOppositeFaces(facet, outerFacet);

					// Link adjacent facets in opposite manifold if it is bounding an empty region.
					if(_tessellation.getUserField(oppositeFacet.first) == HalfEdgeMesh::InvalidIndex) {
						SurfaceMeshData::edge_index edge = _mesh.firstFaceEdge(outerFacet);
						for(int e = 0; e < 3; e++, edge = _mesh.nextFaceEdge(edge)) {
							if(_mesh.hasOppositeEdge(edge)) continue;
							SurfaceMeshData::face_index adjacentFace = findAdjacentFace(oppositeFacet.first, oppositeFacet.second, e);
							SurfaceMeshData::edge_index oppositeEdge = _mesh.findEdge(adjacentFace, _mesh.vertex2(edge), _mesh.vertex1(edge));
							if(oppositeEdge == HalfEdgeMesh::InvalidIndex)
								throw Exception("Cannot construct mesh for this input dataset. Opposite half-edge (2) not found.");
							_mesh.linkOppositeEdges(edge, oppositeEdge);
						}
					}
				}
			}
			++tet;
		}
		OVITO_ASSERT(tet == _tetrahedraFaceList.cend());
		OVITO_ASSERT(_mesh.topology()->isClosed());

		// Set up manifold pointers at edges of the mesh.
		if(CreateTwoSidedMesh) {
			for(SurfaceMeshData::edge_index edge1 = 0; edge1 < _mesh.edgeCount() && !task.isCanceled(); edge1++) {
				// Link surface manifolds.
				SurfaceMeshData::edge_index oppositeEdge = _mesh.oppositeEdge(edge1);
				SurfaceMeshData::face_index adjacentFace = _mesh.adjacentFace(oppositeEdge);
				SurfaceMeshData::face_index oppositeFace = _mesh.oppositeFace(adjacentFace);
				for(SurfaceMeshData::edge_index edge2 = _mesh.firstFaceEdge(oppositeFace); ; edge2 = _mesh.nextFaceEdge(edge2)) {
					if(_mesh.vertex2(edge2) == _mesh.vertex2(edge1)) {
						_mesh.setNextManifoldEdge(edge1, edge2);
						break;
					}
				}
			}
		}

		return !task.isCanceled();
	}

	SurfaceMeshData::face_index findCellFace(const std::pair<DelaunayTessellation::CellHandle,int>& facet)
	{
		// If the cell is a ghost cell, find the corresponding real cell first.
		auto cell = facet.first;
		if(_tessellation.getCellIndex(cell) != -1) {
			OVITO_ASSERT(_tessellation.getCellIndex(cell) >= 0 && _tessellation.getCellIndex(cell) < _tetrahedraFaceList.size());
			return _tetrahedraFaceList[_tessellation.getCellIndex(cell)][facet.second];
		}
		else {
			std::array<size_t,3> faceVerts;
			for(size_t i = 0; i < 3; i++) {
				int vertexIndex = DelaunayTessellation::cellFacetVertexIndex(facet.second, _flipOrientation ? i : (2-i));
				faceVerts[i] = _tessellation.vertexIndex(_tessellation.cellVertex(cell, vertexIndex));
			}
			reorderFaceVertices(faceVerts);
			auto iter = _faceLookupMap.find(faceVerts);
			if(iter != _faceLookupMap.end())
				return iter->second;
			else
				return HalfEdgeMesh::InvalidIndex;
		}
	}

	static void reorderFaceVertices(std::array<size_t,3>& vertexIndices) {
#if !defined(Q_OS_MACOS) && !defined(Q_OS_WASM)
		// Shift the order of vertices so that the smallest index is at the front.
		std::rotate(std::begin(vertexIndices), std::min_element(std::begin(vertexIndices), std::end(vertexIndices)), std::end(vertexIndices));
#else
		// Workaround for compiler bug in Xcode 10.0. Clang hangs when compiling the code above with -O2/-O3 flag.
		auto min_index = std::min_element(vertexIndices.begin(), vertexIndices.end()) - vertexIndices.begin();
		std::rotate(vertexIndices.begin(), vertexIndices.begin() + min_index, vertexIndices.end());
#endif
	}

	/// Computes the volume of the given Delaunay cell that is overlapping with the simulation box.
	FloatType calculateVolumeOverlap(DelaunayTessellation::CellHandle cell, std::array<bool,3>& outsideDir) 
	{
		// Gather the positions of the four Delaunay vertices and check if they are all inside the simulation cell domain.
		bool isCompletelyInsideBox = true;
		Point3 vertexPositions[4];
		Point3 reducedVertexPositions[4];
		for(int v = 0; v < 4; v++) {
			const Point3& vpos = vertexPositions[v] = _tessellation.vertexPosition(_tessellation.cellVertex(cell, v));
			const Point3& rp = reducedVertexPositions[v] = _tessellation.simCell().absoluteToReduced(vpos);
			if(rp.x() < 0.0 || rp.x() > 1.0) {
				isCompletelyInsideBox = false;
				outsideDir[0] = true;
			}
			if(rp.y() < 0.0 || rp.y() > 1.0) {
				isCompletelyInsideBox = false;
				outsideDir[1] = true;
			}
			if(rp.z() < 0.0 || rp.z() > 1.0) {
				isCompletelyInsideBox = false;
				outsideDir[2] = true;
			}
		}

		// Compute volume of the full tetrahedron if it is completely inside the simulation box.
		if(isCompletelyInsideBox) {
			Vector3 ad = vertexPositions[1] - vertexPositions[0];
			Vector3 bd = vertexPositions[2] - vertexPositions[0];
			Vector3 cd = vertexPositions[3] - vertexPositions[0];
			return std::abs(ad.dot(cd.cross(bd))) / FloatType(6);
		}

		// There may be a partial overlap. We need to compute the intersection of the two shapes to determine the overlap volume.

		// We start with six line segments along the six edges of the tetrahedron.
		// The six line segments will be clipped at the boundaries of the simulation box.
		Point3 lineSegments[(6 + 12) * 2];
		static const int edgeVertices[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
		int numPoints = 0;
		for(int e = 0; e < 6; e++) {
			Point3 p1 = reducedVertexPositions[edgeVertices[e][0]];
			Point3 p2 = reducedVertexPositions[edgeVertices[e][1]];

			// Clip the segment at the boundaries of the simulation box.
			bool isDegenerate = false;
			for(size_t dim = 0; dim < 3; dim++) {
				// Clip at lower box boundary.
				if(p1[dim] <= 0 && p2[dim] <= 0) {
					isDegenerate = true;
					break;
				}
				else if(p1[dim] < 0 && p2[dim] >= 0) {
					Vector3 delta = p2 - p1;
					p1 -= delta * (p1[dim] / delta[dim]);
					OVITO_ASSERT(std::abs(p1[dim]) < FLOATTYPE_EPSILON);
				}
				else if(p1[dim] >= 0 && p2[dim] < 0) {
					Vector3 delta = p2 - p1;
					p2 -= delta * (p2[dim] / delta[dim]);
					OVITO_ASSERT(std::abs(p2[dim]) < FLOATTYPE_EPSILON);
				}

				// Clip at upper box boundary.
				if(p1[dim] >= 1.0 && p2[dim] >= 1.0) {
					isDegenerate = true;
					break;
				}
				else if(p1[dim] > 1.0 && p2[dim] <= 1.0) {
					Vector3 delta = p2 - p1;
					p1 -= delta * ((p1[dim] - 1.0) / delta[dim]);
					OVITO_ASSERT(std::abs(p1[dim] - 1.0) < FLOATTYPE_EPSILON);
				}
				else if(p1[dim] <= 1.0 && p2[dim] > 1.0) {
					Vector3 delta = p2 - p1;
					p2 -= delta * ((p2[dim] - 1.0) / delta[dim]);
					OVITO_ASSERT(std::abs(p2[dim] - 1.0) < FLOATTYPE_EPSILON);
				}
			}

			// Add the line segment to the list of output segments if it wasn't clipped away entirely.
			if(!isDegenerate) {
				lineSegments[numPoints++] = p1;
				lineSegments[numPoints++] = p2;
			}
		}

		// Compute the four planes of the tetrahedron faces.
		Plane3 planes[4];
		for(int f = 0; f < 4; f++) {
			const Point3& p0 = reducedVertexPositions[DelaunayTessellation::cellFacetVertexIndex(f, 0)];
			const Point3& p1 = reducedVertexPositions[DelaunayTessellation::cellFacetVertexIndex(f, 1)];
			const Point3& p2 = reducedVertexPositions[DelaunayTessellation::cellFacetVertexIndex(f, 2)];
			planes[f] = Plane3(p0, p1, p2, true);
		}

		// Clip the 12 edges of the unit box at the 4 four planes of the tetrahedron.
		static const Point3 boxVertices[12][2] = {
			{Point3(0,0,0),Point3(1,0,0)}, {Point3(0,1,0),Point3(1,1,0)}, {Point3(0,0,1),Point3(1,0,1)}, {Point3(0,1,1),Point3(1,1,1)},
			{Point3(0,0,0),Point3(0,1,0)}, {Point3(1,0,0),Point3(1,1,0)}, {Point3(0,0,1),Point3(0,1,1)}, {Point3(1,0,1),Point3(1,1,1)},
			{Point3(0,0,0),Point3(0,0,1)}, {Point3(1,0,0),Point3(1,0,1)}, {Point3(0,1,0),Point3(0,1,1)}, {Point3(1,1,0),Point3(1,1,1)}
		};
		for(int e = 0; e < 12; e++) {
			Point3 p1 = boxVertices[e][0];
			Point3 p2 = boxVertices[e][1];

			// Clip at each plane.
			bool isDegenerate = false;
			for(const Plane3& plane : planes) {
				FloatType t1 = plane.pointDistance(p1);
				FloatType t2 = plane.pointDistance(p2);
				if(t1 < 0 && t2 < 0) {
					isDegenerate = true;
					break;
				}
				else if(t1 < 0 && t2 > 0) {
					Vector3 delta = p2 - p1;
					p1 -= delta * (t1 / (t2 - t1));
					OVITO_ASSERT(std::abs(plane.pointDistance(p1)) <= FLOATTYPE_EPSILON);
				}
				else if(t1 > 0 && t2 < 0) {
					Vector3 delta = p2 - p1;
					p2 += delta * (t2 / (t1 - t2));
					OVITO_ASSERT(std::abs(plane.pointDistance(p2)) <= FLOATTYPE_EPSILON);
				}
			}

			// Add the line segment to the list of output segments if it wasn't clipped away entirely.
			if(!isDegenerate) {
				lineSegments[numPoints++] = p1;
				lineSegments[numPoints++] = p2;
			}
		}
		if(numPoints < 4)
			return 0;	// Overlap is degenerate.

		// Construct convex hull of remaining line segments.
		SurfaceMeshData mesh;
		mesh.constructConvexHull(std::vector<Point3>(lineSegments, lineSegments + numPoints));

		// The convex hull may be empty if the input point set is degenerate.
		if(mesh.faceCount() == 0)
			return 0;

		// Compute volume enclosed by the convex hull polyhedron.
		const Point3 apex = mesh.vertexPosition(0);
		FloatType convexVolume = 0;
		for(SurfaceMeshData::edge_index firstEdge : mesh.firstFaceEdges()) {
			SurfaceMeshData::edge_index edge1 = mesh.nextFaceEdge(firstEdge);
			SurfaceMeshData::edge_index edge2 = mesh.nextFaceEdge(edge1);
			Matrix3 tripod;
			tripod.column(0) = mesh.vertexPosition(mesh.vertex1(firstEdge)) - apex;
			tripod.column(2) = mesh.vertexPosition(mesh.vertex2(firstEdge)) - apex;
			while(edge2 != firstEdge) {
				tripod.column(1) = tripod.column(2);
				tripod.column(2) = mesh.vertexPosition(mesh.vertex2(edge1)) - apex;
				convexVolume += tripod.determinant();
				edge1 = edge2;
				edge2 = mesh.nextFaceEdge(edge2);
			}
		}

		return (convexVolume / 6.0) * _tessellation.simCell().volume3D();
	}

private:

	/// The tetrahedral tessellation.
	DelaunayTessellation& _tessellation;

	/// The squared probe sphere radius used to classify tetrahedra as open or solid.
	FloatType _alpha;

	/// Counts the number of tetrehedral cells that belong to the filled regions.
	size_t _numFilledCells = 0;

	/// Number of filled regions that have been identified.
	SurfaceMeshData::size_type _filledRegionCount = 0;
	
	/// Number of empty regions that have been identified.
	SurfaceMeshData::size_type _emptyRegionCount = 0;

	/// The input particle positions.
	ConstPropertyAccess<Point3> _positions;

	/// The output mesh topology.
	SurfaceMeshData& _mesh;

	/// Controls the reversal of the normal orientation of the generated surface facets.
	bool _flipOrientation = false;

	/// Stores the faces of the local tetrahedra that have a least one facet for which a triangle has been created.
	std::vector<std::array<SurfaceMeshData::face_index, 4>> _tetrahedraFaceList;

	/// This map allows looking up surface mesh faces based on their three vertices.
#if 0
	std::unordered_map<std::array<size_t,3>, SurfaceMeshData::face_index, boost::hash<std::array<size_t,3>>> _faceLookupMap;
#else
	std::map<std::array<size_t,3>, SurfaceMeshData::face_index> _faceLookupMap;
#endif

	/// This map allows looking up the tetrahedron that is adjacent to a given triangular face.
#if 0
	std::unordered_map<std::array<size_t,3>, DelaunayTessellation::CellHandle, boost::hash<std::array<size_t,3>>> _cellLookupMap;
#else
	std::map<std::array<size_t,3>, DelaunayTessellation::CellHandle> _cellLookupMap;
#endif
};

}	// End of namespace
}	// End of namespace
