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

#include <ovito/grid/Grid.h>
#include <ovito/grid/objects/VoxelGrid.h>
#include <ovito/grid/modifier/CreateIsosurfaceModifier.h>
#include <ovito/mesh/surface/SurfaceMesh.h>
#include <ovito/mesh/surface/SurfaceMeshVis.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include <ovito/core/dataset/DataSet.h>
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include "VoxelGridSliceModifierDelegate.h"
#include "MarchingCubes.h"

namespace Ovito { namespace Grid {

IMPLEMENT_OVITO_CLASS(VoxelGridSliceModifierDelegate);
DEFINE_REFERENCE_FIELD(VoxelGridSliceModifierDelegate, surfaceMeshVis);

/******************************************************************************
* Constructs the object.
******************************************************************************/
VoxelGridSliceModifierDelegate::VoxelGridSliceModifierDelegate(DataSet* dataset) : SliceModifierDelegate(dataset)
{
	// Create the vis element for rendering the mesh.
	setSurfaceMeshVis(new SurfaceMeshVis(dataset));
	surfaceMeshVis()->setShowCap(false);
	surfaceMeshVis()->setHighlightEdges(false);
	surfaceMeshVis()->setSmoothShading(false);
	surfaceMeshVis()->setSurfaceIsClosed(false);
	surfaceMeshVis()->setObjectTitle(tr("Volume slice"));
}

/******************************************************************************
* Calculates a cross-section of a voxel grid.
******************************************************************************/
PipelineStatus VoxelGridSliceModifierDelegate::apply(Modifier* modifier, PipelineFlowState& state, TimePoint time, ModifierApplication* modApp, const std::vector<std::reference_wrapper<const PipelineFlowState>>& additionalInputs)
{
	OVITO_ASSERT(!dataset()->undoStack().isRecording());

	SliceModifier* mod = static_object_cast<SliceModifier>(modifier);
	QString statusMessage;

	// Obtain modifier parameter values.
	Plane3 plane;
	FloatType sliceWidth;
	std::tie(plane, sliceWidth) = mod->slicingPlane(time, state.mutableStateValidity());
	sliceWidth /= 2;
	bool invert = mod->inverse();
	int numPlanes = 0;
	Plane3 planes[2];
	if(sliceWidth <= 0) {
		planes[numPlanes++] = plane;
	}
	else {
		planes[numPlanes++] = Plane3( plane.normal,  plane.dist + sliceWidth);
		planes[numPlanes++] = Plane3(-plane.normal, -plane.dist + sliceWidth);
	}

	for(const DataObject* obj : state.data()->objects()) {
		if(const VoxelGrid* voxelGrid = dynamic_object_cast<VoxelGrid>(obj)) {
			// Get domain of voxel grid.
			VoxelGrid::GridDimensions gridShape = voxelGrid->shape();
			SimulationCell cell = voxelGrid->domain()->data();
			if(cell.is2D())
				continue;

			// Verify consistency of input property container.
			voxelGrid->verifyIntegrity();

			// Construct cross section mesh using a special version of the marching cubes algorithm.
			SurfaceMeshData mesh(cell);

			// The slice plane does NOT exist in a periodic domain.
			mesh.cell().setPbcFlags(false, false, false);

			// The level of subdivision.
			const int resolution = 2;

			Plane3 planeGridSpace;
			for(int pidx = 0; pidx < numPlanes; pidx++) {
				// Transform plane from world space to orthogonal grid space.
				planeGridSpace = (Matrix3(
					gridShape[0]*resolution, 0, 0,
					0, gridShape[1]*resolution, 0,
					0, 0, gridShape[2]*resolution) * cell.inverseMatrix()) * planes[pidx];

				// Set up callback function returning the field value, which will be passed to the marching cubes algorithm.
				auto getFieldValue = [&](int i, int j, int k) {
					return planeGridSpace.pointDistance(Point3(i,j,k));
				};

				MarchingCubes mc(mesh, gridShape[0]*resolution, gridShape[1]*resolution, gridShape[2]*resolution, false, std::move(getFieldValue), true);
				mc.generateIsosurface(0.0, *SynchronousOperation::createSignal(dataset()->taskManager()).task());
			}

			// Create a manifold by connecting adjacent faces.
			mesh.connectOppositeHalfedges();

			// Form quadrilaterals from pairs of triangles.
			// This only makes sense when the slicing plane is aligned with the grid cells axes such that only quads result from 
			// the marching cubes algorithm.
			if(std::abs(planeGridSpace.normal.x()) <= FLOATTYPE_MAX || std::abs(planeGridSpace.normal.y()) <= FLOATTYPE_MAX || std::abs(planeGridSpace.normal.z()) <= FLOATTYPE_MAX)
				mesh.makeQuadrilateralFaces();

			// Transform from double resolution grid to single resolution.
			mesh.transformVertices(AffineTransformation(
				FloatType(1) / resolution, 0, 0, -0.5 + FloatType(1) / resolution,
				0, FloatType(1) / resolution, 0, -0.5 + FloatType(1) / resolution,
				0, 0, FloatType(1) / resolution, -0.5 + FloatType(1) / resolution));

			// Collect the set of voxel grid properties that should be transferred over to the isosurface mesh vertices.
			std::vector<ConstPropertyPtr> fieldProperties;
			for(const PropertyObject* property : voxelGrid->properties())
				fieldProperties.push_back(property->storage());

			// Copy field values from voxel grid to surface mesh vertices.
			SynchronousOperation operation = SynchronousOperation::createSignal(dataset()->taskManager());
			CreateIsosurfaceModifier::transferPropertiesFromGridToMesh(*operation.task(), mesh, fieldProperties, cell, gridShape);

			// Transform mesh vertices from orthogonal grid space to world space.
			const AffineTransformation tm = cell.matrix() * Matrix3(
				FloatType(1) / gridShape[0], 0, 0,
				0, FloatType(1) / gridShape[1], 0,
				0, 0, FloatType(1) / gridShape[2]) *
				AffineTransformation::translation(Vector3(0.5, 0.5, 0.5));
			mesh.transformVertices(tm);

			// Flip surface orientation if cell matrix is a mirror transformation.
			if(tm.determinant() < 0)
				mesh.flipFaces();
			
			// Create the output mesh data object.
			SurfaceMesh* meshObj = state.createObject<SurfaceMesh>(QStringLiteral("volume-slice"), modApp, tr("Volume slice"));
			mesh.transferTo(meshObj);
			meshObj->setVisElement(surfaceMeshVis());
		}
	}

	return PipelineStatus(PipelineStatus::Success, statusMessage);
}

}	// End of namespace
}	// End of namespace
