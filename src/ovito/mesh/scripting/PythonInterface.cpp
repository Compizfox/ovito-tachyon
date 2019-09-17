///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2018) Alexander Stukowski
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

#include <ovito/mesh/Mesh.h>
#include <ovito/mesh/io/VTKFileImporter.h>
#include <ovito/mesh/io/VTKTriangleMeshExporter.h>
#include <ovito/mesh/io/WavefrontOBJImporter.h>
#include <ovito/mesh/io/STLImporter.h>
#include <ovito/mesh/tri/TriMeshObject.h>
#include <ovito/mesh/tri/TriMeshVis.h>
#include <ovito/mesh/surface/SurfaceMesh.h>
#include <ovito/mesh/surface/SurfaceMeshVis.h>
#include <ovito/pyscript/binding/PythonBinding.h>
#include <ovito/core/app/PluginManager.h>
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include <ovito/core/utilities/io/CompressedTextWriter.h>

namespace Ovito { namespace Mesh {

using namespace PyScript;

PYBIND11_MODULE(MeshPython, m)
{
	// Register the classes of this plugin with the global PluginManager.
	PluginManager::instance().registerLoadedPluginClasses();

	py::options options;
	options.disable_function_signatures();

	ovito_class<TriMeshObject, DataObject>{m}
	;

	ovito_class<TriMeshVis, DataVis>(m, nullptr,
		// Python class name:
		"TriMeshVis")
		.def_property("color", &TriMeshVis::color, &TriMeshVis::setColor)
		.def_property("transparency", &TriMeshVis::transparency, &TriMeshVis::setTransparency)
		.def_property("highlight_edges", &TriMeshVis::highlightEdges, &TriMeshVis::setHighlightEdges,
				"Activates the highlighted rendering of the polygonal edges of the mesh."
				"\n\n"
				":Default: ``False``\n")
	;

	ovito_class<SurfaceMeshVertices, PropertyContainer>{m};
	ovito_class<SurfaceMeshFaces, PropertyContainer>{m};
	ovito_class<SurfaceMeshRegions, PropertyContainer>{m};

	auto SurfaceMesh_py = ovito_class<SurfaceMesh, PeriodicDomainDataObject>(m,
			":Base class: :py:class:`ovito.data.DataObject`"
			"\n\n"
			"This data object type stores a triangle mesh describing a surface or, more precisely, a two-dimensional manifold that is closed and orientable. "
			"Typically, surface meshes are produced by modifiers such as the :py:class:`~ovito.modifiers.ConstructSurfaceModifier`, "
			":py:class:`~ovito.modifiers.CreateIsosurfaceModifier` or :py:class:`~ovito.modifiers.CoordinationPolyhedraModifier`. "
			"See also the corresponding :ovitoman:`page of the user manual <../../scene_objects.surface_mesh>` for more information on surface meshes."
			"\n\n"
			"**Periodic domains**"
			"\n\n"
			"What is particular about surface meshes is that they may be embedded in a periodic domain, i.e. in a simulation cell with periodic boundary conditions applied. "
			"That means triangles of a surface mesh can connect vertices on opposite sides of a simulation box and wrap around correctly. "
			"OVITO takes care of computing the intersections of the periodic surface with the box boundaries and automatically produces a non-periodic representation of the triangle mesh "
			"when it comes to visualizing the surface. "
			"\n\n"
			"The spatial domain the surface mesh is embedded in is represented by a :py:class:`~ovito.data.SimulationCell` object, which is attached to the "
			":py:class:`!SurfaceMesh` instance. You can access it through the :py:attr:`.domain` attribute. "
			"\n\n"
			"**Visual representation**"
			"\n\n"
			"The visual appearance of the surface mesh in rendered images is controlled by an attached :py:class:`~ovito.vis.SurfaceMeshVis` element, which is "
			"accessible through the :py:attr:`~DataObject.vis` base class attribute. "
			"\n\n"
			"**Spatial regions**"
			"\n\n"
			"As surface meshes are closed orientable manifolds, one can define an *interior* and an *exterior* region of space that are separated by the manifold. "
			"For example, if the surface mesh is constructed by the :py:class:`~ovito.modifiers.ConstructSurfaceModifier` from a set of particles, "
			"then the region enclosed by the surface is the \"solid\" region and the outside region is the one containing no particles. "
			"\n\n"
			"It can be that there is no interior region and the exterior region is infinite and fills all space. In this case the surface mesh is degenerate and "
			"comprises no triangles. The opposite extreme is also possible in periodic domains: The interior region extends over the entire domain "
			"and there is no outside region. Again, the surface mesh will consist of zero triangles in this case. "
			"To discriminate between the two situations, the :py:class:`!SurfaceMesh` class provides the :py:attr:`.space_filling_region` field, which is "
			"set when the interior region fills the entire periodic domain. "
			"\n\n"
			"The :py:meth:`locate_point` method can be used to test whether some point in space belongs to the interior or the exterior region. "
			"\n\n"
			"**File export**"
			"\n\n"
			"A surface mesh can be written to a file in the form of a conventional triangle mesh. "
			"For this, a non-periodic version is produced by truncating triangles at the domain boundaries and generating \"cap polygons\" to fill the holes that "
			"occur at the intersection of the interior region with the domain boundaries. To export the mesh, use the :py:func:`ovito.io.export_file` function "
			"and select ``vtk/trimesh`` as output format: "
			"\n\n"
			".. literalinclude:: ../example_snippets/surface_mesh_export.py\n"
			"   :lines: 7-\n"
			"\n\n"
			"**Cutting planes**"
			"\n\n"
			"A set of *cutting planes* can be assigned to a :py:class:`!SurfaceMesh` to cut away parts of the mesh for visualization purposes. "
			"This may be useful to e.g. cut a hole into a closed surface allowing to look inside the enclosed volume. "
			"The :py:class:`!SurfaceMesh` objects manages a list of cutting planes, which are accessible through the :py:meth:`.get_cutting_planes` and :py:meth:`.set_cutting_planes` "
			"methods. Note that the cuts are non-destructive and get performed only on a transient copy of the mesh generated during image rendering or when exporting the mesh to a file. "
			"The original data structures of the mesh are not affected. "
			"The :py:class:`~ovito.modifiers.SliceModifier`, which can act on a :py:class:`!SurfaceMesh`, performs the slice by simply adding a new entry to the :py:class:`!SurfaceMesh`'s "
			"list of cutting planes. "
			"\n\n"
			"**Mesh data access**"
			"\n\n"
			"The methods :py:meth:`.get_vertices`, :py:meth:`.get_faces` and :py:meth:`.get_face_adjacency` methods provide access to the internal data of the "
			"surface mesh. "
		)
		.def("locate_point", &SurfaceMesh::locatePoint,
			"locate_point(pos, eps=1e-6)"
			"\n\n"
			"Determines the index of the spatial region that contains the given location in 3-D space. "
			"Note that region index 0 is typically reserved for the empty/outer region, which doesn't contain any atoms or particles. "
			"Regions starting at index 1 are the filled/solid regions. "
			"\n\n"
			"The parameter *eps* is used as a precision thedholds to detect cases where the query point is positioned exactly on the surface itself, i.e. on the boundary "
			"between two spatial regions. Such a condition is indicates by the special return value -1. You can set *eps* to 0.0 to disable "
			"the point-on-boundary test. Then the method will never yield -1 as a result. "
			"\n\n"
			":param pos: The (x,y,z) coordinates of the query point\n"
			":param eps: Numerical precision threshold for point-on-boundary test\n"
			":return: The ID of the spatial region containing *pos*; or -1 if *pos* is exactly on the dividing surface between two regions\n",
			py::arg("pos"), py::arg("eps") = 1e-6)

		.def("get_vertices", [](const SurfaceMesh& meshObj) {
				meshObj.verifyMeshIntegrity();
				size_t vcount = meshObj.vertices()->elementCount();
				ConstPropertyPtr positions = meshObj.vertices()->expectProperty(SurfaceMeshVertices::PositionProperty)->storage();
				py::array_t<FloatType> array({ vcount, (size_t)3 });
				auto r = array.mutable_unchecked();
				for(size_t i = 0; i < vcount; i++) {
					for(size_t j = 0; j < 3; j++)
						r(i,j) = positions->getFloatComponent(i, j);
				}
				return array;
			},
			"Returns a *N* x 3 array with the xyz coordinates of the *N* vertices in the mesh. "
			"Note that the returned Numpy array is a copy of the internal data stored by the :py:class:`!SurfaceMesh`. ")

		.def("get_faces", [](const SurfaceMesh& meshObj) {
				meshObj.verifyMeshIntegrity();
				HalfEdgeMeshPtr topology = meshObj.topology();
				size_t fcount = topology->faceCount();
				py::array_t<HalfEdgeMesh::vertex_index> array({ fcount, (size_t)3 });
				auto r = array.mutable_unchecked();
				for(size_t i = 0; i < fcount; i++) {
					if(topology->countFaceEdges(i) != 3)
						meshObj.throwException("Mesh contains at least one face that is not a triangle.");
					r(i,0) = topology->firstFaceVertex(i);
					r(i,1) = topology->secondFaceVertex(i);
					r(i,2) = topology->thirdFaceVertex(i);
				}
				return array;
			},
			"Returns a *M* x 3 array with the vertex indices of the *M* triangles in the mesh. "
			"Note that the returned Numpy array is a copy of the internal data stored by the :py:class:`!SurfaceMesh`. "
			"Also keep in mind that triangle faces can cross the domain boundaries if the periodic boundary conditions are used. ")

		.def("get_face_adjacency", [](const SurfaceMesh& meshObj) {
				meshObj.verifyMeshIntegrity();
				HalfEdgeMeshPtr topology = meshObj.topology();
				size_t fcount = topology->faceCount();
				py::array_t<HalfEdgeMesh::face_index> array({ fcount, (size_t)3 });
				auto r = array.mutable_unchecked();
				for(size_t i = 0; i < fcount; i++) {
					if(topology->countFaceEdges(i) != 3)
						meshObj.throwException("Mesh contains at least one face that is not a triangle.");
					auto edge = topology->firstFaceEdge(i);
					for(size_t j = 0; j < 3; j++, edge = topology->nextFaceEdge(edge)) {
						if(!topology->hasOppositeEdge(edge))
							meshObj.throwException("Mesh is not closed. Some faces do not have an adjacent face.");
						r(i,j) = topology->adjacentFace(topology->oppositeEdge(edge));
					}
				}
				return array;
			},
			"Returns a *M* x 3 array listing the indices of the three faces that are adjacent to each of the *M* triangle faces in the mesh. "
			"This information can be used to traverse the neighbors of triangle faces. Every triangle face has exactly three neighbors, because surface "
			"meshes are closed manifolds. ")

		.def("get_cutting_planes", [](const SurfaceMesh& meshObj) {
				py::array_t<FloatType> array({ (size_t)meshObj.cuttingPlanes().size(), (size_t)4 });
				auto r = array.mutable_unchecked();
				for(size_t i = 0; i < meshObj.cuttingPlanes().size(); i++) {
					r(i,0) = meshObj.cuttingPlanes()[i].normal.x();
					r(i,1) = meshObj.cuttingPlanes()[i].normal.y();
					r(i,2) = meshObj.cuttingPlanes()[i].normal.z();
					r(i,3) = meshObj.cuttingPlanes()[i].dist;
				}
				return array;
			},
			"Returns a *N* x 4 array containing the definitions of the *N* cutting planes attached to this :py:class:`!SurfaceMesh`. "
			"\n\n"
			"Each plane is defined by its unit normal vector and a signed displacement magnitude, which determines the plane's distance from the coordinate origin along the normal, "
			"giving four numbers per plane in total. Those parts of the surface mesh which are on the positive side of the plane (in the direction the normal vector) are cut away. "
			"\n\n"
			"Note that the returned Numpy array is a copy of the internal data stored by the :py:class:`!SurfaceMesh`. ")

		.def("set_cutting_planes", [](SurfaceMesh& meshObj, py::array_t<FloatType> array) {
				ensureDataObjectIsMutable(meshObj);
				if(array.ndim() != 2) throw py::value_error("Array must be two-dimensional.");
				if(array.shape()[1] != 4) throw py::value_error("Second array dimension must have length 4.");
				QVector<Plane3> planes(array.shape()[0]);
				auto r = array.unchecked();
				for(size_t i = 0; i < planes.size(); i++) {
					planes[i].normal.x() = r(i,0);
					planes[i].normal.y() = r(i,1);
					planes[i].normal.z() = r(i,2);
					planes[i].dist = r(i,3);
				}
				meshObj.setCuttingPlanes(std::move(planes));
			}, py::arg("planes"),
			"set_cutting_planes(planes)"
			"\n\n"
			"Sets the cutting planes to be applied to this :py:class:`!SurfaceMesh`. "
			"The array *planes* must follow the same format as the one returned by :py:meth:`.get_cutting_planes`. ")
	;
	createDataPropertyAccessors(SurfaceMesh_py, "space_filling_region", &SurfaceMesh::spaceFillingRegion, &SurfaceMesh::setSpaceFillingRegion,
		"Indicates the index of the spatial region that fills the entire domain in case the surface is degenerate, i.e. the mesh has zero faces.");
	createDataSubobjectAccessors(SurfaceMesh_py, "domain", &PeriodicDomainDataObject::domain, &PeriodicDomainDataObject::setDomain,
		"The :py:class:`~ovito.data.SimulationCell` describing the (possibly periodic) domain which this "
		"surface mesh is embedded in. Note that this cell generally is indepenent of and may be different from the :py:attr:`~ovito.data.DataCollection.cell` "
		"found in the :py:class:`~ovito.data.DataCollection`. ");
	createDataSubobjectAccessors(SurfaceMesh_py, "vertices", &SurfaceMesh::vertices, &SurfaceMesh::setVertices,
		"The :py:class:`PropertyContainer` storing the vertex properties of the mesh, including the vertex coordinates. ");
	createDataSubobjectAccessors(SurfaceMesh_py, "faces", &SurfaceMesh::faces, &SurfaceMesh::setFaces,
		"The :py:class:`PropertyContainer` storing the properties of the faces of the mesh. ");
	createDataSubobjectAccessors(SurfaceMesh_py, "regions", &SurfaceMesh::regions, &SurfaceMesh::setRegions,
		"The :py:class:`PropertyContainer` storing the properties of the spatial regions of the mesh. ");

	ovito_class<SurfaceMeshVis, TransformingDataVis>(m,
			":Base class: :py:class:`ovito.vis.DataVis`"
			"\n\n"
			"Controls the visual appearance of a :py:class:`~ovito.data.SurfaceMesh` object, which is typically generated by modifiers such as "
			":py:class:`~ovito.modifiers.ConstructSurfaceModifier` or :py:class:`~ovito.modifiers.CreateIsosurfaceModifier`. "
			"See also the corresponding :ovitoman:`user manual page <../../display_objects.surface_mesh>` for this visual element. "
			,
			// Python class name:
			"SurfaceMeshVis")
		.def_property("surface_color", &SurfaceMeshVis::surfaceColor, &SurfaceMeshVis::setSurfaceColor,
				"The display color of the surface mesh."
				"\n\n"
				":Default: ``(1.0, 1.0, 1.0)``\n")
		.def_property("cap_color", &SurfaceMeshVis::capColor, &SurfaceMeshVis::setCapColor,
				"The display color of the cap polygons at periodic boundaries."
				"\n\n"
				":Default: ``(0.8, 0.8, 1.0)``\n")
		.def_property("show_cap", &SurfaceMeshVis::showCap, &SurfaceMeshVis::setShowCap,
				"Controls the visibility of cap polygons, which are created at the intersection of the surface mesh with periodic box boundaries."
				"\n\n"
				":Default: ``True``\n")
		.def_property("surface_transparency", &SurfaceMeshVis::surfaceTransparency, &SurfaceMeshVis::setSurfaceTransparency,
				"The level of transparency of the displayed surface. Valid range is 0.0 -- 1.0."
				"\n\n"
				":Default: 0.0\n")
		.def_property("cap_transparency", &SurfaceMeshVis::capTransparency, &SurfaceMeshVis::setCapTransparency,
				"The level of transparency of the displayed cap polygons. Valid range is 0.0 -- 1.0."
				"\n\n"
				":Default: 0.0\n")
		.def_property("smooth_shading", &SurfaceMeshVis::smoothShading, &SurfaceMeshVis::setSmoothShading,
				"Enables smooth shading of the triangulated surface mesh."
				"\n\n"
				":Default: ``True``\n")
		.def_property("highlight_edges", &SurfaceMeshVis::highlightEdges, &SurfaceMeshVis::setHighlightEdges,
				"Activates the highlighted rendering of the polygonal edges of the mesh."
				"\n\n"
				":Default: ``False``\n")
		.def_property("reverse_orientation", &SurfaceMeshVis::reverseOrientation, &SurfaceMeshVis::setReverseOrientation,
				"Flips the orientation of the surface. This affects the generation of cap polygons."
				"\n\n"
				":Default: ``False``\n")
	;

	ovito_class<RenderableSurfaceMesh, TransformedDataObject>{m}
	;

	ovito_class<VTKFileImporter, FileSourceImporter>{m}
	;

	ovito_class<VTKTriangleMeshExporter, FileExporter>{m}
	;

	ovito_class<WavefrontOBJImporter, FileSourceImporter>{m}
	;

	ovito_class<STLImporter, FileSourceImporter>{m}
	;
}

OVITO_REGISTER_PLUGIN_PYTHON_INTERFACE(MeshPython);

}	// End of namespace
}	// End of namespace