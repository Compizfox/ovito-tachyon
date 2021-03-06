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

#include <ovito/mesh/Mesh.h>
#include <ovito/core/utilities/io/CompressedTextReader.h>
#include <ovito/mesh/tri/TriMeshObject.h>
#include "STLImporter.h"
#include "TriMeshFrameData.h"

namespace Ovito { namespace Mesh {

IMPLEMENT_OVITO_CLASS(STLImporter);

/******************************************************************************
* Returns whether this importer class supports importing data of the given type.
******************************************************************************/
bool STLImporter::OOMetaClass::supportsDataType(const DataObject::OOMetaClass& dataObjectType) const
{
	return TriMeshObject::OOClass().isDerivedFrom(dataObjectType);
}

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool STLImporter::OOMetaClass::checkFileFormat(const FileHandle& file) const
{
	// Open input file.
	CompressedTextReader stream(file);

	// Read first line.
	stream.readLine(256);
	if(!stream.lineStartsWithToken("solid"))
		return false;

	// Read a couple of more lines until we find the first "facet normal" line, just to make sure.
	while(!stream.eof()) {
		const char* line = stream.readLineTrimLeft();
		if(stream.lineStartsWithToken("facet normal", true))
			return true;
		if(line[0] != '\0')
			return false;
	}

	return false;
}

/******************************************************************************
* Parses the given input file and stores the data in the given container object.
******************************************************************************/
FileSourceImporter::FrameDataPtr STLImporter::FrameLoader::loadFile()
{
	// Open file for reading.
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Reading STL file %1").arg(fileHandle().toString()));
	setProgressMaximum(stream.underlyingSize());

	// Jump to byte offset.
	if(frame().byteOffset != 0)
		stream.seek(frame().byteOffset, frame().lineNumber);

	// Create output data structure.
	auto frameData = std::make_shared<TriMeshFrameData>();
	TriMesh& mesh = frameData->mesh();

	// Read first line and make sure it's an STL file.
	stream.readLine();
	if(!stream.lineStartsWithToken("solid"))
		throw Exception(tr("Invalid STL file. Expected 'solid' keyword in first line but found '%2'").arg(stream.lineNumber()).arg(stream.lineString()));

	// Parse file line by line.
	int nVertices = -1;
	int vindices[3];
	while(!stream.eof()) {
		const char* line = stream.readLineTrimLeft();

		if(line[0] == '\0')
			continue;	// Skip empty lines.

		if(stream.lineStartsWithToken("facet normal", true) || stream.lineStartsWithToken("endfacet", true)) {
			// Ignore these lines.
		}
		else if(stream.lineStartsWithToken("outer loop", true)) {
			// Begin a new face.
			nVertices = 0;
		}
		else if(stream.lineStartsWithToken("vertex", true)) {
			if(nVertices == -1)
				throw Exception(tr("Unexpected vertex specification in line %1 of STL file").arg(stream.lineNumber()));
			// Parse face vertex.
			Point3 xyz;
			if(sscanf(line, "vertex " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING, &xyz.x(), &xyz.y(), &xyz.z()) != 3)
				throw Exception(tr("Invalid vertex specification in line %1 of STL file: %2").arg(stream.lineNumber()).arg(stream.lineString()));
			vindices[std::min(nVertices,2)] = mesh.addVertex(xyz);
			nVertices++;
			// Emit a new face to triangulate the polygon.
			if(nVertices >= 3) {
				TriMeshFace& f = mesh.addFace();
				f.setVertices(vindices[0], vindices[1], vindices[2]);
				if(nVertices == 3)
					f.setEdgeVisibility(true, true, false);
				else
					f.setEdgeVisibility(false, true, false);
				vindices[1] = vindices[2];
			}
		}
		else if(stream.lineStartsWithToken("endloop", true)) {
			// Close the current face.
			if(nVertices >= 3)
				mesh.faces().back().setEdgeVisible(2);
			nVertices = -1;
		}
		else if(stream.lineStartsWithToken("endsolid", true)) {
			break;	// End of file.
		}
		else {
			throw Exception(tr("Unknown keyword encountered in line %1 of STL file: %2").arg(stream.lineNumber()).arg(stream.lineString()));
		}

		if(!setProgressValueIntermittent(stream.underlyingByteOffset()))
			return {};
	}

	// STL files do not use shared vertices.
	// Try to unite identical vertices now.
	mesh.removeDuplicateVertices(1e-8 * mesh.boundingBox().size().length());

	mesh.determineEdgeVisibility();

	frameData->setStatus(tr("%1 vertices, %2 triangles").arg(mesh.vertexCount()).arg(mesh.faceCount()));
	return frameData;
}

}	// End of namespace
}	// End of namespace
