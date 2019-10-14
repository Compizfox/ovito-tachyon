////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2013 Alexander Stukowski
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

uniform mat4 modelview_projection_matrix;
uniform int pickingBaseID;
uniform int vertexIdDivisor;

#if __VERSION__ >= 130

	// The input data:
	in vec3 position;

	// Output passed to fragment shader.
	flat out vec4 vertex_color_fs;

#else

	// The input data:
	attribute float vertexID;

#endif

void main()
{
	// Compute color from object ID.
#if __VERSION__ >= 130
	int objectID = pickingBaseID + gl_VertexID / vertexIdDivisor;
	vertex_color_fs = vec4(
		float(objectID & 0xFF) / 255.0,
		float((objectID >> 8) & 0xFF) / 255.0,
		float((objectID >> 16) & 0xFF) / 255.0,
		float((objectID >> 24) & 0xFF) / 255.0);

	gl_Position = modelview_projection_matrix * vec4(position, 1.0);
#else
	float objectID = pickingBaseID + floor(vertexID / vertexIdDivisor);
	gl_FrontColor = vec4(
		floor(mod(objectID, 256.0)) / 255.0,
		floor(mod(objectID / 256.0, 256.0)) / 255.0,
		floor(mod(objectID / 65536.0, 256.0)) / 255.0,
		floor(mod(objectID / 16777216.0, 256.0)) / 255.0);

	gl_Position = modelview_projection_matrix * gl_Vertex;
#endif
}
