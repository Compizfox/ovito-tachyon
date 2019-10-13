////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2013 Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
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

// Inputs from calling program:
uniform mat4 modelview_matrix;
uniform mat4 modelview_projection_matrix;
uniform float modelview_uniform_scale;
uniform int pickingBaseID;
uniform int verticesPerElement;

#if __VERSION__ >= 130
	in vec3 position;
#else
	#define in attribute
	#define out varying
	#define flat
#endif

// The vertex data
in float vertexID;

// The cylinder data:
in vec3 cylinder_base;				// The position of the cylinder in model coordinates.
in vec3 cylinder_axis;				// The axis of the cylinder in model coordinates.
in float cylinder_radius;			// The radius of the cylinder in model coordinates.

// Outputs to fragment shader
flat out vec4 cylinder_color_fs;		// The base color of the cylinder.
flat out vec3 cylinder_view_base;		// Transformed cylinder position in view coordinates
flat out vec3 cylinder_view_axis;		// Transformed cylinder axis in view coordinates
flat out float cylinder_radius_sq_fs;	// The squared radius of the cylinder
flat out float cylinder_length;			// The length of the cylinder

void main()
{
#if __VERSION__ >= 130

	// Compute color from object ID.
	int objectID = pickingBaseID + (int(vertexID) / verticesPerElement);
	cylinder_color_fs = vec4(
		float(objectID & 0xFF) / 255.0,
		float((objectID >> 8) & 0xFF) / 255.0,
		float((objectID >> 16) & 0xFF) / 255.0,
		float((objectID >> 24) & 0xFF) / 255.0);

	// Transform and project vertex position.
	gl_Position = modelview_projection_matrix * vec4(position, 1.0);

#else

	// Compute color from object ID.
	float objectID = pickingBaseID + floor(vertexID / verticesPerElement);
	cylinder_color_fs = vec4(
		floor(mod(objectID, 256.0)) / 255.0,
		floor(mod(objectID / 256.0, 256.0)) / 255.0,
		floor(mod(objectID / 65536.0, 256.0)) / 255.0,
		floor(mod(objectID / 16777216.0, 256.0)) / 255.0);

	// Transform and project vertex position.
	gl_Position = modelview_projection_matrix * gl_Vertex;

#endif

	// Pass radius to fragment shader.
	cylinder_radius_sq_fs = cylinder_radius * modelview_uniform_scale;
	cylinder_radius_sq_fs *= cylinder_radius_sq_fs;

	// Transform cylinder to eye coordinates.
	cylinder_view_base = vec3(modelview_matrix * vec4(cylinder_base, 1));
	cylinder_view_axis = vec3(modelview_matrix * vec4(cylinder_axis, 0));

	// Pass length to fragment shader.
	cylinder_length = length(cylinder_view_axis);
}
