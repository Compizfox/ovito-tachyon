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

precision highp float;

// Inputs from calling program:
uniform mat4 modelview_matrix;
uniform mat4 projection_matrix;
uniform vec2 imposter_texcoords[6];
uniform vec4 imposter_voffsets[6];
uniform float radius_scalingfactor;

// The particle data:
attribute vec3 position;
attribute vec4 color;
attribute float particle_radius;
attribute float vertexID;

// Output to fragment shader:
varying vec4 particle_color_fs;
varying float particle_radius_fs;	// The particle's radius.
varying float ze0;					// The particle's Z coordinate in eye coordinates.
varying vec2 texcoords;

void main()
{
	// Forward color to fragment shader.
	particle_color_fs = color;

	// Transform and project particle position.
	vec4 eye_position = modelview_matrix * vec4(position, 1.0);

	// Assign texture coordinates.
	int cornerIndex = int(mod(vertexID, 6.0));
	texcoords = imposter_texcoords[cornerIndex];

	// Transform and project particle position.
	gl_Position = projection_matrix * (eye_position + (particle_radius * radius_scalingfactor) * imposter_voffsets[cornerIndex]);

	// Forward particle radius to fragment shader.
	particle_radius_fs = particle_radius * radius_scalingfactor;

	// Pass particle position in eye coordinates to fragment shader.
	ze0 = eye_position.z;
}
