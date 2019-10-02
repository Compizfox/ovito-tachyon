///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2013) Alexander Stukowski
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

#pragma once


#include <ovito/core/Core.h>
#include <ovito/core/rendering/MeshPrimitive.h>
#include <ovito/core/utilities/mesh/TriMesh.h>
#include "OpenGLBuffer.h"

namespace Ovito { OVITO_BEGIN_INLINE_NAMESPACE(Rendering) OVITO_BEGIN_INLINE_NAMESPACE(Internal)

/**
 * \brief Buffer object that stores a triangle mesh to be rendered in the viewports.
 */
class OpenGLMeshPrimitive : public MeshPrimitive, public std::enable_shared_from_this<OpenGLMeshPrimitive>
{
public:

	/// Constructor.
	OpenGLMeshPrimitive(OpenGLSceneRenderer* renderer);

	/// Sets the mesh to be stored in this buffer object.
	virtual void setMesh(const TriMesh& mesh, const ColorA& meshColor, bool emphasizeEdges) override;

	/// \brief Returns the number of triangle faces stored in the buffer.
	virtual int faceCount() override { return _vertexBuffer.elementCount(); }

	/// \brief Returns true if the geometry buffer is filled and can be rendered with the given renderer.
	virtual bool isValid(SceneRenderer* renderer) override;

	/// \brief Renders the geometry.
	virtual void render(SceneRenderer* renderer) override;

	/// Activates rendering of multiple instances of the mesh.
	virtual void setInstancedRendering(std::vector<AffineTransformation> perInstanceTMs, std::vector<ColorA> perInstanceColors) override;

private:

	/// Stores data of a single vertex passed to the OpenGL implementation.
	struct ColoredVertexWithNormal {
		Point_3<float> pos;
		Vector_3<float> normal;
		ColorAT<float> color;
	};

	/// The internal OpenGL vertex buffer that stores the vertex data.
	OpenGLBuffer<ColoredVertexWithNormal> _vertexBuffer;

	/// The GL context group under which the GL vertex buffer has been created.
	QOpenGLContextGroup* _contextGroup;

	/// The OpenGL shader program used to render the triangles.
	QOpenGLShaderProgram* _shader;

	/// The OpenGL shader program used to render the triangles in picking mode.
	QOpenGLShaderProgram* _pickingShader;

	/// The OpenGL shader program used to render the wireframe edges.
	QOpenGLShaderProgram* _lineShader;

	/// Are we rendering a semi-transparent mesh?
	FloatType _alpha = 1.0;

	/// This is required to render translucent triangles in the correct order from back to front.
	std::vector<Point3> _triangleCoordinates;

	/// The internal OpenGL vertex buffer that stores the vertex data for rendering polygon edges.
	OpenGLBuffer<Point_3<float>> _edgeLinesBuffer;

	/// The list of transformation matrices when rendering multiple instances of the mesh.
	std::vector<AffineTransformation> _perInstanceTMs;

	/// The list of colors when rendering multiple instances of the mesh.
	std::vector<ColorA> _perInstanceColors;

	/// Activates the rendering of multiple instances of the same mesh.
	bool _useInstancedRendering = false;
};

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace
