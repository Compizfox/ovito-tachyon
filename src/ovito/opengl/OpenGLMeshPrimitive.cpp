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

#include <ovito/core/Core.h>
#include "OpenGLMeshPrimitive.h"
#include "OpenGLSceneRenderer.h"

namespace Ovito {

/******************************************************************************
* Constructor.
******************************************************************************/
OpenGLMeshPrimitive::OpenGLMeshPrimitive(OpenGLSceneRenderer* renderer) :
	_contextGroup(QOpenGLContextGroup::currentContextGroup())
{
	OVITO_ASSERT(renderer->glcontext()->shareGroup() == _contextGroup);

	// Initialize OpenGL shader.
	_shader = renderer->loadShaderProgram("mesh", ":/openglrenderer/glsl/mesh/mesh.vs", ":/openglrenderer/glsl/mesh/mesh.fs");
	_pickingShader = renderer->loadShaderProgram("mesh.picking", ":/openglrenderer/glsl/mesh/picking/mesh.vs", ":/openglrenderer/glsl/mesh/picking/mesh.fs");
	_lineShader = renderer->loadShaderProgram("wireframe_line", ":/openglrenderer/glsl/lines/line.vs", ":/openglrenderer/glsl/lines/line.fs");
}

/******************************************************************************
* Sets the mesh to be stored in this buffer object.
******************************************************************************/
void OpenGLMeshPrimitive::setMesh(const TriMesh& mesh, const ColorA& meshColor, bool emphasizeEdges, DepthSortingMode depthSortingMode)
{
	OVITO_ASSERT(QOpenGLContextGroup::currentContextGroup() == _contextGroup);

	// Allocate render vertex buffer.
	_vertexBuffer.create(QOpenGLBuffer::StaticDraw, mesh.faceCount(), 3);
	if((mesh.hasVertexColors() || mesh.hasFaceColors()) && meshColor.a() == 1.0)
		_alpha = 1.0;
	else {
		if(materialColors().empty()) {
			_alpha = meshColor.a();
		}
		else {
			_alpha = 1.0;
			for(const ColorA& c : materialColors()) {
				if(c.a() != 1.0) {
					_alpha = c.a();
					break;
				}
			}
		}
	}

	// Discard any previous polygon edges.
	_edgeLinesBuffer.destroy();

	if(mesh.faceCount() == 0)
		return;

	ColoredVertexWithNormal* renderVertices = _vertexBuffer.map();

	if(!mesh.hasNormals()) {
		quint32 allMask = 0;

		// Compute face normals.
		std::vector<Vector_3<float>> faceNormals(mesh.faceCount());
		auto faceNormal = faceNormals.begin();
		for(auto face = mesh.faces().constBegin(); face != mesh.faces().constEnd(); ++face, ++faceNormal) {
			const Point3& p0 = mesh.vertex(face->vertex(0));
			Vector3 d1 = mesh.vertex(face->vertex(1)) - p0;
			Vector3 d2 = mesh.vertex(face->vertex(2)) - p0;
			*faceNormal = static_cast<Vector_3<float>>(d1.cross(d2));
			if(*faceNormal != Vector_3<float>::Zero()) {
				//faceNormal->normalize();
				allMask |= face->smoothingGroups();
			}
		}

		// Initialize render vertices.
		ColoredVertexWithNormal* rv = renderVertices;
		faceNormal = faceNormals.begin();
		ColorAT<float> defaultVertexColor = static_cast<ColorAT<float>>(meshColor);
		for(auto face = mesh.faces().constBegin(); face != mesh.faces().constEnd(); ++face, ++faceNormal) {

			// Initialize render vertices for this face.
			for(size_t v = 0; v < 3; v++, rv++) {
				if(face->smoothingGroups())
					rv->normal = Vector_3<float>::Zero();
				else
					rv->normal = *faceNormal;
				rv->pos = static_cast<Point_3<float>>(mesh.vertex(face->vertex(v)));
				if(mesh.hasVertexColors()) {
					rv->color = static_cast<ColorAT<float>>(mesh.vertexColor(face->vertex(v)));
					if(rv->color.a() != 1) _alpha = rv->color.a();
					else if(meshColor.a() != 1) rv->color.a() = meshColor.a();
				}
				else if(mesh.hasFaceColors()) {
					rv->color = static_cast<ColorAT<float>>(mesh.faceColor(face - mesh.faces().constBegin()));
					if(rv->color.a() != 1) _alpha = rv->color.a();
					else if(meshColor.a() != 1) rv->color.a() = meshColor.a();
				}
				else if(face->materialIndex() < materialColors().size() && face->materialIndex() >= 0) {
					rv->color = static_cast<ColorAT<float>>(materialColors()[face->materialIndex()]);
				}
				else {
					rv->color = defaultVertexColor;
				}
			}
		}

		if(allMask) {
			std::vector<Vector_3<float>> groupVertexNormals(mesh.vertexCount());
			for(int group = 0; group < OVITO_MAX_NUM_SMOOTHING_GROUPS; group++) {
				quint32 groupMask = quint32(1) << group;
				if((allMask & groupMask) == 0)
					continue;	// Group is not used.

				// Reset work arrays.
				std::fill(groupVertexNormals.begin(), groupVertexNormals.end(), Vector_3<float>::Zero());

				// Compute vertex normals at original vertices for current smoothing group.
				faceNormal = faceNormals.begin();
				for(auto face = mesh.faces().constBegin(); face != mesh.faces().constEnd(); ++face, ++faceNormal) {
					// Skip faces that do not belong to the current smoothing group.
					if((face->smoothingGroups() & groupMask) == 0) continue;

					// Add face's normal to vertex normals.
					for(size_t fv = 0; fv < 3; fv++)
						groupVertexNormals[face->vertex(fv)] += *faceNormal;
				}

				// Transfer vertex normals from original vertices to render vertices.
				rv = renderVertices;
				for(const auto& face : mesh.faces()) {
					if(face.smoothingGroups() & groupMask) {
						for(size_t fv = 0; fv < 3; fv++, ++rv)
							rv->normal += groupVertexNormals[face.vertex(fv)];
					}
					else rv += 3;
				}
			}
		}
	}
	else {
		// Use normals stored in the mesh.
		ColoredVertexWithNormal* rv = renderVertices;
		const Vector3* faceNormal = mesh.normals().begin();
		ColorAT<float> defaultVertexColor = static_cast<ColorAT<float>>(meshColor);
		for(auto face = mesh.faces().constBegin(); face != mesh.faces().constEnd(); ++face) {
			// Initialize render vertices for this face.
			for(size_t v = 0; v < 3; v++, rv++) {
				rv->normal = static_cast<Vector_3<float>>(*faceNormal++);
				rv->pos = static_cast<Point_3<float>>(mesh.vertex(face->vertex(v)));
				if(mesh.hasVertexColors()) {
					rv->color = static_cast<ColorAT<float>>(mesh.vertexColor(face->vertex(v)));
					if(rv->color.a() != 1) _alpha = rv->color.a();
					else if(meshColor.a() != 1) rv->color.a() = meshColor.a();
				}
				else if(mesh.hasFaceColors()) {
					rv->color = static_cast<ColorAT<float>>(mesh.faceColor(face - mesh.faces().constBegin()));
					if(rv->color.a() != 1) _alpha = rv->color.a();
					else if(meshColor.a() != 1) rv->color.a() = meshColor.a();
				}
				else if(face->materialIndex() >= 0 && face->materialIndex() < materialColors().size()) {
					rv->color = static_cast<ColorAT<float>>(materialColors()[face->materialIndex()]);
				}
				else {
					rv->color = defaultVertexColor;
				}
			}
		}
	}

	_vertexBuffer.unmap();

	// Save a list of coordinates which will be used to sort faces back-to-front.
	if(_alpha != 1 && depthSortingMode != MeshPrimitive::ConvexShapeMode) {
		_triangleDepthSortData.resize(mesh.faceCount());
		auto tc = _triangleDepthSortData.begin();
		for(auto face = mesh.faces().constBegin(); face != mesh.faces().constEnd(); ++face, ++tc) {
			// Compute centroid of triangle.
			const auto& v1 = mesh.vertex(face->vertex(0));
			const auto& v2 = mesh.vertex(face->vertex(1));
			const auto& v3 = mesh.vertex(face->vertex(2));
			tc->x() = (float)(v1.x() + v2.x() + v3.x()) / 3.0f;
			tc->y() = (float)(v1.y() + v2.y() + v3.y()) / 3.0f;
			tc->z() = (float)(v1.z() + v2.z() + v3.z()) / 3.0f;
		}
	}
	else _triangleDepthSortData.clear();

	// Create buffer for rendering polygon edges.
	if(emphasizeEdges) {
		// Count how many polygon edge are in the mesh.
		int numVisibleEdges = 0;
		for(const TriMeshFace& face : mesh.faces()) {
			for(size_t e = 0; e < 3; e++)
				if(face.edgeVisible(e)) numVisibleEdges++;
		}
		// Allocate storage buffer for line elements.
		_edgeLinesBuffer.create(QOpenGLBuffer::StaticDraw, numVisibleEdges, 2);
		Point_3<float>* lineVertices = _edgeLinesBuffer.map();

		// Generate line elements.
		for(const TriMeshFace& face : mesh.faces()) {
			for(size_t e = 0; e < 3; e++) {
				if(face.edgeVisible(e)) {
					*lineVertices++ = Point_3<float>(mesh.vertex(face.vertex(e)));
					*lineVertices++ = Point_3<float>(mesh.vertex(face.vertex((e+1)%3)));
				}
			}
		}

		_edgeLinesBuffer.unmap();
	}
}

/******************************************************************************
* Activates rendering of multiple instances of the mesh.
******************************************************************************/
void OpenGLMeshPrimitive::setInstancedRendering(std::vector<AffineTransformation> perInstanceTMs, std::vector<ColorA> perInstanceColors)
{
	OVITO_ASSERT(perInstanceTMs.size() == perInstanceColors.size() || perInstanceColors.empty());
	_alpha = std::any_of(perInstanceColors.begin(), perInstanceColors.end(), [](const ColorA& c) { return c.a() != FloatType(1); }) ? 0.5 : 1.0;
	_perInstanceTMs = std::move(perInstanceTMs);
	_perInstanceColors = std::move(perInstanceColors);
	_useInstancedRendering = true;
}

/******************************************************************************
* Returns true if the geometry buffer is filled and can be rendered with the given renderer.
******************************************************************************/
bool OpenGLMeshPrimitive::isValid(SceneRenderer* renderer)
{
	OpenGLSceneRenderer* vpRenderer = qobject_cast<OpenGLSceneRenderer*>(renderer);
	if(!vpRenderer) return false;
	return _vertexBuffer.isCreated() && (_contextGroup == vpRenderer->glcontext()->shareGroup());
}

/******************************************************************************
* Renders the geometry.
******************************************************************************/
void OpenGLMeshPrimitive::render(SceneRenderer* renderer)
{
	OVITO_ASSERT(_contextGroup == QOpenGLContextGroup::currentContextGroup());
	OpenGLSceneRenderer* vpRenderer = dynamic_object_cast<OpenGLSceneRenderer>(renderer);

	if(faceCount() <= 0 || !vpRenderer || (_useInstancedRendering && _perInstanceTMs.empty()))
		return;

	// If object is translucent, don't render it during the first rendering pass.
	// Queue primitive so that it gets rendered during the second pass.
	if(!renderer->isPicking() && _alpha != 1 && vpRenderer->translucentPass() == false) {
		vpRenderer->registerTranslucentPrimitive(shared_from_this());
		return;
	}

	vpRenderer->rebindVAO();

	// Render wireframe edges.
	if(!renderer->isPicking() && _edgeLinesBuffer.isCreated()) {
		if(!_lineShader->bind())
			vpRenderer->throwException(QStringLiteral("Failed to bind OpenGL shader."));
		ColorA wireframeColor(0.1, 0.1, 0.1, _alpha);
		if(vpRenderer->glformat().majorVersion() >= 3) {
			OVITO_CHECK_OPENGL(vpRenderer, _lineShader->setAttributeValue("color", wireframeColor.r(), wireframeColor.g(), wireframeColor.b(), wireframeColor.a()));
		}
#ifndef Q_OS_WASM	
		else if(vpRenderer->oldGLFunctions()) {
			// Older OpenGL implementations cannot take vertex colors through a custom shader attribute.
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->oldGLFunctions()->glColor4f(wireframeColor.r(), wireframeColor.g(), wireframeColor.b(), wireframeColor.a()));
		}
#endif		
		if(_alpha != 1.0) {
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glEnable(GL_BLEND));
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glBlendEquation(GL_FUNC_ADD));
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE_MINUS_DST_COLOR, GL_ONE));
		}
		_edgeLinesBuffer.bindPositions(vpRenderer, _lineShader);
		Matrix4 mvp_matrix = vpRenderer->projParams().projectionMatrix * vpRenderer->modelViewTM();
		if(!_useInstancedRendering) {
			_lineShader->setUniformValue("modelview_projection_matrix", (QMatrix4x4)mvp_matrix);
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDrawArrays(GL_LINES, 0, _edgeLinesBuffer.elementCount() * _edgeLinesBuffer.verticesPerElement()));
		}
		else {
			if(_alpha == 1.0) {
				for(const AffineTransformation& instanceTM : _perInstanceTMs) {
					_lineShader->setUniformValue("modelview_projection_matrix", (QMatrix4x4)(mvp_matrix * instanceTM));
					OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDrawArrays(GL_LINES, 0, _edgeLinesBuffer.elementCount() * _edgeLinesBuffer.verticesPerElement()));
				}
			}
			else {
				OVITO_ASSERT(_perInstanceColors.size() == _perInstanceTMs.size());
				auto instanceColor = _perInstanceColors.cbegin();
				for(const AffineTransformation& instanceTM : _perInstanceTMs) {
					_lineShader->setUniformValue("modelview_projection_matrix", (QMatrix4x4)(mvp_matrix * instanceTM));
					wireframeColor.a() = instanceColor->a();
					++instanceColor;
					if(vpRenderer->glformat().majorVersion() >= 3) {
						OVITO_CHECK_OPENGL(vpRenderer, _lineShader->setAttributeValue("color", wireframeColor.r(), wireframeColor.g(), wireframeColor.b(), wireframeColor.a()));
					}
#ifndef Q_OS_WASM	
					else if(vpRenderer->oldGLFunctions()) {
						// Older OpenGL implementations cannot take vertex colors through a custom shader attribute.
						OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->oldGLFunctions()->glColor4f(wireframeColor.r(), wireframeColor.g(), wireframeColor.b(), wireframeColor.a()));
					}
#endif					
					OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDrawArrays(GL_LINES, 0, _edgeLinesBuffer.elementCount() * _edgeLinesBuffer.verticesPerElement()));
				}
			}
		}
		_edgeLinesBuffer.detachPositions(vpRenderer, _lineShader);
		_lineShader->release();
		OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glEnable(GL_POLYGON_OFFSET_FILL));
		OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glPolygonOffset(1.0f, 1.0f));
		if(_alpha != 1.0)
			vpRenderer->glDisable(GL_BLEND);
	}

	OVITO_REPORT_OPENGL_ERRORS(vpRenderer);

	if(cullFaces()) {
		OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glEnable(GL_CULL_FACE));
		OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glCullFace(GL_BACK));
	}
	else {
		OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDisable(GL_CULL_FACE));
	}

	QOpenGLShaderProgram* shader;
	if(!renderer->isPicking())
		shader = _shader;
	else
		shader = _pickingShader;

	if(!shader->bind())
		renderer->throwException(QStringLiteral("Failed to bind OpenGL shader."));

	_vertexBuffer.bindPositions(vpRenderer, shader, offsetof(ColoredVertexWithNormal, pos));
	if(!renderer->isPicking()) {
		if(_alpha != 1.0) {
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glEnable(GL_BLEND));
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glBlendEquation(GL_FUNC_ADD));
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE_MINUS_DST_COLOR, GL_ONE));
		}
		_vertexBuffer.bindNormals(vpRenderer, shader, offsetof(ColoredVertexWithNormal, normal));
	}
	else {
		vpRenderer->activateVertexIDs(_pickingShader, _vertexBuffer.elementCount() * _vertexBuffer.verticesPerElement());
	}
	OVITO_REPORT_OPENGL_ERRORS(vpRenderer);

	size_t numInstances = !_useInstancedRendering ? 1 : _perInstanceTMs.size();
	for(size_t instance = 0; instance < numInstances; instance++) {

		AffineTransformation mv_matrix;
		if(_useInstancedRendering)
			mv_matrix = vpRenderer->modelViewTM() * _perInstanceTMs[instance];
		else
			mv_matrix = vpRenderer->modelViewTM();
		OVITO_CHECK_OPENGL(vpRenderer, shader->setUniformValue("modelview_projection_matrix", (QMatrix4x4)(vpRenderer->projParams().projectionMatrix * mv_matrix)));
		if(!renderer->isPicking()) {
			OVITO_CHECK_OPENGL(vpRenderer, shader->setUniformValue("normal_matrix", (QMatrix3x3)(mv_matrix.linear().inverse().transposed())));
			if(!_useInstancedRendering || _perInstanceColors.empty()) {
				_vertexBuffer.bindColors(vpRenderer, shader, 4, offsetof(ColoredVertexWithNormal, color));
			}
			else {
				const ColorA& color = _perInstanceColors[instance];
				if(vpRenderer->glformat().majorVersion() >= 3) {
					OVITO_CHECK_OPENGL(vpRenderer, shader->setAttributeValue("color", color.r(), color.g(), color.b(), color.a()));
				}
#ifndef Q_OS_WASM	
				else if(vpRenderer->oldGLFunctions()) {
					// Older OpenGL implementations cannot take colors through a custom shader attribute.
					OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->oldGLFunctions()->glColor4f(color.r(), color.g(), color.b(), color.a()));
				}
#endif				
			}
		}
		else {
			if(!_useInstancedRendering) {
				OVITO_CHECK_OPENGL(vpRenderer, _pickingShader->setUniformValue("pickingBaseID", (GLint)vpRenderer->registerSubObjectIDs(faceCount())));
				OVITO_CHECK_OPENGL(vpRenderer, _pickingShader->setUniformValue("vertexIdDivisor", (GLint)3));
			}
			else {
				OVITO_CHECK_OPENGL(vpRenderer, _pickingShader->setUniformValue("pickingBaseID", (GLint)vpRenderer->registerSubObjectIDs(1)));
				OVITO_CHECK_OPENGL(vpRenderer, _pickingShader->setUniformValue("vertexIdDivisor", (GLint)faceCount()*3));
			}
		}

		if(!renderer->isPicking() && _alpha != 1.0) {
			if(!_triangleDepthSortData.empty()) {
				OVITO_ASSERT(_triangleDepthSortData.size() == faceCount());
				OVITO_ASSERT(_vertexBuffer.verticesPerElement() == 3);
				
				// Render faces in back-to-front order to avoid artifacts at overlapping translucent faces.
				std::vector<GLuint> indices(faceCount());
				std::iota(indices.begin(), indices.end(), 0);

				// First compute distance of each face from the camera along viewing direction (=camera z-axis).
				std::vector<float> distances(faceCount());
				Vector_3<float> direction = (Vector_3<float>)mv_matrix.inverse().column(2);
				std::transform(_triangleDepthSortData.begin(), _triangleDepthSortData.end(), distances.begin(), [direction](const Vector_3<float>& v) {
					return direction.dot(v);
				});
				
				// Now sort face indices with respect to distance (back-to-front order).
				std::sort(indices.begin(), indices.end(), [&](GLuint a, GLuint b) {
					return distances[a] < distances[b];
				});
				
				// Create OpenGL index buffer which can be used with glDrawElements.
				OpenGLBuffer<GLuint> primitiveIndices(QOpenGLBuffer::IndexBuffer);
				primitiveIndices.create(QOpenGLBuffer::StaticDraw, 3 * faceCount());
				GLuint* p = primitiveIndices.map();
				for(size_t i = 0; i < indices.size(); i++, p += 3)
					std::iota(p, p + 3, indices[i]*3);
				primitiveIndices.unmap();
				
				// Render triangles in depth-sorted order.
				primitiveIndices.oglBuffer().bind();
				OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDrawElements(GL_TRIANGLES, _vertexBuffer.elementCount() * _vertexBuffer.verticesPerElement(), GL_UNSIGNED_INT, nullptr));
				primitiveIndices.oglBuffer().release();
			}
			else {
				// Assuming that the input mesh is convex, render semi-transparent triangles in two passes: 
				// First, render triangles facing away from the viewer, then render triangles facing toward the viewer.
				// Each time we pass the entire triangle list to OpenGL and use OpenGL's backface/frontfrace culling
				// option to render the right subset of triangles.
				OVITO_REPORT_OPENGL_ERRORS(vpRenderer);
				OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glEnable(GL_CULL_FACE));
				if(!cullFaces()) {
					// First pass is only needed if backface culling is not active.
					OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glCullFace(GL_FRONT));
					OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDrawArrays(GL_TRIANGLES, 0, _vertexBuffer.elementCount() * _vertexBuffer.verticesPerElement()));
				}
				OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glCullFace(GL_BACK));
				OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDrawArrays(GL_TRIANGLES, 0, _vertexBuffer.elementCount() * _vertexBuffer.verticesPerElement()));
				OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDisable(GL_CULL_FACE));
			}
		}
		else {
			// Render faces in arbitrary order.
			OVITO_REPORT_OPENGL_ERRORS(vpRenderer);
			OVITO_CHECK_OPENGL(vpRenderer, vpRenderer->glDrawArrays(GL_TRIANGLES, 0, _vertexBuffer.elementCount() * _vertexBuffer.verticesPerElement()));
		}
	}

	if(!renderer->isPicking() && _edgeLinesBuffer.isCreated()) {
		vpRenderer->glDisable(GL_POLYGON_OFFSET_FILL);
	}

	_vertexBuffer.detachPositions(vpRenderer, shader);
	if(!renderer->isPicking()) {
		if(!_useInstancedRendering)
			_vertexBuffer.detachColors(vpRenderer, shader);
		_vertexBuffer.detachNormals(vpRenderer, shader);
		if(_alpha != 1.0) 
			vpRenderer->glDisable(GL_BLEND);
	}
	else {
		vpRenderer->deactivateVertexIDs(_pickingShader);
	}
	shader->release();

	OVITO_REPORT_OPENGL_ERRORS(vpRenderer);

	// Restore old state.
	if(cullFaces()) {
		vpRenderer->glDisable(GL_CULL_FACE);
		vpRenderer->glCullFace(GL_BACK);
	}
}

}	// End of namespace
