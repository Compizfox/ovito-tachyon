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

#include <core/Core.h>
#include <core/scene/SceneNode.h>
#include <core/scene/SceneRoot.h>
#include <core/scene/ObjectNode.h>
#include <core/scene/GroupNode.h>
#include <core/scene/pipeline/PipelineObject.h>
#include <core/scene/pipeline/Modifier.h>
#include <core/scene/display/DisplayObject.h>
#include "ViewportSceneRenderer.h"
#include "ViewportSceneRendererEditor.h"

namespace Ovito {

IMPLEMENT_SERIALIZABLE_OVITO_OBJECT(Core, ViewportSceneRenderer, SceneRenderer);
SET_OVITO_OBJECT_EDITOR(ViewportSceneRenderer, ViewportSceneRendererEditor)
DEFINE_PROPERTY_FIELD(ViewportSceneRenderer, _antialiasingLevel, "AntialiasingLevel")
SET_PROPERTY_FIELD_LABEL(ViewportSceneRenderer, _antialiasingLevel, "Antialiasing level")

/******************************************************************************
* This method is called just before renderFrame() is called.
******************************************************************************/
void ViewportSceneRenderer::beginRender()
{
	OVITO_CHECK_OBJECT_POINTER(viewport());
	_glcontext = viewport()->_glcontext;
	OVITO_ASSERT(_glcontext == QOpenGLContext::currentContext());

	// Obtain a functions object that allows to call OpenGL 2.0 functions in a platform-independent way.
	_glFunctions = _glcontext->functions();

	// Obtain a functions object that allows to call OpenGL 3.0 functions in a platform-independent way.
	_glFunctions30 = _glcontext->versionFunctions<QOpenGLFunctions_3_0>();
	if(!_glFunctions30 || !_glFunctions30->initializeOpenGLFunctions())
		_glFunctions30 = nullptr;

	// Obtain a functions object that allows to call OpenGL 3.2 core functions in a platform-independent way.
	_glFunctions32 = _glcontext->versionFunctions<QOpenGLFunctions_3_2_Core>();
	if(!_glFunctions32 || !_glFunctions32->initializeOpenGLFunctions())
		_glFunctions32 = nullptr;

	OVITO_ASSERT(_glFunctions30 || _glFunctions32);

	// Obtain surface format.
	_glformat = _glcontext->format();

	// Set up a vertex array object. This is only required when using OpenGL 3.2 Core Profile.
	if(glformat().majorVersion() >= 3 && glformat().minorVersion() >= 2) {
		if(!_vertexArrayObject.isCreated())
			_vertexArrayObject.create();
		_vertexArrayObject.bind();
	}
}

/******************************************************************************
* This method is called after renderFrame() has been called.
******************************************************************************/
void ViewportSceneRenderer::endRender()
{
	if(_vertexArrayObject.isCreated())
		_vertexArrayObject.release();
	_glcontext = nullptr;
}

/******************************************************************************
* Renders the current animation frame.
******************************************************************************/
void ViewportSceneRenderer::renderFrame()
{
	// Clear background.
	Color backgroundColor = Viewport::viewportColor(ViewportSettings::COLOR_VIEWPORT_BKG);
	glClearColor(backgroundColor.r(), backgroundColor.g(), backgroundColor.b(), 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	renderScene();
}

/******************************************************************************
* Changes the current local to world transformation matrix.
******************************************************************************/
void ViewportSceneRenderer::setWorldTransform(const AffineTransformation& tm)
{
	_modelViewTM = projParams().viewMatrix * tm;
}

/******************************************************************************
* Translates an OpenGL error code to a human-readable message string.
******************************************************************************/
const char* ViewportSceneRenderer::openglErrorString(GLenum errorCode)
{
	switch(errorCode) {
	case GL_NO_ERROR: return "GL_NO_ERROR - No error has been recorded.";
	case GL_INVALID_ENUM: return "GL_INVALID_ENUM - An unacceptable value is specified for an enumerated argument.";
	case GL_INVALID_VALUE: return "GL_INVALID_VALUE - A numeric argument is out of range.";
	case GL_INVALID_OPERATION: return "GL_INVALID_OPERATION - The specified operation is not allowed in the current state.";
	case GL_STACK_OVERFLOW: return "GL_STACK_OVERFLOW - This command would cause a stack overflow.";
	case GL_STACK_UNDERFLOW: return "GL_STACK_UNDERFLOW - This command would cause a stack underflow.";
	case GL_OUT_OF_MEMORY: return "GL_OUT_OF_MEMORY - There is not enough memory left to execute the command.";
	case GL_TABLE_TOO_LARGE: return "GL_TABLE_TOO_LARGE - The specified table exceeds the implementation's maximum supported table size.";
	default: return "Unknown OpenGL error code.";
	}
}

/******************************************************************************
* Renders a single node.
******************************************************************************/
void ViewportSceneRenderer::renderNode(SceneNode* node)
{
    OVITO_CHECK_OBJECT_POINTER(node);

    // Setup transformation matrix.
	TimeInterval interval;
	const AffineTransformation& nodeTM = node->getWorldTransform(time(), interval);
	setWorldTransform(nodeTM);

	if(node->isObjectNode()) {
		ObjectNode* objNode = static_object_cast<ObjectNode>(node);

		// Do not render node if it is the view node of the viewport or
		// if it is the target of the view node.
		if(viewport() && viewport()->viewNode()) {
			if(viewport()->viewNode() == objNode || viewport()->viewNode()->targetNode() == objNode)
				return;
		}

		// Evaluate geometry pipeline of object node and render the results.
		objNode->render(time(), this);

		// Render the modifiers applied to the object.
		PipelineObject* pipelineObj = dynamic_object_cast<PipelineObject>(objNode->sceneObject());
		renderPipelineObject(pipelineObj, objNode);
	}

	// Continue with rendering the child nodes.
	SceneRenderer::renderNode(node);
}

/******************************************************************************
* Renders the selected modifiers.
******************************************************************************/
void ViewportSceneRenderer::renderPipelineObject(PipelineObject* pipelineObj, ObjectNode* objNode)
{
	if(!pipelineObj) return;

#if 0
	// Retrieve the modifier that is currently selected.
	Modifier* selectedModifier = dynamic_object_cast<Modifier>(MAIN_FRAME->commandPanel()->editObject());
	if(!selectedModifier) return;

	// Iterate over all modifier applications until the selected modifier is found.
	// Then render it.
	Q_FOREACH(ModifierApplication* modApp, modObj->modifierApplications()) {
		if(modApp->modifier() == selectedModifier) {

			// Setup transformation.
			TimeInterval iv;
			const AffineTransformation& nodeTM = objNode->getWorldTransform(currentTime, iv);
			viewport->setWorldMatrix(objNode->objectTransform() * nodeTM);

			// Render selected modifier.
			selectedModifier->renderModifier(currentTime, objNode, modApp, viewport);

			return;
		}
	}

	// Not found yet -> go on with nested modified objects.
	for(int i = 0; i < modObj->inputObjectCount(); i++)
		renderModifiedObject(dynamic_object_cast<ModifiedObject>(modObj->inputObject(i)), objNode);
#endif
}

/******************************************************************************
* Loads an OpenGL shader program.
******************************************************************************/
QOpenGLShaderProgram* ViewportSceneRenderer::loadShaderProgram(const QString& id, const QString& vertexShaderFile, const QString& fragmentShaderFile, const QString& geometryShaderFile)
{
	QOpenGLContextGroup* contextGroup = glcontext()->shareGroup();
	OVITO_ASSERT(contextGroup == QOpenGLContextGroup::currentContextGroup());

	OVITO_ASSERT(QOpenGLShaderProgram::hasOpenGLShaderPrograms());
	OVITO_ASSERT(QOpenGLShader::hasOpenGLShaders(QOpenGLShader::Vertex));
	OVITO_ASSERT(QOpenGLShader::hasOpenGLShaders(QOpenGLShader::Fragment));

	// The OpenGL shaders are only created once per OpenGL context group.
	QScopedPointer<QOpenGLShaderProgram> program(contextGroup->findChild<QOpenGLShaderProgram*>(id));
	if(program)
		return program.take();

	program.reset(new QOpenGLShaderProgram(contextGroup));
	program->setObjectName(id);

	// Load and compile vertex shader source.
	loadShader(program.data(), QOpenGLShader::Vertex, vertexShaderFile);

	// Load and compile fragment shader source.
	loadShader(program.data(), QOpenGLShader::Fragment, fragmentShaderFile);

	// Load and compile geometry shader source.
	if(!geometryShaderFile.isEmpty()) {
		loadShader(program.data(), QOpenGLShader::Geometry, geometryShaderFile);
	}

	if(!program->link()) {
		qDebug() << "OpenGL shader log:";
		qDebug() << program->log();
		throw Exception(QString("The OpenGL shader program %1 failed to link. See log for details.").arg(id));
	}

	OVITO_ASSERT(contextGroup->findChild<QOpenGLShaderProgram*>(id) == program.data());
	return program.take();
}

/******************************************************************************
* Loads and compiles a GLSL shader and adds it to the given program object.
******************************************************************************/
void ViewportSceneRenderer::loadShader(QOpenGLShaderProgram* program, QOpenGLShader::ShaderType shaderType, const QString& filename)
{
	// Load shader source.
	QFile shaderSourceFile(filename);
	if(!shaderSourceFile.open(QFile::ReadOnly))
		throw Exception(QString("Unable to open shader source file %1.").arg(filename));
    QByteArray shaderSource = shaderSourceFile.readAll();

    // Insert GLSL version string at the top.
	// Pick GLSL language version based on current OpenGL version.
	if(glformat().profile() == QSurfaceFormat::CoreProfile)
		shaderSource.prepend("#version 150\n");
	else
		shaderSource.prepend("#version 130\n");

	// Load and compile vertex shader source.
    if(!program->addShaderFromSourceCode(shaderType, shaderSource)) {
		qDebug() << "OpenGL shader log:";
		qDebug() << program->log();
		throw Exception(QString("The shader source file %1 failed to compile. See log for details.").arg(filename));
	}
}


};