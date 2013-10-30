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
#include <core/viewport/Viewport.h>
#include <core/viewport/ViewportManager.h>
#include <core/scene/ObjectNode.h>
#include <core/gui/properties/FloatParameterUI.h>
#include <core/gui/properties/BooleanParameterUI.h>
#include <core/rendering/RenderSettings.h>
#include <core/scene/display/camera/CameraDisplayObject.h>
#include "CameraObject.h"
#include "moc_AbstractCameraObject.cpp"

namespace Ovito {

IMPLEMENT_SERIALIZABLE_OVITO_OBJECT(Core, AbstractCameraObject, SceneObject)
IMPLEMENT_SERIALIZABLE_OVITO_OBJECT(Core, CameraObject, AbstractCameraObject)
IMPLEMENT_OVITO_OBJECT(Core, CameraObjectEditor, PropertiesEditor)
SET_OVITO_OBJECT_EDITOR(CameraObject, CameraObjectEditor)
DEFINE_PROPERTY_FIELD(CameraObject, _isPerspective, "IsPerspective")
DEFINE_REFERENCE_FIELD(CameraObject, _fov, "FOV", FloatController)
DEFINE_REFERENCE_FIELD(CameraObject, _zoom, "Zoom", FloatController)
SET_PROPERTY_FIELD_LABEL(CameraObject, _isPerspective, "Perspective projection")
SET_PROPERTY_FIELD_LABEL(CameraObject, _fov, "Field of View")
SET_PROPERTY_FIELD_LABEL(CameraObject, _zoom, "Zoom")
SET_PROPERTY_FIELD_UNITS(CameraObject, _fov, AngleParameterUnit)
SET_PROPERTY_FIELD_UNITS(CameraObject, _zoom, WorldParameterUnit)

/******************************************************************************
* Constructs a camera object.
******************************************************************************/
CameraObject::CameraObject() : _isPerspective(true)
{
	INIT_PROPERTY_FIELD(CameraObject::_isPerspective);
	INIT_PROPERTY_FIELD(CameraObject::_fov);
	INIT_PROPERTY_FIELD(CameraObject::_zoom);

	_fov = ControllerManager::instance().createDefaultController<FloatController>();
	_fov->setValue(0, FLOATTYPE_PI/3.0);
	_zoom = ControllerManager::instance().createDefaultController<FloatController>();
	_zoom->setValue(0, 100);

	setDisplayObject(new CameraDisplayObject());
}

/******************************************************************************
* Asks the object for its validity interval at the given time.
******************************************************************************/
TimeInterval CameraObject::objectValidity(TimePoint time)
{
	TimeInterval interval = TimeInterval::forever();
	if(isPerspective() && _fov) interval.intersect(_fov->validityInterval(time));
	if(!isPerspective() && _zoom) interval.intersect(_zoom->validityInterval(time));
	return interval;
}

/******************************************************************************
* Fills in the missing fields of the camera view descriptor structure.
******************************************************************************/
void CameraObject::projectionParameters(TimePoint time, ViewProjectionParameters& params)
{
	// Transform scene bounding box to camera space.
	Box3 bb = params.boundingBox.transformed(params.viewMatrix).centerScale(1.01);

	// Compute projection matrix.
	params.isPerspective = isPerspective();
	if(params.isPerspective) {
		if(bb.minc.z() < -FLOATTYPE_EPSILON) {
			params.zfar = -bb.minc.z();
			params.znear = std::max(-bb.maxc.z(), params.zfar * 1e-4f);
		}
		else {
			params.zfar = std::max(params.boundingBox.size().length(), FloatType(1));
			params.znear = params.zfar * 1e-4f;
		}
		params.zfar = std::max(params.zfar, params.znear * 1.01f);

		// Get the camera angle.
		_fov->getValue(time, params.fieldOfView, params.validityInterval);
		if(params.fieldOfView < FLOATTYPE_EPSILON) params.fieldOfView = FLOATTYPE_EPSILON;
		if(params.fieldOfView > FLOATTYPE_PI - FLOATTYPE_EPSILON) params.fieldOfView = FLOATTYPE_PI - FLOATTYPE_EPSILON;

		params.projectionMatrix = Matrix4::perspective(params.fieldOfView, 1.0 / params.aspectRatio, params.znear, params.zfar);
	}
	else {
		if(!bb.isEmpty()) {
			params.znear = -bb.maxc.z();
			params.zfar  = std::max(-bb.minc.z(), params.znear + 1.0f);
		}
		else {
			params.znear = 1;
			params.zfar = 100;
		}

		// Get the camera zoom.
		_zoom->getValue(time, params.fieldOfView, params.validityInterval);
		if(params.fieldOfView < FLOATTYPE_EPSILON) params.fieldOfView = FLOATTYPE_EPSILON;

		params.projectionMatrix = Matrix4::ortho(-params.fieldOfView / params.aspectRatio, params.fieldOfView / params.aspectRatio,
							-params.fieldOfView, params.fieldOfView, params.znear, params.zfar);
	}
	params.inverseProjectionMatrix = params.projectionMatrix.inverse();
}

/******************************************************************************
* Returns the field of view of the camera.
******************************************************************************/
FloatType CameraObject::fieldOfView(TimePoint time, TimeInterval& validityInterval)
{
	FloatType fov;
	if(isPerspective())
		_fov->getValue(time, fov, validityInterval);
	else
		_zoom->getValue(time, fov, validityInterval);
	return fov;
}

/******************************************************************************
* Changes the field of view of the camera.
******************************************************************************/
void CameraObject::setFieldOfView(TimePoint time, FloatType newFOV)
{
	if(isPerspective())
		_fov->setValue(time, newFOV);
	else
		_zoom->setValue(time, newFOV);
}

/******************************************************************************
* Constructor that creates the UI controls for the editor.
******************************************************************************/
void CameraObjectEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	// Create the rollout.
	QWidget* rollout = createRollout(tr("Camera"), rolloutParams);

	QGridLayout* layout = new QGridLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(0);
	layout->setColumnStretch(1, 1);

	// Is perspective parameter.
	BooleanParameterUI* isPerspectivePUI = new BooleanParameterUI(this, PROPERTY_FIELD(CameraObject::_isPerspective));
	layout->addWidget(isPerspectivePUI->checkBox(), 0, 0, 1, 2);

	// FOV parameter.
	FloatParameterUI* fovPUI = new FloatParameterUI(this, PROPERTY_FIELD(CameraObject::_fov));
	layout->addWidget(fovPUI->label(), 1, 0);
	layout->addLayout(fovPUI->createFieldLayout(), 1, 1);
	fovPUI->setMinValue(fovPUI->parameterUnit()->userToNative(0.01f));
	fovPUI->setMaxValue(fovPUI->parameterUnit()->userToNative(179.99f));

	// Zoom parameter.
	FloatParameterUI* zoomPUI = new FloatParameterUI(this, PROPERTY_FIELD(CameraObject::_zoom));
	layout->addWidget(zoomPUI->label(), 2, 0);
	layout->addLayout(zoomPUI->createFieldLayout(), 2, 1);
	zoomPUI->setMinValue(0);
}

};