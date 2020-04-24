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

#include <ovito/gui/desktop/GUI.h>
#include <ovito/core/viewport/Viewport.h>
#include <ovito/core/utilities/units/UnitsManager.h>
#include <ovito/gui/desktop/widgets/general/SpinnerWidget.h>
#include <ovito/gui/desktop/mainwin/MainWindow.h>
#include "AdjustViewDialog.h"

namespace Ovito {

/******************************************************************************
* The constructor of the camera dialog.
******************************************************************************/
AdjustViewDialog::AdjustViewDialog(Viewport* viewport, QWidget* parent) :
	QDockWidget(tr("Adjust View"), parent)
{
	setFeatures(QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetClosable);
	setAllowedAreas(Qt::NoDockWidgetArea);
	setFloating(true);
	setAttribute(Qt::WA_DeleteOnClose);
	QWidget* widget = new QWidget();
	setWidget(widget);

	_oldViewType = viewport->viewType();
	_oldCameraTM = viewport->cameraTransformation();
	_oldFOV = viewport->fieldOfView();

	_viewportListener.setTarget(viewport);

	connect(&_viewportListener, &RefTargetListenerBase::notificationEvent, this, [this](const ReferenceEvent& event) {
		// Update the values displayed in the dialog when the viewport camera is moved by the user.
		if(event.type() == ReferenceEvent::TargetChanged)
			updateGUI();
		// Close the dialog when the viewport is deleted.
		else if(event.type() == ReferenceEvent::TargetDeleted)
			close();
	});

	QVBoxLayout* mainLayout = new QVBoxLayout(widget);

	QGroupBox* viewPosBox = new QGroupBox(tr("View position"));
	mainLayout->addWidget(viewPosBox);

	QGridLayout* gridLayout = new QGridLayout(viewPosBox);
	gridLayout->setColumnStretch(1,1);
	gridLayout->setColumnStretch(2,1);
	gridLayout->setColumnStretch(3,1);
	gridLayout->addWidget(new QLabel(tr("XYZ:")), 0, 0);

	QHBoxLayout* fieldLayout;
	QLineEdit* textBox;

	_camPosXSpinner = new SpinnerWidget();
	_camPosYSpinner = new SpinnerWidget();
	_camPosZSpinner = new SpinnerWidget();
	_camPosXSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());
	_camPosYSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());
	_camPosZSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	_camPosXSpinner->setTextBox(textBox);
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_camPosXSpinner);
	gridLayout->addLayout(fieldLayout, 0, 1);
	connect(_camPosXSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	_camPosYSpinner->setTextBox(textBox);
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_camPosYSpinner);
	gridLayout->addLayout(fieldLayout, 0, 2);
	connect(_camPosYSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	_camPosZSpinner->setTextBox(textBox);
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_camPosZSpinner);
	gridLayout->addLayout(fieldLayout, 0, 3);
	connect(_camPosZSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	QGroupBox* viewDirBox = new QGroupBox(tr("View direction"));
	mainLayout->addWidget(viewDirBox);

	gridLayout = new QGridLayout(viewDirBox);
	gridLayout->setColumnStretch(1,1);
	gridLayout->setColumnStretch(2,1);
	gridLayout->setColumnStretch(3,1);
	gridLayout->addWidget(new QLabel(tr("XYZ:")), 0, 0);

	_camDirXSpinner = new SpinnerWidget();
	_camDirYSpinner = new SpinnerWidget();
	_camDirZSpinner = new SpinnerWidget();
	_camDirXSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());
	_camDirYSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());
	_camDirZSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	_camDirXSpinner->setTextBox(textBox);
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_camDirXSpinner);
	gridLayout->addLayout(fieldLayout, 0, 1);
	connect(_camDirXSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	_camDirYSpinner->setTextBox(textBox);
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_camDirYSpinner);
	gridLayout->addLayout(fieldLayout, 0, 2);
	connect(_camDirYSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	_camDirZSpinner->setTextBox(textBox);
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_camDirZSpinner);
	gridLayout->addLayout(fieldLayout, 0, 3);
	connect(_camDirZSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	QGroupBox* upDirBox = new QGroupBox(tr("Up direction (read-only)"));
	mainLayout->addWidget(upDirBox);

	gridLayout = new QGridLayout(upDirBox);
	gridLayout->setColumnStretch(1,1);
	gridLayout->setColumnStretch(2,1);
	gridLayout->setColumnStretch(3,1);
	gridLayout->addWidget(new QLabel(tr("XYZ:")), 0, 0);

	_upDirXSpinner = new SpinnerWidget();
	_upDirYSpinner = new SpinnerWidget();
	_upDirZSpinner = new SpinnerWidget();
	_upDirXSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());
	_upDirYSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());
	_upDirZSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	textBox->setReadOnly(true);
	_upDirXSpinner->setTextBox(textBox);
	_upDirXSpinner->hide();
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_upDirXSpinner);
	gridLayout->addLayout(fieldLayout, 0, 1);
//	connect(_upDirXSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	textBox->setReadOnly(true);
	_upDirYSpinner->setTextBox(textBox);
	_upDirYSpinner->hide();
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_upDirYSpinner);
	gridLayout->addLayout(fieldLayout, 0, 2);
//	connect(_upDirYSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	textBox->setReadOnly(true);
	_upDirZSpinner->setTextBox(textBox);
	_upDirZSpinner->hide();
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_upDirZSpinner);
	gridLayout->addLayout(fieldLayout, 0, 3);
//	connect(_camDirZSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	QGroupBox* projectionBox = new QGroupBox(tr("Projection type"));
	mainLayout->addWidget(projectionBox);

	gridLayout = new QGridLayout(projectionBox);
	gridLayout->setColumnMinimumWidth(0, 30);
	gridLayout->setColumnStretch(3, 1);

	_camPerspective = new QRadioButton(tr("Perspective:"));
	connect(_camPerspective, &QRadioButton::clicked, this, &AdjustViewDialog::onAdjustCamera);
	gridLayout->addWidget(_camPerspective, 0, 0, 1, 3);

	gridLayout->addWidget(new QLabel(tr("View angle:")), 1, 1);
	_camFOVAngleSpinner = new SpinnerWidget();
	_camFOVAngleSpinner->setUnit(viewport->dataset()->unitsManager().angleUnit());
	_camFOVAngleSpinner->setMinValue(FloatType(1e-4));
	_camFOVAngleSpinner->setMaxValue(FLOATTYPE_PI - FloatType(1e-2));
	_camFOVAngleSpinner->setFloatValue(qDegreesToRadians(FloatType(35)));
	_camFOVAngleSpinner->setEnabled(false);
	connect(_camPerspective, &QRadioButton::toggled, _camFOVAngleSpinner, &SpinnerWidget::setEnabled);

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	_camFOVAngleSpinner->setTextBox(textBox);
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_camFOVAngleSpinner);
	gridLayout->addLayout(fieldLayout, 1, 2);
	connect(_camFOVAngleSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	_camParallel = new QRadioButton(tr("Parallel:"));
	connect(_camParallel, &QRadioButton::clicked, this, &AdjustViewDialog::onAdjustCamera);
	gridLayout->addWidget(_camParallel, 2, 0, 1, 3);

	gridLayout->addWidget(new QLabel(tr("Field of view:")), 3, 1);
	_camFOVSpinner = new SpinnerWidget();
	_camFOVSpinner->setUnit(viewport->dataset()->unitsManager().worldUnit());
	_camFOVSpinner->setMinValue(FloatType(1e-4));
	_camFOVSpinner->setFloatValue(200);
	_camFOVSpinner->setEnabled(false);
	connect(_camParallel, &QRadioButton::toggled, _camFOVSpinner, &SpinnerWidget::setEnabled);

	fieldLayout = new QHBoxLayout();
	fieldLayout->setContentsMargins(0,0,0,0);
	fieldLayout->setSpacing(0);
	textBox = new QLineEdit();
	_camFOVSpinner->setTextBox(textBox);
	fieldLayout->addWidget(textBox);
	fieldLayout->addWidget(_camFOVSpinner);
	gridLayout->addLayout(fieldLayout, 3, 2);
	connect(_camFOVSpinner, &SpinnerWidget::spinnerValueChanged, this, &AdjustViewDialog::onAdjustCamera);

	QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel | QDialogButtonBox::Help, Qt::Horizontal, this);
	connect(buttonBox, &QDialogButtonBox::accepted, this, &AdjustViewDialog::close);
	connect(buttonBox, &QDialogButtonBox::rejected, this, &AdjustViewDialog::onCancel);
	connect(buttonBox, &QDialogButtonBox::helpRequested, [this]() {
		MainWindow::openHelpTopic("viewports.adjust_view_dialog.html");
	});
	mainLayout->addWidget(buttonBox);

	updateGUI();
}

/******************************************************************************
* Updates the values displayed in the dialog.
******************************************************************************/
void AdjustViewDialog::updateGUI()
{
	_isUpdatingGUI = true;
	Viewport* viewport = _viewportListener.target();

	const Point3& cameraPos = viewport->cameraPosition();
	_camPosXSpinner->setFloatValue(cameraPos.x());
	_camPosYSpinner->setFloatValue(cameraPos.y());
	_camPosZSpinner->setFloatValue(cameraPos.z());

	Vector3 oldCameraDir(_camDirXSpinner->floatValue(), _camDirYSpinner->floatValue(), _camDirZSpinner->floatValue());
	FloatType oldDirLength = oldCameraDir.length();
	if(oldDirLength == 0) oldDirLength = 1;
	const Vector3& cameraDir = viewport->cameraDirection();
	_camDirXSpinner->setFloatValue(cameraDir.x() * oldDirLength);
	_camDirYSpinner->setFloatValue(cameraDir.y() * oldDirLength);
	_camDirZSpinner->setFloatValue(cameraDir.z() * oldDirLength);

	const Vector3& upDir = viewport->cameraTransformation().column(1);
	_upDirXSpinner->setFloatValue(upDir.x());
	_upDirYSpinner->setFloatValue(upDir.y());
	_upDirZSpinner->setFloatValue(upDir.z());

	if(viewport->isPerspectiveProjection()) {
		_camPerspective->setChecked(true);
		_camFOVAngleSpinner->setFloatValue(viewport->fieldOfView());
	}
	else {
		_camParallel->setChecked(true);
		_camFOVSpinner->setFloatValue(viewport->fieldOfView());
	}
	_isUpdatingGUI = false;
}

/******************************************************************************
* Is called when the user has changed the camera settings.
******************************************************************************/
void AdjustViewDialog::onAdjustCamera()
{
	if(_isUpdatingGUI)
		return;
	Viewport* viewport = _viewportListener.target();
	if(_camPerspective->isChecked()) {
		viewport->setViewType(Viewport::VIEW_PERSPECTIVE);
		viewport->setFieldOfView(_camFOVAngleSpinner->floatValue());
	}
	else {
		viewport->setViewType(Viewport::VIEW_ORTHO);
		viewport->setFieldOfView(_camFOVSpinner->floatValue());
	}

	viewport->setCameraPosition(Point3(_camPosXSpinner->floatValue(), _camPosYSpinner->floatValue(), _camPosZSpinner->floatValue()));
	viewport->setCameraDirection(Vector3(_camDirXSpinner->floatValue(), _camDirYSpinner->floatValue(), _camDirZSpinner->floatValue()));
}

/******************************************************************************
* Event handler for the Cancel button.
******************************************************************************/
void AdjustViewDialog::onCancel()
{
	Viewport* viewport = _viewportListener.target();
	viewport->setViewType(_oldViewType);
	viewport->setCameraTransformation(_oldCameraTM);
	viewport->setFieldOfView(_oldFOV);

	close();
}

}	// End of namespace
