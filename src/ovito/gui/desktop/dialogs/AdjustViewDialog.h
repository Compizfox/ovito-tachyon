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

#pragma once


#include <ovito/gui/desktop/GUI.h>
#include <ovito/core/viewport/Viewport.h>
#include <ovito/core/oo/RefTargetListener.h>

namespace Ovito {

/**
 * This dialog box lets the user adjust the camera settings of the current viewport.
 */
class AdjustViewDialog : public QDockWidget
{
	Q_OBJECT

public:

	/// Constructor.
	AdjustViewDialog(Viewport* viewport, QWidget* parentWindow);

private Q_SLOTS:

	/// Event handler for the Cancel button.
	void onCancel();

	/// Is called when the user has changed the camera settings.
	void onAdjustCamera();

	/// Updates the values displayed in the dialog.
	void updateGUI();

private:

	bool _isUpdatingGUI = false;

	QRadioButton* _camPerspective;
	QRadioButton* _camParallel;

	SpinnerWidget* _camPosXSpinner;
	SpinnerWidget* _camPosYSpinner;
	SpinnerWidget* _camPosZSpinner;

	SpinnerWidget* _camDirXSpinner;
	SpinnerWidget* _camDirYSpinner;
	SpinnerWidget* _camDirZSpinner;

	SpinnerWidget* _upDirXSpinner;
	SpinnerWidget* _upDirYSpinner;
	SpinnerWidget* _upDirZSpinner;

	SpinnerWidget* _camFOVAngleSpinner;
	SpinnerWidget* _camFOVSpinner;

	RefTargetListener<Viewport> _viewportListener;
	Viewport::ViewType _oldViewType;
	AffineTransformation _oldCameraTM;
	FloatType _oldFOV;
};

}	// End of namespace
