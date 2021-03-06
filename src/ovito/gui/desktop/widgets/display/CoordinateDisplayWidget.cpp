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

#include <ovito/gui/desktop/GUI.h>
#include <ovito/gui/desktop/widgets/general/SpinnerWidget.h>
#include <ovito/core/dataset/DataSet.h>
#include <ovito/core/dataset/DataSetContainer.h>
#include <ovito/core/dataset/UndoStack.h>
#include <ovito/core/viewport/ViewportConfiguration.h>
#include "CoordinateDisplayWidget.h"

namespace Ovito {

/******************************************************************************
* Constructor.
******************************************************************************/
CoordinateDisplayWidget::CoordinateDisplayWidget(DataSetContainer& datasetContainer, QWidget* parent) : QFrame(parent), _datasetContainer(datasetContainer)
{
	//setFrameStyle(QFrame::Panel | QFrame::Sunken);
	QHBoxLayout* layout = new QHBoxLayout(this);
	layout->setContentsMargins(2,0,2,0);
	layout->setSpacing(0);
	setEnabled(false);
	hide();

	QLabel* xlabel = new QLabel(tr("X:"), this);
	QLabel* ylabel = new QLabel(tr("Y:"), this);
	QLabel* zlabel = new QLabel(tr("Z:"), this);

	class ShortLineEdit : public QLineEdit {
	public:
		ShortLineEdit(QWidget* parent) : QLineEdit(parent) {}
		virtual QSize sizeHint() const override { return QSize(70, QLineEdit::sizeHint().height()); }
	};

	QLineEdit* xedit = new ShortLineEdit(this);
	QLineEdit* yedit = new ShortLineEdit(this);
	QLineEdit* zedit = new ShortLineEdit(this);

	_spinners[0] = new SpinnerWidget(this, xedit);
	_spinners[1] = new SpinnerWidget(this, yedit);
	_spinners[2] = new SpinnerWidget(this, zedit);

	layout->addWidget(xlabel);
	layout->addWidget(xedit, 1);
	layout->addWidget(_spinners[0]);
	layout->addSpacing(6);
	layout->addWidget(ylabel);
	layout->addWidget(yedit, 1);
	layout->addWidget(_spinners[1]);
	layout->addSpacing(6);
	layout->addWidget(zlabel);
	layout->addWidget(zedit, 1);
	layout->addWidget(_spinners[2]);

	connect(_spinners[0], &SpinnerWidget::spinnerValueChanged, this, &CoordinateDisplayWidget::onSpinnerValueChanged);
	connect(_spinners[1], &SpinnerWidget::spinnerValueChanged, this, &CoordinateDisplayWidget::onSpinnerValueChanged);
	connect(_spinners[2], &SpinnerWidget::spinnerValueChanged, this, &CoordinateDisplayWidget::onSpinnerValueChanged);
	connect(_spinners[0], &SpinnerWidget::spinnerDragStart, this, &CoordinateDisplayWidget::onSpinnerDragStart);
	connect(_spinners[1], &SpinnerWidget::spinnerDragStart, this, &CoordinateDisplayWidget::onSpinnerDragStart);
	connect(_spinners[2], &SpinnerWidget::spinnerDragStart, this, &CoordinateDisplayWidget::onSpinnerDragStart);
	connect(_spinners[0], &SpinnerWidget::spinnerDragStop, this, &CoordinateDisplayWidget::onSpinnerDragStop);
	connect(_spinners[1], &SpinnerWidget::spinnerDragStop, this, &CoordinateDisplayWidget::onSpinnerDragStop);
	connect(_spinners[2], &SpinnerWidget::spinnerDragStop, this, &CoordinateDisplayWidget::onSpinnerDragStop);
	connect(_spinners[0], &SpinnerWidget::spinnerDragAbort, this, &CoordinateDisplayWidget::onSpinnerDragAbort);
	connect(_spinners[1], &SpinnerWidget::spinnerDragAbort, this, &CoordinateDisplayWidget::onSpinnerDragAbort);
	connect(_spinners[2], &SpinnerWidget::spinnerDragAbort, this, &CoordinateDisplayWidget::onSpinnerDragAbort);

	QToolButton* animateButton = new QToolButton(this);
	animateButton->setText(tr("A"));
	animateButton->setFocusPolicy(Qt::NoFocus);
	animateButton->setAutoRaise(true);
	animateButton->setToolButtonStyle(Qt::ToolButtonTextOnly);
	animateButton->setToolTip(tr("Animate transformation..."));
	layout->addSpacing(6);
	layout->addWidget(animateButton);
	connect(animateButton, &QAbstractButton::clicked, this, &CoordinateDisplayWidget::animatePressed);
}

/******************************************************************************
* Shows the coordinate display widget.
******************************************************************************/
void CoordinateDisplayWidget::activate(const QString& undoOperationName)
{
	setEnabled(true);
	_undoOperationName = undoOperationName;
	show();
}

/******************************************************************************
* Deactivates the coordinate display widget.
******************************************************************************/
void CoordinateDisplayWidget::deactivate()
{
	if(isEnabled()) {
		setEnabled(false);
		hide();
		_spinners[0]->setFloatValue(0);
		_spinners[1]->setFloatValue(0);
		_spinners[2]->setFloatValue(0);
	}
}

/******************************************************************************
* Is called when a spinner value has been changed by the user.
******************************************************************************/
void CoordinateDisplayWidget::onSpinnerValueChanged()
{
	if(DataSet* dataset = _datasetContainer.currentSet()) {
		int component;
		if(sender() == _spinners[0]) component = 0;
		else if(sender() == _spinners[1]) component = 1;
		else if(sender() == _spinners[2]) component = 2;
		else return;
		ViewportSuspender noVPUpdate(dataset->viewportConfig());
		if(!dataset->undoStack().isRecording()) {
			UndoableTransaction transaction(dataset->undoStack(), _undoOperationName);
			Q_EMIT valueEntered(component, _spinners[component]->floatValue());
			transaction.commit();
		}
		else {
			dataset->undoStack().resetCurrentCompoundOperation();
			Q_EMIT valueEntered(component, _spinners[component]->floatValue());
		}
	}
}

/******************************************************************************
* Is called when the user has started a spinner drag operation.
******************************************************************************/
void CoordinateDisplayWidget::onSpinnerDragStart()
{
	if(DataSet* dataset = _datasetContainer.currentSet()) {
		dataset->undoStack().beginCompoundOperation(_undoOperationName);
	}
}

/******************************************************************************
* Is called when the user has finished the spinner drag operation.
******************************************************************************/
void CoordinateDisplayWidget::onSpinnerDragStop()
{
	if(DataSet* dataset = _datasetContainer.currentSet()) {
		dataset->undoStack().endCompoundOperation();
	}
}

/******************************************************************************
* Is called when the user has aborted the spinner drag operation.
******************************************************************************/
void CoordinateDisplayWidget::onSpinnerDragAbort()
{
	if(DataSet* dataset = _datasetContainer.currentSet()) {
		dataset->undoStack().endCompoundOperation(false);
	}
}

}	// End of namespace
