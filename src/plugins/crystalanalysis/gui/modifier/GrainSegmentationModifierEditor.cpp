///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2018) Alexander Stukowski
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

#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <plugins/particles/gui/modifier/analysis/StructureListParameterUI.h>
#include <gui/properties/FloatParameterUI.h>
#include <gui/properties/IntegerParameterUI.h>
#include <gui/properties/BooleanParameterUI.h>
#include <gui/properties/VariantComboBoxParameterUI.h>
#include <gui/properties/SubObjectParameterUI.h>
#include <gui/properties/BooleanGroupBoxParameterUI.h>
#include <gui/utilities/concurrent/ProgressDialog.h>
#include <core/dataset/DataSetContainer.h>
#include <plugins/crystalanalysis/modifier/grains/GrainSegmentationEngine.h>
#include <plugins/crystalanalysis/modifier/grains/GrainSegmentationModifier.h>
#include "GrainSegmentationModifierEditor.h"

#include <3rdparty/qwt/qwt_plot.h>
#include <3rdparty/qwt/qwt_plot_curve.h>
#include <3rdparty/qwt/qwt_plot_zoneitem.h>
#include <3rdparty/qwt/qwt_plot_grid.h>


namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

IMPLEMENT_OVITO_CLASS(GrainSegmentationModifierEditor);
SET_OVITO_OBJECT_EDITOR(GrainSegmentationModifier, GrainSegmentationModifierEditor);

/******************************************************************************
* Sets up the UI widgets of the editor.
******************************************************************************/
void GrainSegmentationModifierEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	// Create the rollout.
	QWidget* rollout = createRollout(tr("Grain segmentation"), rolloutParams);

	QVBoxLayout* layout = new QVBoxLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(6);

	QGroupBox* paramsBox = new QGroupBox(tr("Parameters"));
	layout->addWidget(paramsBox);
	QGridLayout* sublayout2 = new QGridLayout(paramsBox);
	sublayout2->setContentsMargins(4,4,4,4);
	sublayout2->setSpacing(4);
	sublayout2->setColumnStretch(1, 1);

	FloatParameterUI* rmsdCutoffUI = new FloatParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::rmsdCutoff));
	sublayout2->addWidget(rmsdCutoffUI->label(), 0, 0);
	sublayout2->addLayout(rmsdCutoffUI->createFieldLayout(), 0, 1);

	FloatParameterUI* mergingThresholdUI = new FloatParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::mergingThreshold));
	sublayout2->addWidget(mergingThresholdUI->label(), 1, 0);
	sublayout2->addLayout(mergingThresholdUI->createFieldLayout(), 1, 1);

	IntegerParameterUI* minGrainAtomCountUI = new IntegerParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::minGrainAtomCount));
	sublayout2->addWidget(minGrainAtomCountUI->label(), 2, 0);
	sublayout2->addLayout(minGrainAtomCountUI->createFieldLayout(), 2, 1);

	QPushButton* grainTrackingButton = new QPushButton(tr("Perform grain tracking..."));
	layout->addWidget(grainTrackingButton);
	connect(grainTrackingButton, &QPushButton::clicked, this, &GrainSegmentationModifierEditor::onPerformGrainTracking);

	// Status label.
	layout->addWidget(statusLabel());

	// Structure list.
	StructureListParameterUI* structureTypesPUI = new StructureListParameterUI(this, true);
	layout->addSpacing(10);
	layout->addWidget(new QLabel(tr("Structure types:")));
	layout->addWidget(structureTypesPUI->tableWidget());

	_plot = new QwtPlot();
	_plot->setMinimumHeight(240);
	_plot->setMaximumHeight(240);
	_plot->setCanvasBackground(Qt::white);
	_plot->setAxisTitle(QwtPlot::xBottom, tr("RMSD"));	
	_plot->setAxisTitle(QwtPlot::yLeft, tr("Count"));	

	layout->addSpacing(10);
	layout->addWidget(_plot);
	connect(this, &GrainSegmentationModifierEditor::contentsReplaced, this, &GrainSegmentationModifierEditor::plotHistogram);
}

/******************************************************************************
* This method is called when a reference target changes.
******************************************************************************/
bool GrainSegmentationModifierEditor::referenceEvent(RefTarget* source, const ReferenceEvent& event)
{
	if(event.sender() == editObject() && (event.type() == ReferenceEvent::ObjectStatusChanged || event.type() == ReferenceEvent::TargetChanged)) {
		plotHistogramLater(this);
	}
	return ModifierPropertiesEditor::referenceEvent(source, event);
}

/******************************************************************************
* Replots the histogram computed by the modifier.
******************************************************************************/
void GrainSegmentationModifierEditor::plotHistogram()
{
#if 0
//TODO: put this back in
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(editObject());
	GrainSegmentationModifierApplication* modApp = dynamic_object_cast<GrainSegmentationModifierApplication>(someModifierApplication());
	
	if(!modifier || !modApp || modApp->rmsdHistogramData().empty()) {
		if(_plotCurve) _plotCurve->hide();
		return;
	}

	QVector<QPointF> plotData(modApp->rmsdHistogramData().size());
	double binSize = modApp->rmsdHistogramBinSize();
	double maxHistogramData = 0;
	for(int i = 0; i < modApp->rmsdHistogramData().size(); i++) {
		plotData[i].rx() = binSize * ((double)i + 0.5);
		plotData[i].ry() = modApp->rmsdHistogramData()[i];
		maxHistogramData = std::max(maxHistogramData, plotData[i].y());
	}

	if(!_plotCurve) {
		_plotCurve = new QwtPlotCurve();
	    _plotCurve->setRenderHint(QwtPlotItem::RenderAntialiased, true);
		_plotCurve->setBrush(QColor(255, 160, 100));
		_plotCurve->attach(_plot);
		QwtPlotGrid* plotGrid = new QwtPlotGrid();
		plotGrid->setPen(Qt::gray, 0, Qt::DotLine);
		plotGrid->attach(_plot);
	}
	_plotCurve->setSamples(plotData);

	if(modifier->rmsdCutoff() > 0) {
		if(!_rmsdRange) {
			_rmsdRange = new QwtPlotZoneItem();
			_rmsdRange->setOrientation(Qt::Vertical);
			_rmsdRange->setZ(_plotCurve->z() + 1);
			_rmsdRange->attach(_plot);
		}
		_rmsdRange->show();
		_rmsdRange->setInterval(0, modifier->rmsdCutoff());
	}
	else if(_rmsdRange) {
		_rmsdRange->hide();
	}

	_plot->replot();
#endif
}

/******************************************************************************
* Is called when the user clicks the 'Perform grain tracking' button.
******************************************************************************/
void GrainSegmentationModifierEditor::onPerformGrainTracking()
{
#if 0
//TODO: put this back in
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(editObject());
	GrainSegmentationModifierApplication* modApp = dynamic_object_cast<GrainSegmentationModifierApplication>(someModifierApplication());
	if(!modifier || !modApp) return;

	undoableTransaction(tr("Grain tracking"), [this,modifier,modApp]() {
		ProgressDialog progressDialog(container(), modifier->dataset()->container()->taskManager(), tr("Grain tracking"));
		modifier->trackGrains(progressDialog.taskManager(), modApp);
	});
#endif
}

}	// End of namespace
}	// End of namespace
}	// End of namespace
