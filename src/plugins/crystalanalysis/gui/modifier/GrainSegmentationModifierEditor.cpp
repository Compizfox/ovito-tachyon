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
#include <plugins/crystalanalysis/modifier/grains/GrainSegmentationEngine.h>
#include <plugins/crystalanalysis/modifier/grains/GrainSegmentationModifier.h>
#include <plugins/particles/gui/modifier/analysis/StructureListParameterUI.h>
#include <gui/properties/FloatParameterUI.h>
#include <gui/properties/IntegerParameterUI.h>
#include <gui/properties/BooleanParameterUI.h>
#include <gui/utilities/concurrent/ProgressDialog.h>
#include <core/dataset/DataSetContainer.h>
#include "GrainSegmentationModifierEditor.h"

#include <3rdparty/qwt/qwt_plot_zoneitem.h>

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

	QGroupBox* debuggingParamsBox = new QGroupBox(tr("Debugging options"));
	layout->addWidget(debuggingParamsBox);
	sublayout2 = new QGridLayout(debuggingParamsBox);
	sublayout2->setContentsMargins(4,4,4,4);
	sublayout2->setSpacing(4);
	sublayout2->setColumnStretch(1, 1);

	// Orphan atom adoption
	BooleanParameterUI* orphanAdoptionUI = new BooleanParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::orphanAdoption));
	sublayout2->addWidget(orphanAdoptionUI->checkBox(), 0, 0, 1, 2);

	BooleanParameterUI* outputBondsUI = new BooleanParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::outputBonds));
	sublayout2->addWidget(outputBondsUI->checkBox(), 1, 0, 1, 2);

#if 0
	QPushButton* grainTrackingButton = new QPushButton(tr("Perform grain tracking..."));
	layout->addWidget(grainTrackingButton);
	connect(grainTrackingButton, &QPushButton::clicked, this, &GrainSegmentationModifierEditor::onPerformGrainTracking);
#endif

	// Status label.
	layout->addWidget(statusLabel());

	// Structure list.
	StructureListParameterUI* structureTypesPUI = new StructureListParameterUI(this, true);
	layout->addSpacing(10);
	layout->addWidget(new QLabel(tr("Structure types:")));
	layout->addWidget(structureTypesPUI->tableWidget());

	// Create plot widget for RMSD distribution.
	_rmsdPlotWidget = new DataSeriesPlotWidget();
	_rmsdPlotWidget->setMinimumHeight(200);
	_rmsdPlotWidget->setMaximumHeight(200);
	_rmsdRangeIndicator = new QwtPlotZoneItem();
	_rmsdRangeIndicator->setOrientation(Qt::Vertical);
	_rmsdRangeIndicator->setZ(1);
	_rmsdRangeIndicator->attach(_rmsdPlotWidget);
	_rmsdRangeIndicator->hide();
	layout->addSpacing(10);
	layout->addWidget(_rmsdPlotWidget);
	connect(this, &GrainSegmentationModifierEditor::contentsReplaced, this, &GrainSegmentationModifierEditor::plotHistogram);
}

/******************************************************************************
* This method is called when a reference target changes.
******************************************************************************/
bool GrainSegmentationModifierEditor::referenceEvent(RefTarget* source, const ReferenceEvent& event)
{
	if(source == modifierApplication() && event.type() == ReferenceEvent::PipelineCacheUpdated) {
		plotHistogramLater(this);
	}
	return ModifierPropertiesEditor::referenceEvent(source, event);
}

/******************************************************************************
* Replots the histogram computed by the modifier.
******************************************************************************/
void GrainSegmentationModifierEditor::plotHistogram()
{
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(editObject());

	if(modifier && modifier->rmsdCutoff() > 0) {
		_rmsdRangeIndicator->setInterval(0, modifier->rmsdCutoff());
		_rmsdRangeIndicator->show();
	}
	else {
		_rmsdRangeIndicator->hide();
	}

	if(modifierApplication()) {
		// Request the modifier's pipeline output.
		const PipelineFlowState& state = getModifierOutput();

		// Look up the data series in the modifier's pipeline output.
		_rmsdPlotWidget->setSeries(state.getObjectBy<DataSeriesObject>(modifierApplication(), QStringLiteral("grains-rmsd")));
	}
	else {
		_rmsdPlotWidget->reset();
	}
}

/******************************************************************************
* Is called when the user clicks the 'Perform grain tracking' button.
******************************************************************************/
void GrainSegmentationModifierEditor::onPerformGrainTracking()
{
#if 0
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(editObject());
	GrainSegmentationModifierApplication* modApp = dynamic_object_cast<GrainSegmentationModifierApplication>(someModifierApplication());
	if(!modifier || !modApp) return;

	undoableTransaction(tr("Grain tracking"), [this,modifier,modApp]() {
		ProgressDialog progressDialog(container(), modifier->dataset()->container()->taskManager(), tr("Grain tracking"));
		modifier->trackGrains(progressDialog.taskManager(), modApp);
	});
#endif
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(editObject());
	if(!modifier) return;

	undoableTransaction(tr("Grain tracking"), [this,modifier]() {
		ProgressDialog progressDialog(container(), modifier->dataset()->container()->taskManager(), tr("Grain tracking"));
	});
}

}	// End of namespace
}	// End of namespace
}	// End of namespace
