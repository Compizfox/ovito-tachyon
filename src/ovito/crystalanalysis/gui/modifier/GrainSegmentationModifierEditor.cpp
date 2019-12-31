////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2019 Alexander Stukowski
//  Copyright 2019 Peter Mahler Larsen
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

#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/crystalanalysis/modifier/grains/GrainSegmentationEngine.h>
#include <ovito/crystalanalysis/modifier/grains/GrainSegmentationModifier.h>
#include <ovito/particles/gui/modifier/analysis/StructureListParameterUI.h>
#include <ovito/gui/properties/FloatParameterUI.h>
#include <ovito/gui/properties/IntegerParameterUI.h>
#include <ovito/gui/properties/BooleanParameterUI.h>
#include <ovito/gui/properties/BooleanRadioButtonParameterUI.h>
#include <ovito/gui/utilities/concurrent/ProgressDialog.h>
#include <ovito/core/dataset/DataSetContainer.h>
#include "GrainSegmentationModifierEditor.h"

#include <3rdparty/qwt/qwt_plot_zoneitem.h>

namespace Ovito { namespace CrystalAnalysis {

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

	// TODO: display the label text
	BooleanRadioButtonParameterUI* algorithmTypeUI = new BooleanRadioButtonParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::algorithmType));
	algorithmTypeUI->buttonFalse()->setText(tr("Node Pair Sampling"));
	algorithmTypeUI->buttonTrue()->setText(tr("Min. Span. Tree"));
	//sublayout2->addWidget(algorithmTypeUI->label(), 1, 0);
	sublayout2->addWidget(algorithmTypeUI->buttonFalse(), 1, 1, 1, -1);
	sublayout2->addWidget(algorithmTypeUI->buttonTrue(), 2, 1, 1, -1);
	algorithmTypeUI->setEnabled(true);

	FloatParameterUI* mergingThresholdUI = new FloatParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::mergingThreshold));
	sublayout2->addWidget(mergingThresholdUI->label(), 3, 0);
	sublayout2->addLayout(mergingThresholdUI->createFieldLayout(), 3, 1);

	IntegerParameterUI* minGrainAtomCountUI = new IntegerParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::minGrainAtomCount));
	sublayout2->addWidget(minGrainAtomCountUI->label(), 4, 0);
	sublayout2->addLayout(minGrainAtomCountUI->createFieldLayout(), 4, 1);

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

	// Create plot widget for merge distances
	_mergePlotWidget = new DataSeriesPlotWidget();
	_mergePlotWidget->setMinimumHeight(200);
	_mergePlotWidget->setMaximumHeight(200);
	_mergeRangeIndicator = new QwtPlotZoneItem();
	_mergeRangeIndicator->setOrientation(Qt::Vertical);
	_mergeRangeIndicator->setZ(1);
	_mergeRangeIndicator->attach(_mergePlotWidget);
	_mergeRangeIndicator->hide();
	layout->addSpacing(10);
	layout->addWidget(_mergePlotWidget);
	connect(this, &GrainSegmentationModifierEditor::contentsReplaced, this, &GrainSegmentationModifierEditor::plotMerges);
}

/******************************************************************************
* This method is called when a reference target changes.
******************************************************************************/
bool GrainSegmentationModifierEditor::referenceEvent(RefTarget* source, const ReferenceEvent& event)
{
	if(source == modifierApplication() && event.type() == ReferenceEvent::PipelineCacheUpdated) {
		plotHistogramLater(this);
		plotLater(this);
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

void GrainSegmentationModifierEditor::plotMerges()
{
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(editObject());

	bool activeIndicator = false;
	if(modifier) {
		// Request the modifier's pipeline output.
		const PipelineFlowState& state = getModifierOutput();
		auto mergeSeries = state.getObjectBy<DataSeriesObject>(modifierApplication(), QStringLiteral("grains-merge"));
		if (mergeSeries) {
			ConstPropertyAccess<FloatType> XProperty = mergeSeries->getX();
			if (XProperty && XProperty.size() > 0) {
				activeIndicator = true;
				_mergeRangeIndicator->setInterval(XProperty[0],	// TODO: this would look nicer if it started at the beginning of the window rather than the first data point
									modifier->mergingThreshold());
				_mergeRangeIndicator->show();
			}
		}
	}

	if (!activeIndicator) {
		_mergeRangeIndicator->hide();
	}

	if(modifierApplication()) {
		// Request the modifier's pipeline output.
		const PipelineFlowState& state = getModifierOutput();

		// Look up the data series in the modifier's pipeline output.
		_mergePlotWidget->setSeries(state.getObjectBy<DataSeriesObject>(modifierApplication(), QStringLiteral("grains-merge")));
	}
	else {
		_mergePlotWidget->reset();
	}
}

}	// End of namespace
}	// End of namespace
