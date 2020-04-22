////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2020 Alexander Stukowski
//  Copyright 2020 Peter Mahler Larsen
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
#include <ovito/gui/desktop/properties/FloatParameterUI.h>
#include <ovito/gui/desktop/properties/IntegerParameterUI.h>
#include <ovito/gui/desktop/properties/BooleanParameterUI.h>
#include <ovito/gui/desktop/properties/IntegerRadioButtonParameterUI.h>
#include <ovito/gui/desktop/mainwin/MainWindow.h>
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
	QWidget* rollout = createRollout(tr("Grain segmentation"), rolloutParams, "particles.modifiers.grain_segmentation.html");

	QVBoxLayout* layout = new QVBoxLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(6);

	QGroupBox* paramsBox = new QGroupBox(tr("Parameters"));
	layout->addWidget(paramsBox);
	QGridLayout* sublayout2 = new QGridLayout(paramsBox);
	sublayout2->setContentsMargins(4,4,4,4);
	sublayout2->setSpacing(4);
	sublayout2->setColumnStretch(1, 1);

	IntegerRadioButtonParameterUI* algorithmTypeUI = new IntegerRadioButtonParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::mergeAlgorithm));
	QGridLayout* sublayout3 = new QGridLayout();
	sublayout3->setContentsMargins(0,0,0,0);
	sublayout3->setSpacing(4);
	sublayout2->setColumnStretch(1, 1);
	sublayout3->addWidget(new QLabel(tr("Algorithm:")), 0, 0);
	QRadioButton* automaticModeButton = algorithmTypeUI->addRadioButton(GrainSegmentationModifier::NodePairSamplingAutomatic, tr("Node Pair Sampling (automatic)"));
	sublayout3->addWidget(automaticModeButton, 0, 1);
	sublayout3->addWidget(algorithmTypeUI->addRadioButton(GrainSegmentationModifier::NodePairSamplingManual, tr("Node Pair Sampling (manual)")), 1, 1);
	sublayout3->addWidget(algorithmTypeUI->addRadioButton(GrainSegmentationModifier::MinimumSpanningTree, tr("Minimum Spanning Tree")), 2, 1);
	sublayout2->addLayout(sublayout3, 0, 0, 1, 2);

	FloatParameterUI* mergingThresholdUI = new FloatParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::mergingThreshold));
	sublayout2->addWidget(mergingThresholdUI->label(), 1, 0);
	sublayout2->addLayout(mergingThresholdUI->createFieldLayout(), 1, 1);
	connect(automaticModeButton, &QAbstractButton::toggled, mergingThresholdUI, &ParameterUI::setDisabled);

	IntegerParameterUI* minGrainAtomCountUI = new IntegerParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::minGrainAtomCount));
	sublayout2->addWidget(minGrainAtomCountUI->label(), 2, 0);
	sublayout2->addLayout(minGrainAtomCountUI->createFieldLayout(), 2, 1);

	QGroupBox* optionsBox = new QGroupBox(tr("Options"));
	layout->addWidget(optionsBox);
	sublayout2 = new QGridLayout(optionsBox);
	sublayout2->setContentsMargins(4,4,4,4);
	sublayout2->setSpacing(4);
	sublayout2->setColumnStretch(1, 1);

	// Orphan atom adoption
	BooleanParameterUI* orphanAdoptionUI = new BooleanParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::orphanAdoption));
	sublayout2->addWidget(orphanAdoptionUI->checkBox(), 0, 0, 1, 2);

	BooleanParameterUI* colorParticlesByGrainUI = new BooleanParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::colorParticlesByGrain));
	sublayout2->addWidget(colorParticlesByGrainUI->checkBox(), 1, 0, 1, 2);

	BooleanParameterUI* outputBondsUI = new BooleanParameterUI(this, PROPERTY_FIELD(GrainSegmentationModifier::outputBonds));
	sublayout2->addWidget(outputBondsUI->checkBox(), 2, 0, 1, 2);

	// Status label.
	layout->addWidget(statusLabel());

	QPushButton* btn = new QPushButton(tr("Show list of grains"));
	connect(btn, &QPushButton::clicked, this, [this]() {
		if(modifierApplication())
			mainWindow()->openDataInspector(modifierApplication(), QStringLiteral("grains"), 1); // Note: Mode hint "1" switches to table view.
	});
	layout->addWidget(btn);

	// Create plot widget for merge distances
	_mergePlotWidget = new DataTablePlotWidget();
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

	// Create plot widget for log distances
	_logPlotWidget = new DataTablePlotWidget();
	_logPlotWidget->setMinimumHeight(200);
	_logPlotWidget->setMaximumHeight(200);
	_logRangeIndicator = new QwtPlotZoneItem();
	_logRangeIndicator->setOrientation(Qt::Vertical);
	_logRangeIndicator->setZ(1);
	_logRangeIndicator->attach(_logPlotWidget);
	_logRangeIndicator->hide();
	layout->addSpacing(10);
	layout->addWidget(_logPlotWidget);
	connect(this, &GrainSegmentationModifierEditor::contentsReplaced, this, &GrainSegmentationModifierEditor::plotMerges);
}

/******************************************************************************
* This method is called when a reference target changes.
******************************************************************************/
bool GrainSegmentationModifierEditor::referenceEvent(RefTarget* source, const ReferenceEvent& event)
{
	if(source == modifierApplication() && event.type() == ReferenceEvent::PipelineCacheUpdated) {
		plotLater(this);
	}
	return ModifierPropertiesEditor::referenceEvent(source, event);
}

/******************************************************************************
* Replots the merge sequence computed by the modifier.
******************************************************************************/
void GrainSegmentationModifierEditor::plotMerges()
{
	GrainSegmentationModifier* modifier = static_object_cast<GrainSegmentationModifier>(editObject());

	if(modifier && modifierApplication()) {
		// Request the modifier's pipeline output.
		const PipelineFlowState& state = getModifierOutput();

		// Look up the data table in the modifier's pipeline output.
		_mergePlotWidget->setTable(state.getObjectBy<DataTable>(modifierApplication(), QStringLiteral("grains-merge")));

		// Indicate the current merge threshold in the plot.
		FloatType mergingThreshold = modifier->mergingThreshold();
		if(modifier->mergeAlgorithm() == GrainSegmentationModifier::NodePairSamplingAutomatic) {
			mergingThreshold = state.getAttributeValue(modifierApplication(), QStringLiteral("GrainSegmentation.auto_merge_threshold"), mergingThreshold).value<FloatType>();
		}
		_mergeRangeIndicator->setInterval(std::numeric_limits<double>::lowest(), mergingThreshold);
		_mergeRangeIndicator->show();


		// Look up the data table in the modifier's pipeline output.
		_logPlotWidget->setTable(state.getObjectBy<DataTable>(modifierApplication(), QStringLiteral("grains-log")));

		// Indicate the current log threshold in the plot.
		//FloatType mergingThreshold = modifier->mergingThreshold();
		//if(modifier->mergeAlgorithm() == GrainSegmentationModifier::NodePairSamplingAutomatic) {
		//	mergingThreshold = state.getAttributeValue(modifierApplication(), QStringLiteral("GrainSegmentation.auto_merge_threshold"), mergingThreshold).value<FloatType>();
		//}
		//_mergeRangeIndicator->setInterval(std::numeric_limits<double>::lowest(), mergingThreshold);
		//_mergeRangeIndicator->show();
	}
	else {
		_mergePlotWidget->reset();
		_mergeRangeIndicator->hide();

		_logPlotWidget->reset();
		_logRangeIndicator->hide();
	}
}

}	// End of namespace
}	// End of namespace
