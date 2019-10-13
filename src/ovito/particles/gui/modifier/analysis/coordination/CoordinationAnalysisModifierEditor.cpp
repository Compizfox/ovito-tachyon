////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2019 Alexander Stukowski
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

#include <ovito/particles/gui/ParticlesGui.h>
#include <ovito/particles/modifier/analysis/coordination/CoordinationAnalysisModifier.h>
#include <ovito/gui/properties/IntegerParameterUI.h>
#include <ovito/gui/properties/FloatParameterUI.h>
#include <ovito/gui/properties/BooleanParameterUI.h>
#include <ovito/gui/mainwin/MainWindow.h>
#include "CoordinationAnalysisModifierEditor.h"

namespace Ovito { namespace Particles { OVITO_BEGIN_INLINE_NAMESPACE(Modifiers) OVITO_BEGIN_INLINE_NAMESPACE(Analysis) OVITO_BEGIN_INLINE_NAMESPACE(Internal)

IMPLEMENT_OVITO_CLASS(CoordinationAnalysisModifierEditor);
SET_OVITO_OBJECT_EDITOR(CoordinationAnalysisModifier, CoordinationAnalysisModifierEditor);

/******************************************************************************
* Sets up the UI widgets of the editor.
******************************************************************************/
void CoordinationAnalysisModifierEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	// Create a rollout.
	QWidget* rollout = createRollout(tr("Coordination analysis"), rolloutParams, "particles.modifiers.coordination_analysis.html");

    // Create the rollout contents.
	QVBoxLayout* layout = new QVBoxLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(4);

	QGridLayout* gridlayout = new QGridLayout();
	gridlayout->setContentsMargins(4,4,4,4);
	gridlayout->setColumnStretch(1, 1);

	// Cutoff parameter.
	FloatParameterUI* cutoffRadiusPUI = new FloatParameterUI(this, PROPERTY_FIELD(CoordinationAnalysisModifier::cutoff));
	gridlayout->addWidget(cutoffRadiusPUI->label(), 0, 0);
	gridlayout->addLayout(cutoffRadiusPUI->createFieldLayout(), 0, 1);

	// Number of bins parameter.
	IntegerParameterUI* numBinsPUI = new IntegerParameterUI(this, PROPERTY_FIELD(CoordinationAnalysisModifier::numberOfBins));
	gridlayout->addWidget(numBinsPUI->label(), 1, 0);
	gridlayout->addLayout(numBinsPUI->createFieldLayout(), 1, 1);
	layout->addLayout(gridlayout);

	// Partial RDFs option.
	BooleanParameterUI* partialRdfPUI = new BooleanParameterUI(this, PROPERTY_FIELD(CoordinationAnalysisModifier::computePartialRDF));
	layout->addWidget(partialRdfPUI->checkBox());

	_rdfPlot = new DataSeriesPlotWidget();
	_rdfPlot->setMinimumHeight(200);
	_rdfPlot->setMaximumHeight(200);

	layout->addSpacing(12);
	layout->addWidget(new QLabel(tr("Radial distribution function:")));
	layout->addWidget(_rdfPlot);

	QPushButton* btn = new QPushButton(tr("Show in data inspector"));
	connect(btn, &QPushButton::clicked, this, [this]() {
		if(modifierApplication())
			mainWindow()->openDataInspector(modifierApplication());
	});
	layout->addWidget(btn);

	// Status label.
	layout->addSpacing(6);
	layout->addWidget(statusLabel());

	// Update data plot whenever the modifier has calculated new results.
	connect(this, &ModifierPropertiesEditor::contentsReplaced, this, &CoordinationAnalysisModifierEditor::plotRDF);
	connect(this, &ModifierPropertiesEditor::modifierEvaluated, this, [this]() {
		plotRDFLater(this);
	});
}

/******************************************************************************
* Updates the plot of the RDF computed by the modifier.
******************************************************************************/
void CoordinationAnalysisModifierEditor::plotRDF()
{
	OORef<DataSeriesObject> series;

	if(modifierApplication()) {
		// Look up the data series in the modifier's pipeline output.
		series = getModifierOutput().getObjectBy<DataSeriesObject>(modifierApplication(), QStringLiteral("coordination-rdf"));

		// Determine X plotting range.
		if(series) {
			const auto& rdfY = series->getYStorage();
			double minX = 0;
			for(size_t i = 0; i < rdfY->size(); i++) {
				for(size_t cmpnt = 0; cmpnt < rdfY->componentCount(); cmpnt++) {
					if(rdfY->getFloatComponent(i, cmpnt) != 0) {
						minX = series->getXStorage()->getFloat(i);
						break;
					}
				}
				if(minX) break;
			}
			_rdfPlot->setAxisScale(QwtPlot::xBottom, std::floor(minX * 9.0 / series->intervalEnd()) / 10.0 * series->intervalEnd(), series->intervalEnd());
		}
	}
	_rdfPlot->setSeries(series);
}

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace
}	// End of namespace
