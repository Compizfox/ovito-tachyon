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

#pragma once


#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <plugins/stdobj/gui/widgets/DataSeriesPlotWidget.h>
#include <gui/properties/ModifierPropertiesEditor.h>
#include <core/utilities/DeferredMethodInvocation.h>

class QwtPlotZoneItem;

namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

/**
 * Properties editor for the GrainSegmentationModifier class.
 */
class GrainSegmentationModifierEditor : public ModifierPropertiesEditor
{
	Q_OBJECT
	OVITO_CLASS(GrainSegmentationModifierEditor)

public:

	/// Default constructor.
	Q_INVOKABLE GrainSegmentationModifierEditor() {}

protected Q_SLOTS:

	/// Replots the histogram computed by the modifier.
	void plotHistogram();

void plotMerges();

	/// Is called when the user clicks the 'Perform grain tracking' button.
	void onPerformGrainTracking();

protected:

	/// Creates the user interface controls for the editor.
	virtual void createUI(const RolloutInsertionParameters& rolloutParams) override;

	/// This method is called when a reference target changes.
	virtual bool referenceEvent(RefTarget* source, const ReferenceEvent& event) override;

private:

	/// The graph widget to display the RMSD histogram.
	DataSeriesPlotWidget* _rmsdPlotWidget;

	/// Marks the RMSD cutoff in the histogram plot.
	QwtPlotZoneItem* _rmsdRangeIndicator;

	/// For deferred invocation of the plot repaint function.
	DeferredMethodInvocation<GrainSegmentationModifierEditor, &GrainSegmentationModifierEditor::plotHistogram> plotHistogramLater;

	/// The graph widget to display the merge size scatter plot.
	DataSeriesPlotWidget* _mergePlotWidget;

	/// Marks the merge distance cutoff in the scatter plot.
	QwtPlotZoneItem* _mergeRangeIndicator;

	/// For deferred invocation of the plot repaint function.
	DeferredMethodInvocation<GrainSegmentationModifierEditor, &GrainSegmentationModifierEditor::plotMerges> plotLater;
};

}	// End of namespace
}	// End of namespace
}	// End of namespace
