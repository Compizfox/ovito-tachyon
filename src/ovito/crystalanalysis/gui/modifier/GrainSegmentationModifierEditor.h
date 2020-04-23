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

#pragma once


#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/stdobj/gui/widgets/DataTablePlotWidget.h>
#include <ovito/gui/desktop/properties/ModifierPropertiesEditor.h>
#include <ovito/core/utilities/DeferredMethodInvocation.h>

class QwtPlotZoneItem;

namespace Ovito { namespace CrystalAnalysis {

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

	/// Replots the merge sequence computed by the modifier.
	void plotMerges();

protected:

	/// Creates the user interface controls for the editor.
	virtual void createUI(const RolloutInsertionParameters& rolloutParams) override;

	/// This method is called when a reference target changes.
	virtual bool referenceEvent(RefTarget* source, const ReferenceEvent& event) override;

private:

	/// The graph widget to display the merge size scatter plot.
	DataTablePlotWidget* _mergePlotWidget;

	/// Marks the merge distance cutoff in the scatter plot.
	QwtPlotZoneItem* _mergeRangeIndicator;

	/// The graph widget to display the log-log scatter plot.
	DataTablePlotWidget* _logPlotWidget;

	/// Marks the merge distance cutoff in the log-log scatter plot.
	QwtPlotZoneItem* _logRangeIndicator;

	/// For deferred invocation of the plot repaint function.
	DeferredMethodInvocation<GrainSegmentationModifierEditor, &GrainSegmentationModifierEditor::plotMerges> plotLater;
};

}	// End of namespace
}	// End of namespace
