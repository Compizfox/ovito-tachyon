///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2017) Alexander Stukowski
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


#include <gui/GUI.h>
#include <gui/properties/ModifierPropertiesEditor.h>
#include <plugins/stdobj/gui/widgets/DataSeriesPlotWidget.h>
#include <core/utilities/DeferredMethodInvocation.h>

class QwtPlotMarker;

namespace Ovito { namespace Grid { OVITO_BEGIN_INLINE_NAMESPACE(Internal)

/**
 * \brief A properties editor for the CreateIsosurfaceModifier class.
 */
class CreateIsosurfaceModifierEditor : public ModifierPropertiesEditor
{
	Q_OBJECT
	OVITO_CLASS(CreateIsosurfaceModifierEditor)

public:

	/// Default constructor.
	Q_INVOKABLE CreateIsosurfaceModifierEditor() {}

protected:

	/// Creates the user interface controls for the editor.
	virtual void createUI(const RolloutInsertionParameters& rolloutParams) override;

protected Q_SLOTS:

	/// Replots the value histogram computed by the modifier.
	void plotHistogram();

private:

	/// The graph widget to display the histogram.
	StdObj::DataSeriesPlotWidget* _plotWidget;

	/// The plot item for indicating the current iso level value.
	QwtPlotMarker* _isoLevelIndicator;

	/// For deferred invocation of the plot repaint function.
	DeferredMethodInvocation<CreateIsosurfaceModifierEditor, &CreateIsosurfaceModifierEditor::plotHistogram> plotHistogramLater;
};

OVITO_END_INLINE_NAMESPACE
}	// End of namespace
}	// End of namespace
