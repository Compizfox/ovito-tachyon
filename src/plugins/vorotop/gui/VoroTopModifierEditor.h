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


#include <plugins/vorotop/VoroTopPlugin.h>
#include <plugins/vorotop/VoroTopModifier.h>
#include <gui/properties/ModifierPropertiesEditor.h>

namespace Ovito { namespace VoroTop {

/**
 * A properties editor for the VoroTopModifier class.
 */
class VoroTopModifierEditor : public ModifierPropertiesEditor
{
public:

	/// Default constructor.
	Q_INVOKABLE VoroTopModifierEditor() {}

protected:

	/// Creates the user interface controls for the editor.
	virtual void createUI(const RolloutInsertionParameters& rolloutParams) override;

private Q_SLOTS:

	/// Is called when the user presses the 'Load filter' button.
	void onLoadFilter();

private:

	Q_OBJECT
	OVITO_CLASS
};

}	// End of namespace
}	// End of namespace
