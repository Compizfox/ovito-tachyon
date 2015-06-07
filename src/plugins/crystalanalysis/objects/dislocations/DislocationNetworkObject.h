///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2015) Alexander Stukowski
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

#ifndef __OVITO_CA_DISLOCATION_NETWORK_OBJECT_H
#define __OVITO_CA_DISLOCATION_NETWORK_OBJECT_H

#include <plugins/crystalanalysis/CrystalAnalysis.h>
#include <plugins/crystalanalysis/data/DislocationNetwork.h>
#include <core/scene/objects/DataObjectWithSharedStorage.h>
#include <core/gui/properties/PropertiesEditor.h>

namespace Ovito { namespace Plugins { namespace CrystalAnalysis {

/**
 * \brief Stores a collection of dislocation segments.
 */
class OVITO_CRYSTALANALYSIS_EXPORT DislocationNetworkObject : public DataObjectWithSharedStorage<DislocationNetwork>
{
public:

	/// \brief Constructor.
	Q_INVOKABLE DislocationNetworkObject(DataSet* dataset, DislocationNetwork* network = nullptr);

	/// Returns the title of this object.
	virtual QString objectTitle() override { return tr("Dislocations"); }

	/// Returns the list of dislocation segments.
	const std::vector<DislocationSegment*>& segments() const { return storage()->segments(); }

	/// Returns the list of dislocation segments.
	const std::vector<DislocationSegment*>& modifiableSegments() { return modifiableStorage()->segments(); }

private:

	Q_OBJECT
	OVITO_OBJECT
};

/******************************************************************************
* A properties editor for the DislocationNetworkObject class.
******************************************************************************/
class OVITO_CRYSTALANALYSIS_EXPORT DislocationNetworkObjectEditor : public PropertiesEditor
{
public:

	/// Default constructor.
	Q_INVOKABLE DislocationNetworkObjectEditor() {}

protected:

	/// Creates the user interface controls for the editor.
	virtual void createUI(const RolloutInsertionParameters& rolloutParams) override;

protected Q_SLOTS:

	/// Is called when the user presses the "Open Inspector" button.
	void onOpenInspector();

private:

	Q_OBJECT
	OVITO_OBJECT
};

}	// End of namespace
}	// End of namespace
}	// End of namespace

#endif // __OVITO_CA_DISLOCATION_NETWORK_OBJECT_H
