////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2017 Alexander Stukowski
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
#include <ovito/stdmod/modifiers/ReplicateModifier.h>

namespace Ovito { namespace CrystalAnalysis {

using namespace Ovito::StdMod;

/**
 * \brief Delegate for the ReplicateModifier that operates on dislocation lines.
 */
class DislocationReplicateModifierDelegate : public ReplicateModifierDelegate
{
	/// Give the modifier delegate its own metaclass.
	class OOMetaClass : public ReplicateModifierDelegate::OOMetaClass
	{
	public:

		/// Inherit constructor from base class.
		using ReplicateModifierDelegate::OOMetaClass::OOMetaClass;

		/// Indicates which data objects in the given input data collection the modifier delegate is able to operate on.
		virtual QVector<DataObjectReference> getApplicableObjects(const DataCollection& input) const override;

		/// The name by which Python scripts can refer to this modifier delegate.
		virtual QString pythonDataName() const override { return QStringLiteral("dislocations"); }
	};

	Q_OBJECT
	OVITO_CLASS_META(DislocationReplicateModifierDelegate, OOMetaClass)

	Q_CLASSINFO("DisplayName", "Dislocations");

public:

	/// Constructor.
	Q_INVOKABLE DislocationReplicateModifierDelegate(DataSet* dataset) : ReplicateModifierDelegate(dataset) {}

	/// Applies the modifier operation to the data in a pipeline flow state.
	virtual PipelineStatus apply(Modifier* modifier, PipelineFlowState& state, TimePoint time, ModifierApplication* modApp, const std::vector<std::reference_wrapper<const PipelineFlowState>>& additionalInputs) override;
};

}	// End of namespace
}	// End of namespace
