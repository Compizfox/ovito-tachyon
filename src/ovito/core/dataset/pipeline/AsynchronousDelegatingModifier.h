////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2020 Alexander Stukowski
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


#include <ovito/core/Core.h>
#include <ovito/core/dataset/pipeline/AsynchronousModifier.h>
#include <ovito/core/dataset/pipeline/DelegatingModifier.h>

namespace Ovito {

/**
 * \brief Base class for modifiers that delegate work to a ModifierDelegate object.
 */
class OVITO_CORE_EXPORT AsynchronousDelegatingModifier : public AsynchronousModifier
{
public:

	/// The abstract base class of delegates used by this modifier type.
	using DelegateBaseType = ModifierDelegate;

	/// Give this modifier class its own metaclass.
	class OVITO_CORE_EXPORT DelegatingModifierClass : public ModifierClass
	{
	public:

		/// Inherit constructor from base class.
		using ModifierClass::ModifierClass;

		/// Asks the metaclass whether the modifier can be applied to the given input data.
		virtual bool isApplicableTo(const DataCollection& input) const override;

		/// Return the metaclass of delegates for this modifier type.
		virtual const ModifierDelegate::OOMetaClass& delegateMetaclass() const {
			OVITO_ASSERT_MSG(false, "AsynchronousDelegatingModifier::OOMetaClass::delegateMetaclass()",
				qPrintable(QStringLiteral("Delegating modifier class %1 does not define a corresponding delegate metaclass. "
				"You must override the delegateMetaclass() method in the modifier's metaclass.").arg(name())));
			return DelegateBaseType::OOClass();
		}
	};

	OVITO_CLASS_META(AsynchronousDelegatingModifier, DelegatingModifierClass)
	Q_OBJECT

public:

	/// Constructor.
	AsynchronousDelegatingModifier(DataSet* dataset);

	/// \brief Determines the time interval over which a computed pipeline state will remain valid.
	virtual TimeInterval validityInterval(const PipelineEvaluationRequest& request, const ModifierApplication* modApp) const override;

protected:

	/// Creates a default delegate for this modifier.
	/// This should be called from the modifier's constructor.
	void createDefaultModifierDelegate(const OvitoClass& delegateType, const QString& defaultDelegateTypeName);

protected:

	/// The modifier's delegate.
	DECLARE_MODIFIABLE_REFERENCE_FIELD_FLAGS(ModifierDelegate, delegate, setDelegate, PROPERTY_FIELD_ALWAYS_CLONE);
};

}	// End of namespace
