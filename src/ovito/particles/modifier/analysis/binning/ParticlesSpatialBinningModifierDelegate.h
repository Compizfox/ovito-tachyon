////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2018 Alexander Stukowski
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


#include <ovito/particles/Particles.h>
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/grid/modifier/SpatialBinningModifier.h>

namespace Ovito { namespace Particles { OVITO_BEGIN_INLINE_NAMESPACE(Modifiers) OVITO_BEGIN_INLINE_NAMESPACE(Analysis)

/**
 * \brief Delegate plugin for the SpatialBinningModifier that operates on particles.
 */
class OVITO_PARTICLES_EXPORT ParticlesSpatialBinningModifierDelegate : public SpatialBinningModifierDelegate
{
	/// Give the modifier delegate its own metaclass.
	class OOMetaClass : public SpatialBinningModifierDelegate::OOMetaClass
	{
	public:

		/// Inherit constructor from base class.
		using SpatialBinningModifierDelegate::OOMetaClass::OOMetaClass;

		/// Indicates which data objects in the given input data collection the modifier delegate is able to operate on.
		virtual QVector<DataObjectReference> getApplicableObjects(const DataCollection& input) const override;

		/// Indicates which class of data objects the modifier delegate is able to operate on.
		virtual const DataObject::OOMetaClass& getApplicableObjectClass() const override { return ParticlesObject::OOClass(); }

		/// The name by which Python scripts can refer to this modifier delegate.
		virtual QString pythonDataName() const override { return QStringLiteral("particles"); }
	};

	Q_OBJECT
	OVITO_CLASS_META(ParticlesSpatialBinningModifierDelegate, OOMetaClass);

	Q_CLASSINFO("DisplayName", "Particles");

public:

	/// Constructor.
	Q_INVOKABLE ParticlesSpatialBinningModifierDelegate(DataSet* dataset) : SpatialBinningModifierDelegate(dataset) {}

	/// Creates a computation engine that will perform the actual binning of elements.
	virtual std::shared_ptr<SpatialBinningModifierDelegate::SpatialBinningEngine> createEngine(
				TimePoint time,
				const PipelineFlowState& input,
				const SimulationCell& cell,
				int binDirection,
				ConstPropertyPtr sourceProperty,
				size_t sourceComponent,
				ConstPropertyPtr selectionProperty,
				PropertyPtr binData,
				const Vector3I& binCount,
				const Vector3I& binDir,
				int reductionOperation,
                bool computeFirstDerivative) override;

private:

	/// Asynchronous compute engine that does the actual work in a separate thread.
	class ComputeEngine : public SpatialBinningModifierDelegate::SpatialBinningEngine
	{
	public:

		/// Constructor.
		ComputeEngine(
				const TimeInterval& validityInterval,
				const SimulationCell& cell,
				int binningDirection,
				ConstPropertyPtr sourceProperty,
				size_t sourceComponent,
				ConstPropertyPtr selection,
				ConstPropertyPtr positions,
				PropertyPtr binData,
				const Vector3I& binCount,
				const Vector3I& binDir,
				int reductionOperation,
				bool computeFirstDerivative) :
			SpatialBinningEngine(validityInterval, cell, binningDirection, std::move(sourceProperty), sourceComponent, std::move(selection), std::move(binData), binCount, binDir, reductionOperation, computeFirstDerivative),
			_positions(std::move(positions)) {}

		/// This method is called by the system after the computation was successfully completed.
		virtual void cleanup() override {
			_positions.reset();
			SpatialBinningEngine::cleanup();
		}

		/// Computes the modifier's results.
		virtual void perform() override;

		/// Returns the property storage that contains the input particle positions.
		const ConstPropertyPtr& positions() const { return _positions; }

	private:

		ConstPropertyPtr _positions;
	};
};

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace
}	// End of namespace
