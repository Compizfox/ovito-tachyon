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
#include <ovito/core/utilities/concurrent/AsynchronousTask.h>
#include "Modifier.h"

namespace Ovito {

/**
 * \brief Base class for modifiers that compute their results in a background thread.
 */
class OVITO_CORE_EXPORT AsynchronousModifier : public Modifier
{
	Q_OBJECT
	OVITO_CLASS(AsynchronousModifier)

public:

	/**
	 * Abstract base class for algorithm engines performing the modifier's computation in a background thread.
	 */
	class OVITO_CORE_EXPORT Engine : public AsynchronousTaskBase
	{
	public:

		/// Constructor.
		explicit Engine(const TimeInterval& validityInterval = TimeInterval::infinite()) :
			_validityInterval(validityInterval) {}

#ifdef Q_OS_LINUX
		/// Destructor.
		virtual ~Engine();
#endif

		/// Injects the computed results into the data pipeline.
		virtual void applyResults(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state) = 0;

		/// This method is called by the system whenever a parameter of the modifier changes.
		/// The method can be overriden by subclasses to indicate to the caller whether the engine object should be 
		/// discarded (false) or may be kept in the cache, because the computation results are not affected by the changing parameter (true). 
		virtual bool modifierChanged(const PropertyFieldEvent& event) { return false; }

		/// This method is called by the system whenever the preliminary pipeline input changes.
		/// The method should indicate to the caller whether the cached engine object can be 
		/// kept around in a transient phase until a full evaluation is started (return true) or 
		/// should rather be immediately discarded (return false).
		virtual bool pipelineInputChanged() { return true; }

		/// Creates another engine that performs the next stage of the computation. 
		virtual std::shared_ptr<Engine> createContinuationEngine(ModifierApplication* modApp, const PipelineFlowState& input) { return {}; }

		/// Returns the validity interval of the stored computation results.
		const TimeInterval& validityInterval() const { return _validityInterval; }

		/// Changes the validity interval of the computation results.
		void setValidityInterval(const TimeInterval& iv) { _validityInterval = iv; }

	private:

		/// The validity time interval of the stored computation results.
		TimeInterval _validityInterval;
	};

	/// A managed pointer to an Engine instance.
	using EnginePtr = std::shared_ptr<Engine>;

public:

	/// Constructor.
	AsynchronousModifier(DataSet* dataset);

	/// Modifies the input data synchronously.
	virtual void evaluateSynchronous(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state) override;

	/// Suppress preliminary viewport updates when a parameter of the asynchronous modifier changes.
	virtual bool performPreliminaryUpdateAfterChange() override { return false; }

protected:

	/// Saves the class' contents to the given stream.
	virtual void saveToStream(ObjectSaveStream& stream, bool excludeRecomputableData) override;
	
	/// Loads the class' contents from the given stream.
	virtual void loadFromStream(ObjectLoadStream& stream) override;

	/// Asks the object for the result of the data pipeline.
	virtual Future<PipelineFlowState> evaluate(const PipelineEvaluationRequest& request, ModifierApplication* modApp, const PipelineFlowState& input) override;

	/// Creates a computation engine that will compute the modifier's results.
	virtual Future<EnginePtr> createEngine(const PipelineEvaluationRequest& request, ModifierApplication* modApp, const PipelineFlowState& input) = 0;
};

// Export this class template specialization from the DLL under Windows.
extern template class OVITO_CORE_EXPORT Future<AsynchronousModifier::EnginePtr>;

}	// End of namespace
