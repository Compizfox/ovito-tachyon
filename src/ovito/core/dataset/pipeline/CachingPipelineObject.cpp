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

#include <ovito/core/Core.h>
#include <ovito/core/dataset/DataSet.h>
#include <ovito/core/dataset/animation/AnimationSettings.h>
#include <ovito/core/dataset/pipeline/CachingPipelineObject.h>

namespace Ovito { OVITO_BEGIN_INLINE_NAMESPACE(ObjectSystem) OVITO_BEGIN_INLINE_NAMESPACE(Scene)

IMPLEMENT_OVITO_CLASS(CachingPipelineObject);

/******************************************************************************
* Constructor.
******************************************************************************/
CachingPipelineObject::CachingPipelineObject(DataSet* dataset) : PipelineObject(dataset)
{
}

/******************************************************************************
* Throws away the cached pipeline state.
******************************************************************************/
void CachingPipelineObject::invalidatePipelineCache(TimeInterval keepInterval)
{
	// Reduce the cache validity to the interval to be kept.
	_pipelineCache.invalidate(false, keepInterval);

	// Abort any pipeline evaluation currently in progress unless it
	// falls inside the time interval that should be kept.
	if(!keepInterval.contains(_inProgressEvalTime)) {
		_inProgressEvalFuture.reset();
		_inProgressEvalTime = TimeNegativeInfinity();
	}
}

/******************************************************************************
* Asks the object for the result of the data pipeline.
******************************************************************************/
SharedFuture<PipelineFlowState> CachingPipelineObject::evaluate(TimePoint time, bool breakOnError)
{
	// Check if we can immediately serve the request from the internal cache.
	//
	// Workaround for bug #150: Force FileSource to reload the frame data after the current
	// animation frame has changed, even if it the data is available in the cache.
	if(_pipelineCache.contains(time, time == dataset()->animationSettings()->time()))
		return _pipelineCache.getAt(time);

	// Check if there is already an evaluation in progress whose shared future we can return to the caller.
	if(_inProgressEvalTime == time) {
		SharedFuture<PipelineFlowState> sharedFuture = _inProgressEvalFuture.lock();
		if(sharedFuture.isValid() && !sharedFuture.isCanceled()) {
			return sharedFuture;
		}
	}

	// Let the subclass perform the actual pipeline evaluation.
	Future<PipelineFlowState> stateFuture = evaluateInternal(time, breakOnError);

	// Cache the results in our local pipeline cache.
	if(_pipelineCache.insert(stateFuture, time, this)) {
		// If the cache was updated, we also have a new preliminary state.
		// Inform the pipeline about it.
		if(performPreliminaryUpdateAfterEvaluation() && time == dataset()->animationSettings()->time()) {
			stateFuture = stateFuture.then(executor(), [this](PipelineFlowState&& state) {
				notifyDependents(ReferenceEvent::PreliminaryStateAvailable);
				return std::move(state);
			});
		}
	}
	OVITO_ASSERT(stateFuture.isValid());

	// Keep a weak reference to the future to be able to serve several simultaneous requests.
	SharedFuture<PipelineFlowState> sharedFuture(std::move(stateFuture));
	_inProgressEvalFuture = sharedFuture;
	_inProgressEvalTime = time;

	return sharedFuture;
}

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace
