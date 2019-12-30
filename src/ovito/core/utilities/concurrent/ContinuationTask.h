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


#include <ovito/core/Core.h>
#include "Task.h"
#include "FutureDetail.h"

namespace Ovito { OVITO_BEGIN_INLINE_NAMESPACE(Util) OVITO_BEGIN_INLINE_NAMESPACE(Concurrency)

/******************************************************************************
* This shared state is returned by the Future::then() method.
******************************************************************************/
template<typename tuple_type>
class ContinuationTask : public TaskWithResultStorage<Task, tuple_type>
{
public:

	/// Constructor.
	ContinuationTask(TaskDependency creatorState) :
		TaskWithResultStorage<Task, tuple_type>(typename TaskWithResultStorage<Task, tuple_type>::no_result_init_t()),
		_creatorState(std::move(creatorState)) {}

	/// Cancels this promise.
	virtual void cancel() noexcept override {
		if(!this->isCanceled()) {
			TaskWithResultStorage<Task, tuple_type>::cancel();
			this->setStarted();
			this->setFinished();
		}
	}

	/// Marks this promise as fulfilled.
	virtual void setFinished() override {
		// Our reference to the creator state is no longer needed.
		_creatorState.reset();
		TaskWithResultStorage<Task, tuple_type>::setFinished();
	}

	/// Returns the promise that created this promise as a continuation.
	const TaskPtr& creatorState() const { return _creatorState.get(); }

	template<typename FC, typename Args>
	auto fulfillWith(FC&& cont, Args&& params) noexcept
		-> std::enable_if_t<Ovito::detail::is_void_continuation_func<FC,Args>::value>
	{
		try {
			// Call the continuation function with the results of the fulfilled promise.
			this->setStarted();
			Ovito::detail::apply(std::forward<FC>(cont), std::forward<Args>(params));
			this->setFinished();
		}
		catch(...) {
			this->captureException();
			this->setFinished();
		}
	}

	template<typename FC, typename Args>
	auto fulfillWith(FC&& cont, Args&& params) noexcept
		-> std::enable_if_t<!Ovito::detail::is_void_continuation_func<FC,Args>::value>
	{
		try {
			// Call the continuation function with the results of the fulfilled promise.
			this->setStarted();
			setResultsDirect(Ovito::detail::apply(std::forward<FC>(cont), std::forward<Args>(params)));
			this->setFinished();
		}
		catch(...) {
			this->captureException();
			this->setFinished();
		}
	}

protected:

	/// Assigns a result to this shared state.
	template<typename source_tuple_type>
	auto setResultsDirect(source_tuple_type&& results) -> typename std::enable_if<std::tuple_size<source_tuple_type>::value>::type {
		static_assert(std::tuple_size<tuple_type>::value != 0, "Must not be an empty tuple");
		static_assert(std::is_same<tuple_type, std::decay_t<source_tuple_type>>::value, "Must assign a compatible tuple");
		this->template setResults<tuple_type>(std::forward<source_tuple_type>(results));
	}

	/// Assigns a result to this shared state.
	template<typename value_type>
	void setResultsDirect(value_type&& result) {
		static_assert(std::tuple_size<tuple_type>::value == 1, "Must be a tuple of size 1");
		this->template setResults<tuple_type>(std::forward_as_tuple(std::forward<value_type>(result)));
	}

	/// The shared state that created this shared state as a continuation.
	TaskDependency _creatorState;

	template<typename... R2> friend class Future;
};

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace


