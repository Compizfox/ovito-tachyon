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
#include "ProgressiveTask.h"

namespace Ovito {

/**
 * \brief Type of Task to be used within a single-thread context.
 *
 * The methods of this class are not thread-safe and may only be from the main thread.
 * A new main-thread task is typically created using the TaskManager::createMainThreadOperation() method.
 *
 * \sa ThreadSafeTask
 */
class OVITO_CORE_EXPORT MainThreadTask : public ProgressiveTask
{
public:

	/// Constructor.
	explicit MainThreadTask(State initialState, TaskManager* taskManager) : ProgressiveTask(initialState, taskManager) {}

	/// Sets the current progress value (must be in the range 0 to progressMaximum()).
	/// Returns false if the promise has been canceled.
    virtual bool setProgressValue(qlonglong value) override;

	/// Increments the progress value by 1.
	/// Returns false if the promise has been canceled.
    virtual bool incrementProgressValue(qlonglong increment = 1) override;

	/// Changes the status text of this promise.
	virtual void setProgressText(const QString& progressText) override;
};

}	// End of namespace
