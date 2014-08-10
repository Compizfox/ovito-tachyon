///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2014) Alexander Stukowski
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

#include <plugins/pyscript/PyScript.h>
#include <core/animation/TimeInterval.h>
#include <core/animation/AnimationSettings.h>
#include "PythonBinding.h"

namespace PyScript {

using namespace boost::python;
using namespace Ovito;

void setupAnimationBinding()
{
	class_<TimeInterval>("TimeInterval", init<>())
		.def(init<TimePoint>())
		.def(init<TimePoint, TimePoint>())
		.add_property("start", &TimeInterval::start, &TimeInterval::setStart)
		.add_property("end", &TimeInterval::end, &TimeInterval::setEnd)
		.add_property("isEmpty", &TimeInterval::isEmpty)
		.add_property("isInfinite", &TimeInterval::isInfinite)
		.add_property("duration", &TimeInterval::duration, &TimeInterval::setDuration)
		.def("setInfinite", &TimeInterval::setInfinite)
		.def("setEmpty", &TimeInterval::setEmpty)
		.def("setInstant", &TimeInterval::setInstant)
		.def("contains", &TimeInterval::contains)
		.def("intersect", &TimeInterval::intersect)
		.def("timeToSeconds", &timeToSeconds)
		.def("secondsToTime", &secondsToTime)
		.staticmethod("timeToSeconds")
		.staticmethod("secondsToTime")
		.add_static_property("infinite", &TimeInterval::infinite)
		.add_static_property("empty", &TimeInterval::empty)
		.setattr("TimeNegativeInfinity", TimeNegativeInfinity())
		.setattr("TimePositiveInfinity", TimePositiveInfinity())
		.def(self == TimeInterval())
		.def(self != TimeInterval())
	;

	ovito_class<AnimationSettings, RefTarget>()
		.add_property("time", &AnimationSettings::time, &AnimationSettings::setTime)
		.add_property("animationInterval", make_function(&AnimationSettings::animationInterval, return_value_policy<copy_const_reference>()), &AnimationSettings::setAnimationInterval)
		.add_property("framesPerSecond", &AnimationSettings::framesPerSecond, &AnimationSettings::setFramesPerSecond)
		.add_property("ticksPerFrame", &AnimationSettings::ticksPerFrame, &AnimationSettings::setTicksPerFrame)
		.add_property("currentFrame", &AnimationSettings::currentFrame, &AnimationSettings::setCurrentFrame)
		.add_property("lastFrame", &AnimationSettings::lastFrame, &AnimationSettings::setLastFrame)
		.add_property("firstFrame", &AnimationSettings::firstFrame, &AnimationSettings::setFirstFrame)
		.add_property("playbackSpeed", &AnimationSettings::playbackSpeed, &AnimationSettings::setPlaybackSpeed)
		.add_property("isAnimating", &AnimationSettings::isAnimating)
		.add_property("autoKeyMode", &AnimationSettings::autoKeyMode, &AnimationSettings::setAutoKeyMode)
		.add_property("isTimeChanging", &AnimationSettings::isTimeChanging)
		.def("frameToTime", &AnimationSettings::frameToTime)
		.def("timeToFrame", &AnimationSettings::timeToFrame)
		.def("snapTime", &AnimationSettings::snapTime)
		.def("timeToString", &AnimationSettings::timeToString)
		.def("stringToTime", &AnimationSettings::stringToTime)
		.def("jumpToAnimationStart", &AnimationSettings::jumpToAnimationStart)
		.def("jumpToAnimationEnd", &AnimationSettings::jumpToAnimationEnd)
		.def("jumpToNextFrame", &AnimationSettings::jumpToNextFrame)
		.def("jumpToPreviousFrame", &AnimationSettings::jumpToPreviousFrame)
		.def("startAnimationPlayback", &AnimationSettings::startAnimationPlayback)
		.def("stopAnimationPlayback", &AnimationSettings::stopAnimationPlayback)
	;
}

};
