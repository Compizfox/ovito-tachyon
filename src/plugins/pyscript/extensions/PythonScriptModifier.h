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

#ifndef __OVITO_PYTHON_SCRIPT_MODIFIER_H
#define __OVITO_PYTHON_SCRIPT_MODIFIER_H

#include <plugins/pyscript/PyScript.h>
#include <core/utilities/concurrent/FutureInterface.h>
#include <core/scene/pipeline/Modifier.h>
#include <core/scene/objects/CompoundObject.h>
#include <plugins/pyscript/engine/ScriptEngine.h>

#include <boost/optional.hpp>

namespace PyScript {

using namespace Ovito;

/**
 * \brief A modifier that executes a Python script.
 */
class OVITO_PYSCRIPT_EXPORT PythonScriptModifier : public Modifier
{
public:

	// A helper class that provides the progress callback interface for Python scripts.
	class ProgressHelper : public FutureInterface<void>
	{
	public:

		virtual void cancel() override {
			FutureInterface<void>::cancel();
			// Since this task runs in the main application thread, canceling it
			// means immediate termination.
			reportFinished();
		}
	};

public:

	/// \brief Constructor.
	Q_INVOKABLE PythonScriptModifier(DataSet* dataset);

	/// This modifies the input data.
	virtual PipelineStatus modifyObject(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state) override;

	/// Returns the current status of the modifier.
	virtual PipelineStatus status() const override { return _modifierStatus; }

	/// Sets the status returned by the modifier and generates a ReferenceEvent::ObjectStatusChanged event.
	void setStatus(const PipelineStatus& status);

	/// Returns the Python script.
	const QString& script() const { return _script; }

	/// Sets the Python script.
	void setScript(const QString& script) { _script = script; }

	/// Returns the Python script function executed by the modifier.
	boost::python::object scriptFunction() {
		if(_modifyScriptFunction) return *_modifyScriptFunction;
		else return boost::python::object();
	}

	/// Sets the Python script function to be executed by the modifier.
	void setScriptFunction(const boost::python::object& func) {
		_modifyScriptFunction = func;
		invalidateCachedResults(false);
	}

	/// Returns the log output generated by the script.
	const QString& scriptLogOutput() const { return _scriptLogOutput; }

	/// Interrupts a running script engine if there is one.
	void stopRunningScript();

	/// Asks this object to delete itself. Calls stopRunningEngine() first.
	virtual void deleteReferenceObject() override;

	/// Returns whether the modifier can be applied to the given input data.
	virtual bool isApplicableTo(const PipelineFlowState& input) override { return true; }

protected:

	/// This method is called by the system when the upstream modification pipeline has changed.
	virtual void upstreamPipelineChanged(ModifierApplication* modApp) override;

	/// Is called when the value of a property of this object has changed.
	virtual void propertyChanged(const PropertyFieldDescriptor& field) override;

	/// Invalidates the modifier's result cache so that the results will be recomputed
	/// next time the modifier is evaluated.
	void invalidateCachedResults(bool discardCache);

	/// Executes the Python script function to compute the modifier results.
	Q_INVOKABLE void runScriptFunction();

	/// Compiles the script entered by the user.
	void compileScript();

	/// This is called when the script function was successfully completed.
	void scriptCompleted();

private Q_SLOTS:

	/// Is called when the script generates some output.
	void onScriptOutput(const QString& text);

private:

	/// The Python script.
	PropertyField<QString> _script;

	/// The Python engine.
	std::unique_ptr<ScriptEngine> _scriptEngine;

	/// The compiled modify() function.
	boost::optional<boost::python::object> _modifyScriptFunction;

	/// The log output generated by the script.
	QString _scriptLogOutput;

	/// The cached computation results.
	PipelineFlowState _outputCache;

	/// The cached input data.
	PipelineFlowState _inputCache;

	/// The validity interval of the currently running computing task.
	TimeInterval _computingInterval;

	/// The status returned by the modifier.
	PipelineStatus _modifierStatus;

	/// Flag that indicates that the script is going to be run soon.
	bool _scriptExecutionQueued;

	/// The running script task.
	std::shared_ptr<ProgressHelper> _runningTask;

	/// The generator object returned by the script function.
	boost::optional<boost::python::object> _generatorObject;

	/// The namespace (scope) the script will be executed in.
	boost::optional<boost::python::dict> _mainNamespacePrototype;

	/// The DataCollection passed to the Python script.
	OORef<CompoundObject> _dataCollection;

	Q_CLASSINFO("DisplayName", "Python script");
	Q_CLASSINFO("ModifierCategory", "Modification");

	Q_OBJECT
	OVITO_OBJECT

	DECLARE_PROPERTY_FIELD(_script);
};

}	// End of namespace

#endif // __OVITO_PYTHON_SCRIPT_MODIFIER_H
