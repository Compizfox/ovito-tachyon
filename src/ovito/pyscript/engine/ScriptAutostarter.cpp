///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2017) Alexander Stukowski
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

#include <ovito/pyscript/PyScript.h>
#include <ovito/core/app/Application.h>
#include <ovito/core/app/StandaloneApplication.h>
#include <ovito/core/dataset/DataSetContainer.h>
#include "ScriptAutostarter.h"
#include "ScriptEngine.h"

namespace PyScript {

IMPLEMENT_OVITO_CLASS(ScriptAutostarter);

/******************************************************************************
* Destructor, which is called at program exit.
******************************************************************************/
ScriptAutostarter::~ScriptAutostarter()
{
	// Shut down Python interpreter.
	// This will run the Python functions registered with the 'atexit' module.
	if(Py_IsInitialized()) {
		py::finalize_interpreter();
	}
}

/******************************************************************************
* Registers plugin-specific command line options.
******************************************************************************/
void ScriptAutostarter::registerCommandLineOptions(QCommandLineParser& cmdLineParser)
{
	// Register the --script command line option.
	cmdLineParser.addOption(QCommandLineOption("script", tr("Runs a Python script file."), tr("FILE")));

	// Register the --scriptarg command line option.
	cmdLineParser.addOption(QCommandLineOption("scriptarg", tr("Passes a command line option to the Python script."), tr("ARG")));

	// Register the --exec command line option.
	cmdLineParser.addOption(QCommandLineOption("exec", tr("Executes a single Python statement."), tr("CMD")));
}

/******************************************************************************
* Is called after the application has been completely initialized.
******************************************************************************/
void ScriptAutostarter::applicationStarted()
{
	// Execute the script commands and files passed on the command line.
	QStringList scriptCommands = StandaloneApplication::instance()->cmdLineParser().values("exec");
	QStringList scriptFiles = StandaloneApplication::instance()->cmdLineParser().values("script");

	if((!scriptCommands.empty() || !scriptFiles.empty()) && Application::instance()->datasetContainer()) {

		// Get the current dataset.
		DataSet* dataset = Application::instance()->datasetContainer()->currentSet();

		// Suppress undo recording. Actions performed by startup scripts cannot be undone.
		UndoSuspender noUndo(dataset);

		// Pass command line parameters to the script.
		QStringList scriptArguments = StandaloneApplication::instance()->cmdLineParser().values("scriptarg");

		// This is the async operation object used when calling script functions in the following.
		AsyncOperation scriptOperation(dataset->taskManager());

		// Execute script commands.
		for(int index = scriptCommands.size() - 1; index >= 0; index--) {
			const QString& command = scriptCommands[index];
			try {
				ScriptEngine::executeCommands(command, dataset, scriptOperation.task(), nullptr, true, scriptArguments);
			}
			catch(Exception& ex) {
				ex.prependGeneralMessage(tr("Error during Python script execution."));
				throw;
			}
		}

		// Execute script files.
		for(int index = scriptFiles.size() - 1; index >= 0; index--) {
			ScriptEngine::executeFile(scriptFiles[index], dataset, scriptOperation.task(), nullptr, true, scriptArguments);
		}
	}
}

};