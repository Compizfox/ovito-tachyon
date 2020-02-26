////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2019 Alexander Stukowski
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

#include <ovito/gui/desktop/GUI.h>
#include <ovito/gui/desktop/app/GuiApplication.h>

/**
 * This is the main entry point for the graphical desktop application.
 *
 * Note that most of the application logic is found in the Core and the Gui
 * modules of OVITO.
 */
int main(int argc, char** argv)
{
	// Initialize the application.
	Ovito::GuiApplication app;
	if(!app.initialize(argc, argv))
		return 1;

	// Enter event loop.
	int result = app.runApplication();

	// Shut application down.
	app.shutdown();

	return result;
}
