////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2014 Alexander Stukowski
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

/**
 * \file
 * \brief Contains forward declarations of OVITO's GUI classes and namespaces.
 */

#pragma once


// Sub-namespaces are only used for the documentation generated by Doxygen,
// because current C++ compilers do not fully support C++11 inline namespaces yet.
#define OVITO_BEGIN_INLINE_NAMESPACE(ns)
#define OVITO_END_INLINE_NAMESPACE

namespace Ovito {

	OVITO_BEGIN_INLINE_NAMESPACE(PluginSystem)
		class UtilityApplet;
		class GuiAutoStartObject;
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(Rendering)
		class OpenGLSceneRenderer;
		class ViewportSceneRenderer;
		class StandardSceneRenderer;
		OVITO_BEGIN_INLINE_NAMESPACE(Internal)
		OVITO_END_INLINE_NAMESPACE
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(View)
		OVITO_BEGIN_INLINE_NAMESPACE(Internal)
			class PickingSceneRenderer;
		OVITO_END_INLINE_NAMESPACE
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(Gui)
		class MainWindow;
		class GuiApplication;
		class ActionManager;
		class DataInspectionApplet;
		OVITO_BEGIN_INLINE_NAMESPACE(Widgets)
			class PropertiesPanel;
			class SpinnerWidget;
			class ColorPickerWidget;
			class RolloutContainer;
			class FrameBufferWindow;
			class FrameBufferWidget;
			class AutocompleteTextEdit;
			class AutocompleteLineEdit;
			class ElidedTextLabel;
			class HtmlListWidget;
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(Params)
			class PropertiesEditor;
			class AffineTransformationParameterUI;
			class BooleanParameterUI;
			class BooleanActionParameterUI;
			class BooleanGroupBoxParameterUI;
			class BooleanRadioButtonParameterUI;
			class ColorParameterUI;
			class CustomParameterUI;
			class FilenameParameterUI;
			class FloatParameterUI;
			class FontParameterUI;
			class IntegerRadioButtonParameterUI;
			class IntegerParameterUI;
			class ParameterUI;
			class RefTargetListParameterUI;
			class StringParameterUI;
			class SubObjectParameterUI;
			class VariantComboBoxParameterUI;
			class Vector3ParameterUI;
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(ViewportInput)
			class ViewportInputManager;
			class ViewportInputMode;
			class ViewportModeAction;
			class ViewportGizmo;
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(Dialogs)
			class FileExporterSettingsDialog;
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(Internal)
			class CoordinateDisplayWidget;
			class CommandPanel;
			class DataInspectorPanel;
			class ModifyCommandPage;
			class RenderCommandPage;
			class OverlayCommandPage;
			class UtilityCommandPage;
			class ViewportMenu;
			class ViewportWindow;
		OVITO_END_INLINE_NAMESPACE
	OVITO_END_INLINE_NAMESPACE

	// This should only be visible to Doxygen:
#ifdef ONLY_FOR_DOXYGEN
	using namespace Gui;
	using namespace Gui::Params;
	using namespace Gui::Dialogs;
	using namespace Gui::ViewportInput;
	using namespace Gui::Widgets;
#endif
}


