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
 * \brief Contains forward declarations of OVITO's core classes and namespaces.
 */

#pragma once


// Sub-namespaces are only used for the documentation generated by Doxygen,
// because current C++ compilers do not fully support C++11 inline namespaces yet.
#define OVITO_BEGIN_INLINE_NAMESPACE(ns)
#define OVITO_END_INLINE_NAMESPACE

namespace Ovito {

	class Application;
	OVITO_BEGIN_INLINE_NAMESPACE(Util)
		OVITO_BEGIN_INLINE_NAMESPACE(IO)
			class FileManager;
			class ObjectSaveStream;
			class ObjectLoadStream;
			class CompressedTextReader;
			class CompressedTextWriter;
			OVITO_BEGIN_INLINE_NAMESPACE(Internal)
				class VideoEncoder;
				class SftpJob;
			OVITO_END_INLINE_NAMESPACE
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(Concurrency)
			class PromiseBase;
			class FutureBase;
			class Task;
			class TaskManager;
			class TaskWatcher;
			class AsynchronousTaskBase;
			class TrackingTask;
			class MainThreadTask;
			class AsyncOperation;
			template<typename tuple_type> class ContinuationTask;
			template<typename... R> class Future;
			template<typename... R> class SharedFuture;
			template<typename... R> class Promise;
			template<class BaseState, class tuple_type> class TaskWithResultStorage;
			using TaskPtr = std::shared_ptr<Task>;
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(Mesh)
			class TriMesh;
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(Math)
		OVITO_END_INLINE_NAMESPACE
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(Anim)
		class Controller;
		class AnimationSettings;
		class LookAtController;
		class KeyframeController;
		class PRSTransformationController;
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(PluginSystem)
		class Plugin;
		class PluginManager;
		class ApplicationService;
		OVITO_BEGIN_INLINE_NAMESPACE(Internal)
			class NativePlugin;
		OVITO_END_INLINE_NAMESPACE
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(ObjectSystem)
		class OvitoObject;
		class OvitoClass;
		using OvitoClassPtr = const OvitoClass*;
		class CloneHelper;
		class RefMaker;
		class RefMakerClass;
		class RefTarget;
		class PropertyFieldDescriptor;
		class PropertyFieldBase;
		template<typename property_data_type> class RuntimePropertyField;
		template<typename property_data_type> class PropertyField;
		class SingleReferenceFieldBase;
		template<typename RefTargetType> class ReferenceField;
		class VectorReferenceFieldBase;
		template<typename RefTargetType> class VectorReferenceField;
		class DataSet;
		class DataSetContainer;
		OVITO_BEGIN_INLINE_NAMESPACE(Units)
			class ParameterUnit;
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(Undo)
			class UndoStack;
			class UndoableOperation;
		OVITO_END_INLINE_NAMESPACE
		OVITO_BEGIN_INLINE_NAMESPACE(Scene)
			class SceneNode;
			class DataObject;
			class DataObjectReference;
			class AttributeDataObject;
			class RootSceneNode;
			class SelectionSet;
			class Modifier;
			class ModifierClass;
			using ModifierClassPtr = const ModifierClass*;
			class ModifierApplication;
			class PipelineSceneNode;
			class PipelineFlowState;
			class DataCollection;
			class PipelineObject;
			class PipelineCache;
			class CachingPipelineObject;
			class DataVis;
			class TransformingDataVis;
			class StaticSource;
			class ModifierDelegate;
			class DelegatingModifier;
			class MultiDelegatingModifier;
			class AsynchronousModifier;
			class AsynchronousModifierDelegate;
			class AsynchronousDelegatingModifier;
			OVITO_BEGIN_INLINE_NAMESPACE(StdObj)
				class AbstractCameraObject;
			OVITO_END_INLINE_NAMESPACE
		OVITO_END_INLINE_NAMESPACE
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(Rendering)
		class SceneRenderer;
		class ObjectPickInfo;
		class ViewportPickResult;
		class RenderSettings;
		class FrameBuffer;
		OVITO_BEGIN_INLINE_NAMESPACE(Internal)
		OVITO_END_INLINE_NAMESPACE
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(View)
		class Viewport;
		class ViewportConfiguration;
		class ViewportSettings;
		struct ViewProjectionParameters;
		class ViewportOverlay;
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(DataIO)
		class FileImporter;
		class FileImporterClass;
		class FileExporter;
		class FileExporterClass;
		class FileSource;
		class FileSourceImporter;
	OVITO_END_INLINE_NAMESPACE
	OVITO_BEGIN_INLINE_NAMESPACE(Gui)
		class MainWindow;
	OVITO_END_INLINE_NAMESPACE

	// This should only be visible to Doxygen:
#ifdef ONLY_FOR_DOXYGEN
	using namespace Util;
	using namespace Util::IO;
	using namespace Util::Math;
	using namespace Util::Mesh;
	using namespace Util::Concurrency;
	using namespace Rendering;
	using namespace View;
	using namespace DataIO;
	using namespace Anim;
	using namespace PluginSystem;
	using namespace ObjectSystem;
	using namespace ObjectSystem::Units;
	using namespace ObjectSystem::Undo;
	using namespace ObjectSystem::Scene;
	using namespace ObjectSystem::Scene::StdObj;
#endif
}


