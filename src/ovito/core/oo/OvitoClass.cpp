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

#include <ovito/core/Core.h>
#include <ovito/core/app/PluginManager.h>
#include <ovito/core/utilities/io/ObjectSaveStream.h>
#include <ovito/core/utilities/io/ObjectLoadStream.h>
#include <ovito/core/oo/OvitoObject.h>
#include <ovito/core/oo/RefTarget.h>
#include <ovito/core/dataset/DataSet.h>
#include "OvitoClass.h"

namespace Ovito {

// Head of linked list of native meta-classes.
OvitoClass* OvitoClass::_firstMetaClass = nullptr;

/******************************************************************************
* Constructor.
******************************************************************************/
OvitoClass::OvitoClass(const QString& name, OvitoClassPtr superClass, const char* pluginId, const QMetaObject* qtClassInfo) :
	_name(name),
	_displayName(name),
	_plugin(nullptr),
	_superClass(superClass),
	_isAbstract(false),
	_pluginId(pluginId),
	_qtClassInfo(qtClassInfo)
{
	OVITO_ASSERT(superClass != nullptr || name == QStringLiteral("OvitoObject"));

	// Insert into linked list of all object types.
	_nextMetaclass = _firstMetaClass;
	_firstMetaClass = this;
}

/******************************************************************************
* Is called by the system after construction of the meta-class instance.
******************************************************************************/
void OvitoClass::initialize()
{
	// Remove namespace qualifier from Qt's class name.
	if(qtMetaObject()) {
		// Mark classes as abstract that don't have an invokable constructor.
		setAbstract(qtMetaObject()->constructorCount() == 0);

		_pureClassName = qtMetaObject()->className();
		for(const char* p = _pureClassName; *p != '\0'; p++) {
			if(p[0] == ':' && p[1] == ':') {
				p++;
				_pureClassName = p+1;
			}
		}

		// Interpret Qt class info fields.
		for(int i = qtMetaObject()->classInfoOffset(); i < qtMetaObject()->classInfoCount(); i++) {
			if(qstrcmp(qtMetaObject()->classInfo(i).name(), "DisplayName") == 0) {
				// Fetch display name assigned to the Qt object class.
				setDisplayName(QString::fromLocal8Bit(qtMetaObject()->classInfo(i).value()));
			}
			else if(qstrcmp(qtMetaObject()->classInfo(i).name(), "ClassNameAlias") == 0) {
				// Load name alias assigned to the Qt object class.
				setNameAlias(QString::fromLocal8Bit(qtMetaObject()->classInfo(i).value()));
			}
		}
	}
	else {
		setAbstract(true);
	}
}

/******************************************************************************
* Returns a human-readable string describing this class.
******************************************************************************/
QString OvitoClass::descriptionString() const
{
	for(int i = qtMetaObject()->classInfoOffset(); i < qtMetaObject()->classInfoCount(); i++) {
		if(qstrcmp(qtMetaObject()->classInfo(i).name(), "Description") == 0) {
			return QString::fromUtf8(qtMetaObject()->classInfo(i).value());
		}
	}
	return QString();
}

/******************************************************************************
* Determines if an object is an instance of the class or one of its subclasses.
******************************************************************************/
bool OvitoClass::isMember(const OvitoObject* obj) const
{
	if(!obj) return false;
	return obj->getOOClass().isDerivedFrom(*this);
}

/******************************************************************************
* Creates an instance of this object class.
* Throws an exception if the containing plugin failed to load.
******************************************************************************/
OORef<OvitoObject> OvitoClass::createInstance(DataSet* dataset) const
{
	if(plugin()) {
		OVITO_CHECK_POINTER(plugin());
		if(!plugin()->isLoaded()) {
			// Load plugin first.
			try {
				plugin()->loadPlugin();
			}
			catch(Exception& ex) {
				ex.prependGeneralMessage(OvitoObject::tr("Could not create instance of class %1. Failed to load plugin '%2'").arg(name()).arg(plugin()->pluginId()));
				throw ex;
			}
		}
	}
	if(isAbstract())
		throw Exception(OvitoObject::tr("Cannot instantiate abstract class '%1'.").arg(name()), dataset);

	OVITO_ASSERT_MSG(!isDerivedFrom(RefTarget::OOClass()) || dataset != nullptr || *this == DataSet::OOClass(), "OvitoClass::createInstance()", "Tried to create instance of RefTarget derived class without passing a DatSet.");
	OVITO_ASSERT_MSG(isDerivedFrom(RefTarget::OOClass()) || dataset == nullptr, "OvitoClass::createInstance()", "Passed a DatSet to the constructor of a class that is not derived from RefTarget.");

	return createInstanceImpl(dataset);
}

/******************************************************************************
* Creates an instance of this object class.
******************************************************************************/
OvitoObject* OvitoClass::createInstanceImpl(DataSet* dataset) const
{
#ifdef OVITO_DEBUG
	// Check if class hierarchy is consistent.
	OvitoClassPtr ovitoSuperClass = superClass();
	while(ovitoSuperClass && ovitoSuperClass->qtMetaObject() == nullptr)
		ovitoSuperClass = ovitoSuperClass->superClass();
	OVITO_ASSERT(ovitoSuperClass != nullptr);
	const QMetaObject* qtSuperClass = qtMetaObject()->superClass();
	while(qtSuperClass && qtSuperClass != ovitoSuperClass->qtMetaObject())
		qtSuperClass = qtSuperClass->superClass();
	OVITO_ASSERT_MSG(qtSuperClass != nullptr, "OvitoClass::createInstanceImpl", qPrintable(QString("Class %1 is not derived from base class %2 as specified by the object type descriptor. Qt super class is %3.").arg(name()).arg(superClass()->name()).arg(qtMetaObject()->superClass()->className())));
#endif

	OvitoObject* obj;

	if(isDerivedFrom(RefTarget::OOClass()) && *this != DataSet::OOClass()) {
		OVITO_ASSERT(dataset != nullptr);
		UndoSuspender noUndo(dataset->undoStack());
		obj = qobject_cast<OvitoObject*>(qtMetaObject()->newInstance(Q_ARG(DataSet*, dataset)));
	}
	else {
		obj = qobject_cast<OvitoObject*>(qtMetaObject()->newInstance());
	}

	if(!obj)
		throw Exception(OvitoObject::tr("Failed to instantiate class '%1'.").arg(name()), dataset);

	return obj;
}

/******************************************************************************
* Writes a class descriptor to the stream. This is for internal use of the core only.
******************************************************************************/
void OvitoClass::serializeRTTI(SaveStream& stream, OvitoClassPtr type)
{
	stream.beginChunk(0x10000000);
	if(type) {
		stream << type->plugin()->pluginId();
		stream << type->name();
	}
	else {
		stream << QString() << QString();
	}
	stream.endChunk();
}

/******************************************************************************
* Loads a class descriptor from the stream. This is for internal use of the core only.
* Throws an exception if the class is not defined or the required plugin is not installed.
******************************************************************************/
OvitoClassPtr OvitoClass::deserializeRTTI(LoadStream& stream)
{
	QString pluginId, className;
	stream.expectChunk(0x10000000);
	stream >> pluginId;
	stream >> className;
	stream.closeChunk();

	if(pluginId.isEmpty() && className.isEmpty())
		return nullptr;

	try {
		// Look plugin.
		Plugin* plugin = PluginManager::instance().plugin(pluginId);
		if(!plugin) {

			// If plugin does not exist anymore, fall back to searching other plugins for the requested class.
			for(Plugin* otherPlugin : PluginManager::instance().plugins()) {
				if(OvitoClassPtr clazz = otherPlugin->findClass(className))
					return clazz;
			}

			throw Exception(OvitoObject::tr("A required plugin is not installed: %1").arg(pluginId));
		}
		OVITO_CHECK_POINTER(plugin);

		// Look up class descriptor.
		OvitoClassPtr clazz = plugin->findClass(className);
		if(!clazz) {

			// If class does not exist in the plugin anymore, fall back to searching other plugins for the requested class.
			for(Plugin* otherPlugin : PluginManager::instance().plugins()) {
				if(OvitoClassPtr clazz = otherPlugin->findClass(className))
					return clazz;
			}

			throw Exception(OvitoObject::tr("Required class '%1' not found in plugin '%2'.").arg(className, pluginId));
		}

		return clazz;
	}
	catch(Exception& ex) {
		ex.prependGeneralMessage(OvitoObject::tr("File cannot be loaded, because it contains object types that are not (or no longer) available in this program version."));
		throw ex;
	}
}

/******************************************************************************
* Encodes the plugin ID and the class name in a string.
******************************************************************************/
QString OvitoClass::encodeAsString(OvitoClassPtr type)
{
	OVITO_CHECK_POINTER(type);
	return type->plugin()->pluginId() + QStringLiteral("::") + type->name();
}

/******************************************************************************
* Decodes a class descriptor from a string, which has been generated by encodeAsString().
******************************************************************************/
OvitoClassPtr OvitoClass::decodeFromString(const QString& str)
{
	QStringList tokens = str.split(QStringLiteral("::"));
	if(tokens.size() != 2)
		throw Exception(OvitoObject::tr("Invalid type or encoding: %1").arg(str));

	// Look up plugin.
	Plugin* plugin = PluginManager::instance().plugin(tokens[0]);
	if(!plugin) {

		// If plugin does not exist anymore, fall back to searching other plugins for the requested class.
		for(Plugin* otherPlugin : PluginManager::instance().plugins()) {
			if(OvitoClassPtr clazz = otherPlugin->findClass(tokens[1]))
				return clazz;
		}

		throw Exception(OvitoObject::tr("A required plugin is not installed: %1").arg(tokens[0]));
	}
	OVITO_CHECK_POINTER(plugin);

	// Look up class descriptor.
	OvitoClassPtr clazz = plugin->findClass(tokens[1]);
	if(!clazz) {

		// If class does not exist in the plugin anymore, fall back to searching other plugins for the requested class.
		for(Plugin* otherPlugin : PluginManager::instance().plugins()) {
			if(OvitoClassPtr clazz = otherPlugin->findClass(tokens[1]))
				return clazz;
		}

		throw Exception(OvitoObject::tr("Required class '%1' not found in plugin '%2'.").arg(tokens[1], tokens[0]));
	}

	return clazz;
}

}	// End of namespace
