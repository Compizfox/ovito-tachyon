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


#include <ovito/stdobj/StdObj.h>
#include <ovito/stdobj/properties/PropertyStorage.h>
#include <ovito/stdobj/properties/PropertyAccess.h>

namespace Ovito { namespace StdObj {

/**
 * Helper structure storing the data for a PropertyContainer loaded by a file importer.
 */
class OVITO_STDOBJ_EXPORT PropertyContainerImportData
{
public:

	/// Describes a single element type.
	struct TypeDefinition {
		int id;
		QString name;
		std::string name8bit;
		bool preserveNumericId = false;
		QVariantMap attributes;
	};

	/// Used to store the lists of element types.
	class OVITO_STDOBJ_EXPORT TypeList
	{
	public:

		/// Constructor, which takes the class of type elements defined by this list.
		explicit TypeList(const OvitoClass& elementClass) : _elementClass(elementClass) {}

		/// Returns the class of element types defined in this list.
		const OvitoClass& elementClass() const { return _elementClass; }

		/// Defines a new element type with the given numeric id and no name.
		void addTypeId(int id) {
			if(hasTypeId(id)) return;
			_types.push_back({id});
		}

		/// Defines a new type with the given numeric id and name.
		void addNamedTypeId(int id, const QString& name, bool preserveNumericId, const Color& color = Color(0,0,0), FloatType radius = 0, FloatType mass = 0) {
			if(hasTypeId(id)) return;
			_types.push_back({ id, name, name.toStdString(), preserveNumericId });
			if(color != Color(0,0,0)) _types.back().attributes.insert(QStringLiteral("color"), QVariant::fromValue(color));
			if(radius) _types.back().attributes.insert(QStringLiteral("radius"), QVariant::fromValue(radius));
			if(mass) _types.back().attributes.insert(QStringLiteral("mass"), QVariant::fromValue(mass));
		}

		/// Checks if a type with the given numeric id already exists.
		bool hasTypeId(int id) const {
			for(const auto& type : _types) {
				if(type.id == id)
					return true;
			}
			return false;
		}

		/// Changes the name of an existing type.
		void setTypeName(int id, const QString& name, bool preserveNumericId) {
			for(auto& type : _types) {
				if(type.id == id) {
					type.name = name;
					type.name8bit = name.toStdString();
					type.preserveNumericId = preserveNumericId;
					break;
				}
			}
		}

		/// Gives an existing type a name unless it already has one.
		void setTypeName(int id, const char* name, size_t length, bool preserveNumericId) {
			for(auto& type : _types) {
				if(type.id == id) {
					if(type.name8bit.compare(0, length, name, length)) {
						type.name8bit.assign(name, length);
						type.name = QString::fromUtf8(name, length);
						type.preserveNumericId = preserveNumericId;
					}
					break;
				}
			}
		}

		/// Assigns an attribute to an existing type.
		void setTypeAttribute(int id, const QString& attrName, QVariant value) {
			for(auto& type : _types) {
				if(type.id == id) {
					type.attributes.insert(attrName, std::move(value));
					break;
				}
			}
		}

		/// Changes the mass of an existing particle type.
		void setTypeMass(int id, FloatType mass) {
			setTypeAttribute(id, QStringLiteral("mass"), mass);
		}

		/// Changes the radius of an existing type.
		void setTypeRadius(int id, FloatType radius) {
			setTypeAttribute(id, QStringLiteral("radius"), radius);
		}

		/// Defines a new type with the given name.
		inline int addTypeName(const QString& name) {
			for(const auto& type : _types) {
				if(type.name == name)
					return type.id;
			}
			int id = _types.size() + 1;
			_types.push_back({ id, name, name.toStdString() });
			return id;
		}

		/// Defines a new type with the given name, color, and radius.
		int addTypeName(const char* name, const char* name_end, const Color& color = Color(0,0,0), FloatType radius = 0, FloatType mass = 0) {
			size_t nameLen = (name_end ? (name_end - name) : qstrlen(name));
			for(const auto& type : _types) {
				if(type.name8bit.compare(0, type.name8bit.size(), name, nameLen) == 0)
					return type.id;
			}
			int id = _types.size() + 1;
			_types.push_back({ id, QString::fromLocal8Bit(name, nameLen), std::string(name, nameLen) });
			if(color != Color(0,0,0)) _types.back().attributes.insert(QStringLiteral("color"), QVariant::fromValue(color));
			if(radius) _types.back().attributes.insert(QStringLiteral("radius"), QVariant::fromValue(radius));
			if(mass) _types.back().attributes.insert(QStringLiteral("mass"), QVariant::fromValue(mass));
			return id;
		}

		/// Returns the list of element types.
		const std::vector<TypeDefinition>& types() const { return _types; }

		/// Returns the list of element types.
		std::vector<TypeDefinition>& types() { return _types; }

		/// Sorts the element types w.r.t. their name. Reassigns the per-element type IDs.
		/// This method is used by file parsers that create element types on the go while the read the data.
		/// In such a case, the assignment of IDs to types depends on the storage order of data element in the file, which is not desirable.
		void sortTypesByName(PropertyAccess<int>& typeProperty);

		/// Sorts element types according to numeric identifier.
		void sortTypesById();

	private:

		/// The list of defined element types.
		std::vector<TypeDefinition> _types;

		/// The kind of type elements defined in this list (particles types, bond types, etc.).
		const OvitoClass& _elementClass; 
	};

public:

	/// Returns the list of properties.
	const std::vector<PropertyPtr>& properties() const { return _properties; }

	/// Returns a standard property if already defined.
	PropertyPtr findStandardProperty(int which) const {
		OVITO_ASSERT(which != PropertyStorage::GenericUserProperty);
		for(const auto& prop : _properties)
			if(prop->type() == which)
				return prop;
		return {};
	}

	/// Finds a property by name.
	PropertyPtr findProperty(const QString& name) const {
		for(const auto& prop : _properties)
			if(prop->name() == name)
				return prop;
		return {};
	}

	/// Adds a new property.
	PropertyStorage* addProperty(PropertyPtr property) {
		_properties.push_back(std::move(property));
		return _properties.back().get();
	}

	/// Create a standard property.
	template<class ContainerClass>
	PropertyStorage* createStandardProperty(size_t elementCount, int typeId, bool initializeMemory) {
		return addProperty(ContainerClass::OOClass().createStandardStorage(elementCount, typeId, initializeMemory));
	}

	/// Removes a property from the list.
	void removeProperty(int index) {
		OVITO_ASSERT(index >= 0 && index < _properties.size());
		_typeLists.erase(_properties[index].get());
		_properties.erase(_properties.begin() + index);
	}

	/// Removes a property from the list.
	void removeProperty(const PropertyPtr& property) {
		auto iter = std::find(_properties.begin(), _properties.end(), property);
		OVITO_ASSERT(iter != _properties.end());
		_typeLists.erase(property.get());
		_properties.erase(iter);
	}

	/// Returns the list of types defined for a property.
	TypeList* propertyTypesList(const PropertyStorage* property) {
		auto typeList = _typeLists.find(property);
		if(typeList == _typeLists.end())
			return nullptr;
		return typeList->second.get();
	}

	/// Returns the list of types defined for a property.
	TypeList* propertyTypesList(const PropertyPtr& property) {
		return propertyTypesList(property.get());
	}

	/// Returns the list of types defined for a property.
	TypeList* propertyTypesList(const PropertyAccess<int>& property) {
		return propertyTypesList(property.storage());
	}

	/// Creates a types list for a property.
	TypeList* createPropertyTypesList(const PropertyStorage* property, const OvitoClass& elementClass) {
		auto typeList = _typeLists.find(property);
		if(typeList == _typeLists.end())
			typeList = _typeLists.emplace(property, std::make_unique<TypeList>(elementClass)).first;
		return typeList->second.get();
	}

	/// Creates a types list for a property.
	TypeList* createPropertyTypesList(const PropertyPtr& property, const OvitoClass& elementClass) {
		return createPropertyTypesList(property.get(), elementClass);
	}

	/// Creates a types list for a property.
	TypeList* createPropertyTypesList(const PropertyAccess<int>& property, const OvitoClass& elementClass) {
		return createPropertyTypesList(property.storage(), elementClass);
	}	

	/// Sets the list of types defined for a property.
	void setPropertyTypesList(const PropertyStorage* property, std::unique_ptr<TypeList> list) {
		_typeLists.emplace(property, std::move(list));
	}

	/// Sets the list of types defined for a property.
	void setPropertyTypesList(const PropertyPtr& property, std::unique_ptr<TypeList> list) {
		setPropertyTypesList(property.get(), std::move(list));
	}

	/// Sets the list of types defined for a property.
	void setPropertyTypesList(const PropertyAccess<int>& property, std::unique_ptr<TypeList> list) {
		setPropertyTypesList(property.storage(), std::move(list));
	}

	/// Sorts the data elements with respect to their unique IDs.
	/// Does nothing if data element do not have IDs.
	std::vector<size_t> sortElementsById();

	/// Selects the type of visual element to create for rendering the data elements.
	void setVisElementClass(OvitoClassPtr visElementClass) { _visElementClass = visElementClass; }

	/// Transfers the internal property data to the target PropertyContainer.
	void transferToContainer(const PropertyContainer* existingContainer, PropertyContainer* targetContainer, bool isNewFile, CloneHelper& cloneHelper) const;

private:

	/// Inserts the stored element types into the given property object.
	void insertTypes(PropertyObject* propertyObj, TypeList* typeList, bool isNewFile) const;

private:

	/// List of properties.
	std::vector<PropertyPtr> _properties;

#if defined(Q_CC_MSVC)
	/// Needed as a workaround for a compiler bug in MVSC 2017.
	/// See https://stackoverflow.com/questions/21056872/c-stdunique-ptr-wont-compile-in-map
	/// Problem might be fixed in MSVC 2019.
	std::unique_ptr<TypeList> _dummyWorkaroundForMVSCCompilerBug;
#endif

	/// Stores the lists of types for typed properties (particle/bond/angle/dihedral/improper properties).
	std::map<const PropertyStorage*, std::unique_ptr<TypeList>> _typeLists;

	/// The type of visual element to create for rendering the data elements.
	OvitoClassPtr _visElementClass = {};
};

}	// End of namespace
}	// End of namespace
