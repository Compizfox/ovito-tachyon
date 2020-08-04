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

#include <ovito/stdobj/StdObj.h>
#include <ovito/stdobj/properties/PropertyContainer.h>
#include <ovito/stdobj/properties/ElementType.h>
#include <ovito/core/app/Application.h>
#include <ovito/core/dataset/DataSet.h>
#include "PropertyContainerImportData.h"

#include <boost/container/flat_map.hpp>

namespace Ovito { namespace StdObj {

/******************************************************************************
* Sorts the types w.r.t. their name. Reassigns the per-element type IDs too.
* This method is used by file parsers that create element types on the
* go while the read the data. In such a case, the assignment of IDs to types
* depends on the storage order of data elements in the file, which is not desirable.
******************************************************************************/
void PropertyContainerImportData::TypeList::sortTypesByName(PropertyAccess<int>& typeProperty)
{
	// Check if type IDs form a consecutive sequence starting at 1.
	// If not, we leave the type order as it is.
	for(size_t index = 0; index < _types.size(); index++) {
		if(_types[index].id != index + 1)
			return;
	}

	// Check if types are already in the correct order.
	auto compare = [](const TypeDefinition& a, const TypeDefinition& b) { return a.name.compare(b.name) < 0; };
	if(std::is_sorted(_types.begin(), _types.end(), compare))
		return;

	// Reorder types by name.
	std::sort(_types.begin(), _types.end(), compare);

	// Build map of IDs.
	std::vector<int> mapping(_types.size() + 1);
	for(size_t index = 0; index < _types.size(); index++) {
		mapping[_types[index].id] = index + 1;
		_types[index].id = index + 1;
	}

	// Remap type IDs.
	if(typeProperty) {
		for(int& t : typeProperty) {
			OVITO_ASSERT(t >= 1 && t < mapping.size());
			t = mapping[t];
		}
	}
}

/******************************************************************************
* Sorts particle/bond types according numeric identifier.
******************************************************************************/
void PropertyContainerImportData::TypeList::sortTypesById()
{
	auto compare = [](const TypeDefinition& a, const TypeDefinition& b) { return a.id < b.id; };
	std::sort(_types.begin(), _types.end(), compare);
}

/******************************************************************************
* Transfers the internal property data to the target PropertyContainer.
******************************************************************************/
void PropertyContainerImportData::transferToContainer(const PropertyContainer* existingContainer, PropertyContainer* targetContainer, bool isNewFile, CloneHelper& cloneHelper) const
{
	if(!existingContainer) {
		// Create the vis element requested by the file importer.
		if(_visElementClass && (!targetContainer->visElement() || _visElementClass != &targetContainer->visElement()->getOOMetaClass()))
			targetContainer->setVisElement(static_object_cast<DataVis>(_visElementClass->createInstance(targetContainer->dataset())));
		else if(!_visElementClass)
			targetContainer->setVisElement(nullptr);

		// Initialize the property container and its vis element to default values.
		if(Application::instance()->executionContext() == Application::ExecutionContext::Interactive)
			targetContainer->loadUserDefaults();
	}
	else {
		// Adopt the existing vis element, or create the right vis element requested by the file importer.
		if(_visElementClass && (!existingContainer->visElement() || _visElementClass != &existingContainer->visElement()->getOOMetaClass())) {
			targetContainer->setVisElement(static_object_cast<DataVis>(_visElementClass->createInstance(targetContainer->dataset())));
			if(Application::instance()->executionContext() == Application::ExecutionContext::Interactive)
				targetContainer->visElement()->loadUserDefaults();
		}
		else if(!_visElementClass)
			targetContainer->setVisElement(nullptr);
		else
			targetContainer->setVisElement(existingContainer->visElement());
	}

	// Transfer properties.
	for(auto& property : _properties) {

		// Look up existing property object.
		const PropertyObject* existingPropertyObj = existingContainer ? 
			((property->type() != PropertyStorage::GenericUserProperty) ? existingContainer->getProperty(property->type()) : existingContainer->getProperty(property->name())) 
			: nullptr;

		OORef<PropertyObject> propertyObj;
		if(existingPropertyObj) {
			propertyObj = cloneHelper.cloneObject(existingPropertyObj, false);
			propertyObj->setStorage(std::move(property));
			targetContainer->addProperty(propertyObj);
		}
		else {
			propertyObj = targetContainer->createProperty(std::move(property));
		}

		// Transfer element types.
		auto typeList = _typeLists.find(propertyObj->storage().get());
		insertTypes(propertyObj, (typeList != _typeLists.end()) ? typeList->second.get() : nullptr, isNewFile);
	}

	targetContainer->verifyIntegrity();
}

/******************************************************************************
* Inserts the element types into the given destination property object.
******************************************************************************/
void PropertyContainerImportData::insertTypes(PropertyObject* typeProperty, TypeList* typeList, bool isNewFile) const
{
	QSet<ElementType*> activeTypes;
	boost::container::flat_map<int, int> typeRemapping;

	if(typeList) {
		// Add the new element types one by one to the property object.
		for(auto& item : typeList->types()) {
			// Look up existing element type.
			OORef<ElementType> elementType;
			if(item.name.isEmpty()) {
				elementType = typeProperty->elementType(item.id);
			}
			else {
				elementType = typeProperty->elementType(item.name);
				if(elementType) {
					if(item.id != elementType->numericId()) {
						typeRemapping.emplace(item.id, elementType->numericId());
					}
				}
				else {
					elementType = typeProperty->elementType(item.id);
					if(elementType && elementType->name() != item.name) {
						elementType = nullptr;
						if(!isNewFile) {
							int mappedId = typeProperty->generateUniqueElementTypeId(item.id + typeList->types().size());
							typeRemapping.emplace(item.id, mappedId);
							item.id = mappedId;
						}
					}
				}
			}
			// Create element type (unless it already exists).
			if(!elementType) {
				elementType = static_object_cast<ElementType>(typeList->elementClass().createInstance(typeProperty->dataset()));
				elementType->setNumericId(item.id);
				elementType->setName(item.name);
				typeProperty->addElementType(elementType);
				elementType->initialize(true, item.attributes, typeProperty->type());
			}
			else {
				// Update attributes of existing element type. 
				// Note that this requires a mutable copy of the ElementType instance,
				// which will be created if the first attempt to update the existing instance fails.
				if(!elementType->initialize(false, item.attributes, typeProperty->type())) {
					// Create a mutable copy of the original ElementType.
					elementType = typeProperty->makeMutable<ElementType>(elementType);
					elementType->initialize(false, item.attributes, typeProperty->type());
				}
			}
			activeTypes.insert(elementType);
		}
	}

	if(isNewFile) {
		// Remove existing element types from the property that are no longer needed.
		for(int index = typeProperty->elementTypes().size() - 1; index >= 0; index--) {
			if(!activeTypes.contains(typeProperty->elementTypes()[index]))
				typeProperty->removeElementType(index);
		}
	}

	// Remap type IDs.
	if(!typeRemapping.empty()) {
		for(int& t : PropertyAccess<int>(typeProperty)) {
			auto iter = typeRemapping.find(t);
			if(iter != typeRemapping.end())
				t = iter->second;
		}
	}
}

/******************************************************************************
* Sorts the data elements with respect to their unique IDs.
* Does nothing if data elements do not have IDs.
******************************************************************************/
std::vector<size_t> PropertyContainerImportData::sortElementsById()
{
	ConstPropertyAccess<qlonglong> ids = findStandardProperty(PropertyStorage::GenericIdentifierProperty);
	if(!ids) return {};

	// Determine new permutation of data elements which sorts them by ascending ID.
	std::vector<size_t> permutation(ids.size());
	std::iota(permutation.begin(), permutation.end(), (size_t)0);
	std::sort(permutation.begin(), permutation.end(), [&](size_t a, size_t b) { return ids[a] < ids[b]; });
	std::vector<size_t> invertedPermutation(ids.size());
	bool isAlreadySorted = true;
	for(size_t i = 0; i < permutation.size(); i++) {
		invertedPermutation[permutation[i]] = i;
		if(permutation[i] != i) isAlreadySorted = false;
	}
	if(isAlreadySorted) return {};

	// Re-order all values in the property arrays.
	for(const PropertyPtr& prop : properties()) {
		PropertyStorage copy(*prop);
		prop->mappedCopyFrom(copy, invertedPermutation);
	}

	return invertedPermutation;
}

}	// End of namespace
}	// End of namespace
