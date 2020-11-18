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
#include <ovito/core/dataset/pipeline/PipelineFlowState.h>
#include "PropertyObject.h"
#include "PropertyReference.h"
#include "PropertyContainerClass.h"

namespace Ovito { namespace StdObj {

/******************************************************************************
* This helper method returns a standard property (if present) from the
* given pipeline state.
******************************************************************************/
void PropertyContainerClass::registerStandardProperty(int typeId, QString name, int dataType, QStringList componentNames, OvitoClassPtr typedPropertyElementClass, QString title)
{
	OVITO_ASSERT_MSG(typeId > 0, "PropertyContainerClass::registerStandardProperty", "Invalid standard property type ID");
	OVITO_ASSERT_MSG(_standardPropertyIds.find(name) == _standardPropertyIds.end(), "PropertyContainerClass::registerStandardProperty", "Duplicate standard property name");
	OVITO_ASSERT_MSG(_standardPropertyNames.find(typeId) == _standardPropertyNames.end(), "PropertyContainerClass::registerStandardProperty", "Duplicate standard property type ID");
	OVITO_ASSERT_MSG(dataType == PropertyStorage::Int || dataType == PropertyStorage::Int64 || dataType == PropertyStorage::Float, "PropertyContainerClass::registerStandardProperty", "Invalid standard property data type");
	OVITO_ASSERT_MSG(!typedPropertyElementClass || typedPropertyElementClass->isDerivedFrom(ElementType::OOClass()), "PropertyContainerClass::registerStandardProperty", "Element type class is not derived from ElementType base");

	_standardPropertyList.push_back(typeId);
	if(!name.isEmpty())
		_standardPropertyIds.insert(name, typeId);
	_standardPropertyNames.emplace(typeId, std::move(name));
	_standardPropertyTitles.emplace(typeId, std::move(title));
	_standardPropertyComponents.emplace(typeId, std::move(componentNames));
	_standardPropertyDataTypes.emplace(typeId, dataType);
	if(typedPropertyElementClass)
		_standardPropertyElementTypes.emplace(typeId, typedPropertyElementClass);
}

/******************************************************************************
* Factory function that creates a property object based on an existing storage.
******************************************************************************/
OORef<PropertyObject> PropertyContainerClass::createFromStorage(DataSet* dataset, PropertyPtr storage) const
{
	OORef<PropertyObject> property = new PropertyObject(dataset, std::move(storage));
	if(property->type() != 0)
		property->setTitle(standardPropertyTitle(property->type()));
	prepareNewProperty(property);
	return property;
}

}	// End of namespace
}	// End of namespace
