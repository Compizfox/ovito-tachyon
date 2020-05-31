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

#include <ovito/particles/Particles.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include "AnglesObject.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(AnglesObject);

/******************************************************************************
* Constructor.
******************************************************************************/
AnglesObject::AnglesObject(DataSet* dataset) : PropertyContainer(dataset)
{
	// Assign the default data object identifier.
	setIdentifier(OOClass().pythonName());
}

/******************************************************************************
* Creates a storage object for standard properties.
******************************************************************************/
PropertyPtr AnglesObject::OOMetaClass::createStandardStorage(size_t elementCount, int type, bool initializeMemory, const ConstDataObjectPath& containerPath) const
{
	int dataType;
	size_t componentCount;
	size_t stride;

	switch(type) {
	case TypeProperty:
		dataType = PropertyStorage::Int;
		componentCount = 1;
		stride = sizeof(int);
		break;
	case TopologyProperty:
		dataType = PropertyStorage::Int64;
		componentCount = 3;
		stride = componentCount * sizeof(qlonglong);
		break;
	default:
		OVITO_ASSERT_MSG(false, "AnglesObject::createStandardStorage", "Invalid standard property type");
		throw Exception(tr("This is not a valid standard angle property type: %1").arg(type));
	}
	const QStringList& componentNames = standardPropertyComponentNames(type);
	const QString& propertyName = standardPropertyName(type);

	OVITO_ASSERT(componentCount == standardPropertyComponentCount(type));

	PropertyPtr property = std::make_shared<PropertyStorage>(elementCount, dataType, componentCount, stride,
								propertyName, false, type, componentNames);

	if(initializeMemory) {
		// Default-initialize property values with zeros.
		property->fillZero();
	}

	return property;
}

/******************************************************************************
* Registers all standard properties with the property traits class.
******************************************************************************/
void AnglesObject::OOMetaClass::initialize()
{
	PropertyContainerClass::initialize();

	setPropertyClassDisplayName(tr("Angles"));
	setElementDescriptionName(QStringLiteral("angles"));
	setPythonName(QStringLiteral("angles"));

	const QStringList emptyList;
	const QStringList abcList = QStringList() << "A" << "B" << "C";

	registerStandardProperty(TypeProperty, tr("Angle Type"), PropertyStorage::Int, emptyList, tr("Angle types"));
	registerStandardProperty(TopologyProperty, tr("Topology"), PropertyStorage::Int64, abcList);
}

}	// End of namespace
}	// End of namespace
