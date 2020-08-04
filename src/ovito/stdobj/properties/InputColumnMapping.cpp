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
#include <ovito/core/utilities/io/NumberParsing.h>
#include "InputColumnMapping.h"

namespace Ovito { namespace StdObj {

/******************************************************************************
 * Maps a file column to a standard property unless there is already another 
 * column mapped to the same property.
 *****************************************************************************/
bool InputColumnMapping::mapStandardColumn(int column, int typeId, int vectorComponent) 
{
	OVITO_ASSERT(column >= 0 && column < this->size());
	OVITO_ASSERT(typeId != PropertyStorage::GenericUserProperty);
	OVITO_ASSERT(containerClass());

	// Check if there is another file column already mapped to the same target property.
	for(const InputColumnInfo& columnInfo : *this) {
		if(columnInfo.property.type() == typeId && columnInfo.property.vectorComponent() == vectorComponent)
			return false;
	}

	// If not, record the mapping.
	(*this)[column].mapStandardColumn(containerClass(), typeId, vectorComponent);
	return true;
}

/******************************************************************************
 * Maps this column to a user-defined property unless there is already another 
 * column mapped to the same property.
 *****************************************************************************/
bool InputColumnMapping::mapCustomColumn(int column, const QString& propertyName, int dataType, int vectorComponent) 
{
	OVITO_ASSERT(column >= 0 && column < this->size());
	OVITO_ASSERT(containerClass());

	// Check if there is another file column already mapped to the same target property.
	for(const InputColumnInfo& columnInfo : *this) {
		if(columnInfo.property.type() == PropertyStorage::GenericUserProperty && columnInfo.property.name() == propertyName && columnInfo.property.vectorComponent() == vectorComponent)
			return false;
	}

	// If not, record the mapping.
	(*this)[column].mapCustomColumn(containerClass(), propertyName, dataType, vectorComponent);
	return true;
}

/******************************************************************************
 * Saves the mapping to the given stream.
 *****************************************************************************/
SaveStream& operator<<(SaveStream& stream, const InputColumnMapping& m)
{
	stream.beginChunk(0x02);
	stream << m.containerClass();
	stream.writeSizeT(m.size());
	for(const InputColumnInfo& col : m) {
		stream << col.property;
		stream << col.columnName;
		stream << col.dataType;
	}
	stream.endChunk();
	return stream;
}

/******************************************************************************
 * Loads the mapping from the given stream.
 *****************************************************************************/
LoadStream& operator>>(LoadStream& stream, InputColumnMapping& m)
{
	int version = stream.expectChunkRange(0x0, 0x02);

	// For backward compatibility with OVITO 3.1:
	if(version == 1) {
		int numColumns;
		stream >> numColumns;
		m.resize(numColumns);
		for(InputColumnInfo& col : m) {
			stream >> col.columnName;
			int propertyType;
			stream >> propertyType;
			QString propertyName;
			stream >> propertyName;
			stream >> col.dataType;
			if(col.dataType == qMetaTypeId<float>() || col.dataType == qMetaTypeId<double>())
				col.dataType = PropertyStorage::Float;
			int vectorComponent;
			stream >> vectorComponent;
			if(col.dataType != QMetaType::Void) {
				if(propertyType == PropertyStorage::GenericUserProperty)
					col.property = PropertyReference(m.containerClass(), propertyName, vectorComponent);
				else
					col.property = PropertyReference(m.containerClass(), propertyType, vectorComponent);
			}
		}
	}
	else {
		stream >> m._containerClass;
		m.resize(stream.readSizeT());
		for(InputColumnInfo& col : m) {
			stream >> col.property;
			stream >> col.columnName;
			stream >> col.dataType;
			if(col.dataType == qMetaTypeId<float>() || col.dataType == qMetaTypeId<double>())
				col.dataType = PropertyStorage::Float;
		}
	}
	stream.closeChunk();
	return stream;
}

/******************************************************************************
 * Saves the mapping into a byte array.
 *****************************************************************************/
QByteArray InputColumnMapping::toByteArray(TaskManager& taskManager) const
{
	QByteArray buffer;
	QDataStream dstream(&buffer, QIODevice::WriteOnly);
	SaveStream stream(dstream, SynchronousOperation::createSignal(taskManager));
	stream << *this;
	stream.close();
	return buffer;
}

/******************************************************************************
 * Loads the mapping from a byte array.
 *****************************************************************************/
void InputColumnMapping::fromByteArray(const QByteArray& array, TaskManager& taskManager)
{
	QDataStream dstream(array);
	LoadStream stream(dstream, SynchronousOperation::createSignal(taskManager));
	stream >> *this;
	stream.close();
}

/******************************************************************************
 * Checks if the mapping is valid; throws an exception if not.
 *****************************************************************************/
void InputColumnMapping::validate() const
{
	OVITO_ASSERT(containerClass());

	// Let the property container class perform custom checks.
	containerClass()->validateInputColumnMapping(*this);

	// Check for conflicting mappings, i.e. several file columns being mapped to the same particle property.
	int numMapped = 0;
	for(auto m1 = begin(); m1 != end(); ++m1) {
		if(!m1->isMapped()) continue;
		numMapped++;
		OVITO_ASSERT(m1->property.containerClass() == containerClass());
		for(auto m2 = std::next(m1); m2 != end(); ++m2) {
			if(m1->property == m2->property)
				throw Exception(InputColumnReader::tr("Invalid file column mapping: File columns %1 and %2 cannot both be mapped to the same property '%3'.")
					.arg(std::distance(begin(), m1) + 1)
					.arg(std::distance(begin(), m2) + 1)
					.arg(m1->property.nameWithComponent()));
		}
	}

	if(numMapped == 0)
		throw Exception(InputColumnReader::tr("File column mapping is empty. Please specify how data columns of the input file should be mapped to the properties of %1.").arg(containerClass()->elementDescriptionName()));
}

/******************************************************************************
 * Initializes the object.
 *****************************************************************************/
InputColumnReader::InputColumnReader(const InputColumnMapping& mapping, PropertyContainerImportData& destination, size_t elementCount)
	: _mapping(mapping), _destination(destination)
{
	mapping.validate();

	// Create target properties as defined by the mapping.
	for(int i = 0; i < (int)mapping.size(); i++) {

		PropertyPtr property;
		const PropertyReference& pref = mapping[i].property;

		int vectorComponent = std::max(0, pref.vectorComponent());
		int dataType = mapping[i].dataType;

		TargetPropertyRecord rec;

		if(dataType != QMetaType::Void) {
			if(dataType != PropertyStorage::Int && dataType != PropertyStorage::Int64 && dataType != PropertyStorage::Float)
				throw Exception(tr("Invalid user-defined target property (data type %1) for input file column %2").arg(dataType).arg(i+1));

			if(pref.type() != PropertyStorage::GenericUserProperty) {
				// Look for existing standard property.
				for(const auto& p : destination.properties()) {
					if(p->type() == pref.type()) {
						property = p;
						break;
					}
				}
				if(!property) {
					// Create standard property.
					property = pref.containerClass()->createStandardStorage(elementCount, pref.type(), true);
					destination.addProperty(property);

					// Also create a type list if it is a typed property.
					if(OvitoClassPtr elementTypeClass = pref.containerClass()->typedPropertyElementClass(pref.type()))
						rec.typeList = destination.createPropertyTypesList(property, *elementTypeClass);
				}
			}
			else {
				// Look for existing user-defined property with the same name.
                PropertyPtr oldProperty = nullptr;
				for(int j = 0; j < (int)destination.properties().size(); j++) {
					const auto& p = destination.properties()[j];
					if(p->name() == pref.name()) {
						if(p->dataType() == dataType && (int)p->componentCount() > vectorComponent) {
							property = p;
                        }
						else {
                            oldProperty = p;
                        }
						break;
					}
				}
				if(!property) {
					// Create a new user-defined property for the column.
					property = std::make_shared<PropertyStorage>(elementCount, dataType, vectorComponent + 1, 0, pref.name(), true);
					destination.addProperty(property);
					if(oldProperty) {
						// We need to replace all old properties (with lower vector component count) with this one.
						for(TargetPropertyRecord& rec2 : _properties) {
							if(rec2.property == oldProperty)
								rec2.property = property;
						}
						// Remove old property.
						destination.removeProperty(oldProperty);
					}
				}
			}
			if(property)
				property->setName(pref.name());

			OVITO_ASSERT(vectorComponent < (int)property->componentCount());
			rec.vectorComponent = vectorComponent;
		}

		// Build list of target properties for fast look up during parsing.
		rec.property = property;
		_properties.push_back(rec);
	}

	// Finalize the target property records.
	for(TargetPropertyRecord& rec : _properties) {
		if(rec.property) {
			rec.count = rec.property->size();
			rec.numericElementTypes = true;
			rec.dataType = rec.property->dataType();
			rec.stride = rec.property->stride();
			rec.propertyArray = PropertyAccess<void,true>(rec.property);
			rec.data = reinterpret_cast<uint8_t*>(rec.propertyArray.data(rec.vectorComponent));
		}
	}
}

/******************************************************************************
 * Parses the string tokens from one line of the input file and stores the values
 * in the target properties.
 *****************************************************************************/
const char* InputColumnReader::readElement(size_t elementIndex, const char* s, const char* s_end)
{
	OVITO_ASSERT(_properties.size() == _mapping.size());
	OVITO_ASSERT(s <= s_end);

	int columnIndex = 0;
	while(columnIndex < _properties.size()) {
		// Skip initial whitespace.
		while(s != s_end && (*s == ' ' || *s == '\t' || *s == '\r'))
			++s;
		if(s == s_end || *s == '\n') break;
		const char* token = s;
		// Go to end of token.
		while(s != s_end && (*s > ' ' || *s < 0))
			++s;
		if(s != token) {
			parseField(elementIndex, columnIndex, token, s);
			columnIndex++;
		}
		if(s == s_end) break;
	}
	if(columnIndex < _properties.size())
		throw Exception(tr("Data line in input file does not contain enough columns. Expected %1 file columns, but found only %2.").arg(_properties.size()).arg(columnIndex));

	// Skip to end of line.
	while(s != s_end && *s != '\n')
		++s;
	if(s != s_end) ++s;
	return s;
}

/******************************************************************************
 * Parses the string tokens from one line of the input file and stores the values
 * in the target properties.
 *****************************************************************************/
void InputColumnReader::readElement(size_t elementIndex, const char* s)
{
	OVITO_ASSERT(_properties.size() == _mapping.size());

	int columnIndex = 0;
	while(columnIndex < _properties.size()) {
		while(*s == ' ' || *s == '\t')
			++s;
		const char* token = s;
		while(*s > ' ' || *s < 0)
			++s;
		if(s != token) {
			parseField(elementIndex, columnIndex, token, s);
			columnIndex++;
		}
		if(*s == '\0') break;
		s++;
	}
	if(columnIndex < _properties.size())
		throw Exception(tr("Data line in input file does not contain enough columns. Expected %1 file columns, but found only %2.").arg(_properties.size()).arg(columnIndex));
}

/******************************************************************************
 * Parse a single field from a text line.
 *****************************************************************************/
void InputColumnReader::parseField(size_t elementIndex, int columnIndex, const char* token, const char* token_end)
{
	TargetPropertyRecord& prec = _properties[columnIndex];
	if(!prec.property || !prec.data) return;

	if(elementIndex >= prec.count)
		throw Exception(tr("Too many data lines in input file. Expected only %1 lines.").arg(prec.count));

	if(prec.dataType == PropertyStorage::Float) {
		if(!parseFloatType(token, token_end, *reinterpret_cast<FloatType*>(prec.data + elementIndex * prec.stride)))
			throw Exception(tr("Invalid floating-point value in column %1 (%2): \"%3\"").arg(columnIndex+1).arg(prec.property->name()).arg(QString::fromLocal8Bit(token, token_end - token)));
	}
	else if(prec.dataType == PropertyStorage::Int) {
		int& d = *reinterpret_cast<int*>(prec.data + elementIndex * prec.stride);
		bool ok = parseInt(token, token_end, d);
		if(prec.typeList == nullptr) {
			if(!ok) {
				ok = parseBool(token, token_end, d);
				if(!ok)
					throw Exception(tr("Invalid integer/bool value in column %1 (%2): \"%3\"").arg(columnIndex+1).arg(prec.property->name()).arg(QString::fromLocal8Bit(token, token_end - token)));
			}
		}
		else {
			// Automatically register a new element type if a new type identifier is encountered.
			if(ok) {
				prec.typeList->addTypeId(d);
			}
			else {
				d = prec.typeList->addTypeName(token, token_end);
				prec.numericElementTypes = false;
			}
		}
	}
	else if(prec.dataType == PropertyStorage::Int64) {
		qlonglong& d = *reinterpret_cast<qlonglong*>(prec.data + elementIndex * prec.stride);
		if(!parseInt64(token, token_end, d))
			throw Exception(tr("Invalid 64-bit integer value in column %1 (%2): \"%3\"").arg(columnIndex+1).arg(prec.property->name()).arg(QString::fromLocal8Bit(token, token_end - token)));
	}
}

/******************************************************************************
 * Processes the values from one line of the input file and stores them
 * in the target properties.
 *****************************************************************************/
void InputColumnReader::readElement(size_t elementIndex, const double* values, int nvalues)
{
	OVITO_ASSERT(_properties.size() == _mapping.size());
	if(nvalues < _properties.size())
		throw Exception(tr("Data record in input file does not contain enough columns. Expected %1 file columns, but found only %2.").arg(_properties.size()).arg(nvalues));

	auto prec = _properties.cbegin();
	const double* token = values;

	for(int columnIndex = 0; prec != _properties.cend(); ++columnIndex, ++token, ++prec) {
		if(!prec->property) continue;

		if(elementIndex >= prec->count)
			throw Exception(tr("Too many data lines in input file. Expected only %1 lines.").arg(prec->count));

		if(prec->data) {
			if(prec->dataType == PropertyStorage::Float) {
				*reinterpret_cast<FloatType*>(prec->data + elementIndex * prec->stride) = (FloatType)*token;
			}
			else if(prec->dataType == PropertyStorage::Int) {
				int ival = (int)*token;
				if(prec->typeList) {
					// Automatically register a new element type if a new type identifier is encountered.
					prec->typeList->addTypeId(ival);
				}
				*reinterpret_cast<int*>(prec->data + elementIndex * prec->stride) = ival;
			}
			else if(prec->dataType == PropertyStorage::Int64) {
				*reinterpret_cast<qlonglong*>(prec->data + elementIndex * prec->stride) = (qlonglong)*token;
			}
		}
	}
}

/******************************************************************************
 * Sorts the created element types either by numeric ID or by name,
 * depending on how they were stored in the input file.
 *****************************************************************************/
void InputColumnReader::sortElementTypes()
{
	for(const TargetPropertyRecord& p : _properties) {
		if(p.typeList && p.property) {
			// Since we created element types on the go while reading the elements, the assigned type IDs
			// depend on the storage order of data elements in the file. We rather want a well-defined type ordering, that's
			// why we sort them here according to their names or numeric IDs.
			if(p.numericElementTypes) {
				p.typeList->sortTypesById();
			}
			else {
				PropertyAccess<int> typeArray(p.property);
				p.typeList->sortTypesByName(typeArray);
			}
			break;
		}
	}
}

}	// End of namespace
}	// End of namespace
