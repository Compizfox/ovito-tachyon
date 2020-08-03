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
#include <ovito/stdobj/properties/PropertyReference.h>
#include <ovito/stdobj/properties/PropertyContainer.h>
#include <ovito/stdobj/properties/PropertyContainerClass.h>
#include <ovito/stdobj/properties/PropertyContainerImportData.h>

namespace Ovito { namespace StdObj {

/**
 * \brief Defines the mapping of a column in an imported data file to a target property in OVITO.
 *
 * An InputColumnMapping is composed of a list of these structures, one for each
 * column in the data file.
 */
class OVITO_STDOBJ_EXPORT InputColumnInfo
{
public:

	/// \brief Default constructor mapping the column to no property.
	InputColumnInfo() = default;

	/// \brief Constructor mapping the column to a standard property.
	InputColumnInfo(PropertyContainerClassPtr pclass, int typeId, int vectorComponent = 0) {
		mapStandardColumn(pclass, typeId, vectorComponent);
	}

	/// \brief Constructor mapping the column to a user-defined property.
	InputColumnInfo(PropertyContainerClassPtr pclass, const QString& propertyName, int dataType, int vectorComponent = 0) {
		mapCustomColumn(pclass, propertyName, dataType, vectorComponent);
	}

	/// \brief Maps this column to a user-defined property.
	/// \param propertyName The name of target property.
	/// \param dataType The data type of the property to create.
	/// \param vectorComponent The component index if the target property is a vector property.
	void mapCustomColumn(PropertyContainerClassPtr pclass, const QString& propertyName, int dataType, int vectorComponent = 0) {
		OVITO_ASSERT(pclass);
		this->property = PropertyReference(pclass, propertyName, vectorComponent);
		this->dataType = dataType;
	}

	/// \brief Maps this column to a standard property.
	/// \param typeId Specifies the standard property type.
	/// \param vectorComponent The component index if the target property is a vector property.
	void mapStandardColumn(PropertyContainerClassPtr pclass, int typeId, int vectorComponent = 0) {
		OVITO_ASSERT(pclass);
		OVITO_ASSERT(typeId != PropertyStorage::GenericUserProperty);
		this->property = PropertyReference(pclass, typeId, vectorComponent);
		this->dataType = pclass->standardPropertyDataType(typeId);
	}

	/// \brief Returns true if the file column is mapped to a target property; false otherwise (file column will be ignored during import).
	bool isMapped() const { return dataType != QMetaType::Void; }

	/// \brief Compares two column records for equality.
	bool operator==(const InputColumnInfo& other) const {
		return property == other.property && dataType == other.dataType && columnName == other.columnName;
	}

	/// The target property this column is mapped to.
	PropertyReference property;

	/// The data type of the property if this column is mapped to a user-defined property.
	/// This field can be set to QMetaType::Void to indicate that the column should be ignored during file import.
	int dataType = QMetaType::Void;

	/// The name of the column in the input file. This information is read from the input file (if available).
	QString columnName;
};

/**
 * \brief Defines a mapping between the columns in a column-based input particle file
 *        and OVITO's internal particle properties.
 */
class OVITO_STDOBJ_EXPORT InputColumnMapping : public std::vector<InputColumnInfo>
{
public:

	/// Default constructor.
	InputColumnMapping() = default;

	/// Constructor.
	explicit InputColumnMapping(PropertyContainerClassPtr containerClass) : _containerClass(containerClass) {
		OVITO_ASSERT(containerClass->isDerivedFrom(PropertyContainer::OOClass()));
	}

	/// Return the class of the property container.
	PropertyContainerClassPtr containerClass() const { return _containerClass; }

	/// \brief Saves the mapping into a byte array.
	QByteArray toByteArray(TaskManager& taskManager) const;

	/// \brief Loads the mapping from a byte array.
	void fromByteArray(const QByteArray& array, TaskManager& taskManager);

	/// \brief Checks if the mapping is valid; throws an exception if not.
	void validate() const;

	/// \brief Returns the first few lines of the file, which can help the user to figure out
	///        the column mapping.
	const QString& fileExcerpt() const { return _fileExcerpt; }

	/// \brief Stores the first few lines of the file, which can help the user to figure out
	///        the column mapping.
	void setFileExcerpt(const QString& text) { _fileExcerpt = text; }

	/// \brief Returns whether at least some of the file columns have names.
	bool hasFileColumnNames() const {
		return std::any_of(begin(), end(), [](const InputColumnInfo& column) {
			return column.columnName.isEmpty() == false;
		});
	}

	/// \brief Maps a file column to a standard property unless there is already another column mapped to the same property.
	/// \param column The file column index.
	/// \param typeId Specifies the standard property.
	/// \param vectorComponent The component index if the target property is a vector property.
	bool mapStandardColumn(int column, int typeId, int vectorComponent = 0);

	/// \brief Maps this column to a user-defined property unless there is already another column mapped to the same property.
	/// \param column The file column index.
	/// \param propertyName The name of target particle property.
	/// \param dataType The data type of the property to create.
	/// \param vectorComponent The component index if the target property is a vector property.
	bool mapCustomColumn(int column, const QString& propertyName, int dataType, int vectorComponent = 0);

	/// \brief Compares two mapping for equality.
	bool operator==(const InputColumnMapping& other) const {
		if(containerClass() != other.containerClass()) return false;
		if(size() != other.size()) return false;
		if(!std::equal(begin(), end(), other.begin())) return false;
		return fileExcerpt() == other.fileExcerpt();
	}

	/// \brief Compares two mappings for inequality.
	bool operator!=(const InputColumnMapping& other) const { return !(*this == other); }

	friend OVITO_STDOBJ_EXPORT SaveStream& operator<<(SaveStream& stream, const InputColumnMapping& m);
	friend OVITO_STDOBJ_EXPORT LoadStream& operator>>(LoadStream& stream, InputColumnMapping& m);

private:

	/// A string with the first few lines of the file, which is meant as a hint for the user to figure out
	/// the correct column mapping.
	QString _fileExcerpt;

	/// The property container type.
	PropertyContainerClassPtr _containerClass = nullptr;
};

/**
 * Encapsulates a column mapping for a specific type of property container.
 */
template<class PropertyContainerType>
class TypedInputColumnMapping : public InputColumnMapping
{
public:

	/// \brief Default constructor.
	TypedInputColumnMapping() : InputColumnMapping(&PropertyContainerType::OOClass()) {}

	/// \brief Constructor for conversion from a generic InputColumnMapping.
	TypedInputColumnMapping(const InputColumnMapping& m) : InputColumnMapping(&PropertyContainerType::OOClass()) {
		OVITO_ASSERT(m.containerClass() == containerClass());
		InputColumnMapping::operator=(m);
	}

	friend SaveStream& operator<<(SaveStream& stream, const TypedInputColumnMapping& m) {
		return stream << static_cast<const InputColumnMapping&>(m);
	}

	friend LoadStream& operator>>(LoadStream& stream, TypedInputColumnMapping& m) {
		return stream >> static_cast<InputColumnMapping&>(m);
	}
};

/**
 * \brief Helper class that reads column-based data from an input file and
 *        stores the parsed values in a set of target properties as specified by an InputColumnMapping.
 */
class OVITO_STDOBJ_EXPORT InputColumnReader : public QObject
{
public:

	/// \brief Initializes the object.
	/// \param mapping Defines the mapping of columns of the input file
	///        to the target properties.
	/// \param destination The container where the parsed data will be stored in.
	/// \param elementCount The number of data elements that will be read from the input file.
	/// \throws Exception if the mapping is not valid.
	InputColumnReader(const InputColumnMapping& mapping, PropertyContainerImportData& destination, size_t elementCount);

	/// \brief Parses the string tokens from one line of the input file and stores the values
	///        in the target properties.
	/// \param elementIndex The line index starting at 0 that specifies the data element whose properties
	///                  are read from the input file.
	/// \param dataLine The text line read from the input file containing the field values.
	void readElement(size_t elementIndex, const char* dataLine);

	/// \brief Parses the string tokens from one line of the input file and stores the values
	///        in the target properties.
	/// \param elementIndex The line index starting at 0 that specifies the data element whose properties
	///                  are read from the input file.
	/// \param dataLine The text line read from the input file containing the field values.
	const char* readElement(size_t elementIndex, const char* dataLine, const char* dataLineEnd);

	/// \brief Processes the values from one line of the input file and stores them in the target properties.
	void readElement(size_t elementIndex, const double* values, int nvalues);

	/// \brief Sorts the created element types either by numeric ID or by name, depending on how they were stored in the input file.
	void sortElementTypes();

private:

	/// Parse a single field from a text line.
	void parseField(size_t elementIndex, int columnIndex, const char* token, const char* token_end);

	/// Specifies the mapping of file columns to target properties.
	InputColumnMapping _mapping;

	/// The container that receives the parsed data.
	PropertyContainerImportData& _destination;

	struct TargetPropertyRecord {
		PropertyPtr property;
		PropertyAccess<void,true> propertyArray;
		uint8_t* data;
		size_t stride;
		size_t count;
		int vectorComponent;
		int dataType;
		PropertyContainerImportData::TypeList* typeList = nullptr;
		bool numericElementTypes;
	};

	/// Specifies the destinations for the column data.
	QVector<TargetPropertyRecord> _properties;
};

}	// End of namespace
}	// End of namespace

Q_DECLARE_METATYPE(Ovito::StdObj::InputColumnInfo);
Q_DECLARE_METATYPE(Ovito::StdObj::InputColumnMapping);
