////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2017 Alexander Stukowski
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
#include <ovito/core/dataset/data/DataObject.h>
#include <ovito/stdobj/properties/ElementType.h>
#include "PropertyStorage.h"

namespace Ovito { namespace StdObj {

/**
 * \brief Stores a property data array.
 */
class OVITO_STDOBJ_EXPORT PropertyObject : public DataObject
{
	Q_OBJECT
	OVITO_CLASS(PropertyObject)

public:

	/// \brief Creates a property object.
	Q_INVOKABLE PropertyObject(DataSet* dataset, const PropertyPtr& storage = nullptr);

	/// \brief Gets the property's name.
	/// \return The name of property, which is shown to the user.
	const QString& name() const { return storage()->name(); }

	/// \brief Sets the property's name.
	/// \param name The new name string.
	void setName(const QString& name);

	/// \brief Returns the number of elements.
	size_t size() const { return storage()->size(); }

	/// \brief Resizes the property storage.
	/// \param newSize The new number of elements.
	/// \param preserveData Controls whether the existing per-element data is preserved.
	///                     This also determines whether newly allocated memory is initialized to zero.
	void resize(size_t newSize, bool preserveData);

	/// \brief Returns the type of this property.
	int type() const { return storage()->type(); }

	/// \brief Returns the data type of the property.
	/// \return The identifier of the data type used for the elements stored in
	///         this property storage according to the Qt meta type system.
	int dataType() const { return storage()->dataType(); }

	/// \brief Returns the number of bytes per value.
	/// \return Number of bytes used to store a single value of the data type
	///         specified by dataType().
	size_t dataTypeSize() const { return storage()->dataTypeSize(); }

	/// \brief Returns the number of bytes used per particle.
	size_t stride() const { return storage()->stride(); }

	/// \brief Returns the number of values per element.
	/// \return The number of data values stored per element in this storage object.
	size_t componentCount() const { return storage()->componentCount(); }

	/// \brief Returns the human-readable names for the components of one element.
	/// \return The names of the components if this property contains more than one value per element.
	///         The list may be empty for scalar properties.
	const QStringList& componentNames() const { return storage()->componentNames(); }

	/// \brief Returns the display name of the property including the name of the given
	///        vector component.
	QString nameWithComponent(int vectorComponent) const {
		if(componentCount() <= 1 || vectorComponent < 0)
			return name();
		else if(vectorComponent < componentNames().size())
			return QStringLiteral("%1.%2").arg(name()).arg(componentNames()[vectorComponent]);
		else
			return QStringLiteral("%1.%2").arg(name()).arg(vectorComponent + 1);
	}

	/// Returns the data encapsulated by this object after making sure it is not shared with other owners.
	const PropertyPtr& modifiableStorage();

	/// Extends the data array and replicates the existing data N times.
	void replicate(size_t n, bool replicateValues = true);

	/// Reduces the size of the storage array, removing elements for which
	/// the corresponding bits in the bit array are set.
	void filterResize(const boost::dynamic_bitset<>& mask) {
		modifiableStorage()->filterResize(mask);
		notifyTargetChanged();
	}

	/// \brief Returns a read-only pointer to the elements stored in this property object.
	template<typename T>
	const T* cdata() const {
		return storage()->cdata<T>();
	}

	/// \brief Returns a read-only pointer to the i-th element stored in this property object.
	template<typename T>
	const T* cdata(size_t i) const {
		return storage()->cdata<T>(i);
	}

	/// Returns a read-only pointer to the j-components of the i-th element in the property array.
	template<typename T>
	const T* cdata(size_t i, size_t j) const {
		return storage()->cdata<T>(i, j);
	}

	/// \brief Returns the value of the i-th element from the array.
	template<typename T>
	const T& get(size_t i) const {
		return storage()->get<T>(i);
	}

	/// \brief Returns the value of the j-th component of the i-th element from the array.
	template<typename T>
	const T& get(size_t i, size_t j) const {
		return storage()->get<T>(i, j);
	}

	/// \brief Returns a range of const iterators over the elements stored in this array.
	template<typename T>
	boost::iterator_range<const T*> crange() const {
		return storage()->crange<T>();
	}

	/// Returns a read-write pointer to the raw elements in the property storage.
	void* data() {
		return modifiableStorage()->data();
	}

	/// \brief Returns a read-write pointer to the first integer element stored in this object.
	/// \note This method may only be used if this property is of data type integer.
	int* dataInt() {
		return modifiableStorage()->dataInt();
	}

	/// \brief Returns a read-write pointer to the first integer element stored in this object.
	/// \note This method may only be used if this property is of data type 64-bit integer.
	qlonglong* dataInt64() {
		return modifiableStorage()->dataInt64();
	}

	/// \brief Returns a read-only pointer to the first float element in the property storage.
	/// \note This method may only be used if this property is of data type float.
	FloatType* dataFloat() {
		return modifiableStorage()->dataFloat();
	}

	/// \brief Returns a read-write pointer to the first vector element in the property storage.
	/// \note This method may only be used if this property is of data type Vector3 or a FloatType channel with 3 components.
	Vector3* dataVector3() {
		return modifiableStorage()->dataVector3();
	}

	/// \brief Returns a read-write pointer to the first point element in the property storage.
	/// \note This method may only be used if this property is of data type Point3 or a FloatType channel with 3 components.
	Point3* dataPoint3() {
		return modifiableStorage()->dataPoint3();
	}

	/// \brief Returns a read-write pointer to the first point element in the property storage.
	/// \note This method may only be used if this property is of data type Vector3I or an integer channel with 3 components.
	Vector3I* dataVector3I() {
		return modifiableStorage()->dataVector3I();
	}

	/// \brief Returns a read-write pointer to the first point element in the property storage.
	/// \note This method may only be used if this property is of data type Point3I or an integer channel with 3 components.
	Point3I* dataPoint3I() {
		return modifiableStorage()->dataPoint3I();
	}

	/// \brief Returns a read-write pointer to the first point element in the property storage.
	/// \note This method may only be used if this property is of data type Color or a FloatType channel with 3 components.
	Color* dataColor() {
		return modifiableStorage()->dataColor();
	}

	/// \brief Returns a read-write pointer to the first symmetric tensor element in the property storage.
	/// \note This method may only be used if this property is of data type SymmetricTensor2 or a FloatType channel with 6 components.
	SymmetricTensor2* dataSymmetricTensor2() {
		return modifiableStorage()->dataSymmetricTensor2();
	}

	/// \brief Returns a read-write pointer to the first quaternion element in the property storage.
	/// \note This method may only be used if this property is of data type Quaternion or a FloatType channel with 4 components.
	Quaternion* dataQuaternion() {
		return modifiableStorage()->dataQuaternion();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<int*> intRange() {
		return modifiableStorage()->intRange();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<qlonglong*> int64Range() {
		return modifiableStorage()->int64Range();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<FloatType*> floatRange() {
		return modifiableStorage()->floatRange();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<Point3*> point3Range() {
		return modifiableStorage()->point3Range();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<Vector3*> vector3Range() {
		return modifiableStorage()->vector3Range();
	}

	/// \brief Returns a range of const iterators over the elements stored in this object.
	boost::iterator_range<Color*> colorRange() {
		return modifiableStorage()->colorRange();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<Vector3I*> vector3IRange() {
		return modifiableStorage()->vector3IRange();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<Point3I*> point3IRange() {
		return modifiableStorage()->point3IRange();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<SymmetricTensor2*> symmetricTensor2Range() {
		return modifiableStorage()->symmetricTensor2Range();
	}

	/// \brief Returns a range of iterators over the elements stored in this object.
	boost::iterator_range<Quaternion*> quaternionRange() {
		return modifiableStorage()->quaternionRange();
	}
	
	/// Sets the value of an integer element at the given index (if this is an integer property).
	void setInt(size_t index, int newValue) {
		modifiableStorage()->setInt(index, newValue);
	}

	/// Sets the value of an integer element at the given index (if this is a 64-bit integer property).
	void setInt64(size_t index, qlonglong newValue) {
		modifiableStorage()->setInt64(index, newValue);
	}

	/// Sets the value of a float element at the given index (if this is a float property).
	void setFloat(size_t index, FloatType newValue) {
		modifiableStorage()->setFloat(index, newValue);
	}

	/// Sets the value of an integer element at the given index (if this is an integer property).
	void setIntComponent(size_t index, size_t componentIndex, int newValue) {
		modifiableStorage()->setIntComponent(index, componentIndex, newValue);
	}

	/// Sets the value of an integer element at the given index (if this is a 64-bit integer property).
	void setInt64Component(size_t index, size_t componentIndex, qlonglong newValue) {
		modifiableStorage()->setInt64Component(index, componentIndex, newValue);
	}

	/// Sets the value of a float element at the given index (if this is a float property).
	void setFloatComponent(size_t index, size_t componentIndex, FloatType newValue) {
		modifiableStorage()->setFloatComponent(index, componentIndex, newValue);
	}

	/// Sets the value of a Vector3 element at the given index (if this is a vector property).
	void setVector3(size_t index, const Vector3& newValue) {
		modifiableStorage()->setVector3(index, newValue);
	}

	/// Sets the value of a Point3 element at the given index (if this is a point property).
	void setPoint3(size_t index, const Point3& newValue) {
		modifiableStorage()->setPoint3(index, newValue);
	}

	/// Sets the value of a Vector3I element at the given index (if this is a point property).
	void setVector3I(size_t index, const Vector3I& newValue) {
		modifiableStorage()->setVector3I(index, newValue);
	}

	/// Sets the value of a Point3I element at the given index (if this is a point property).
	void setPoint3I(size_t index, const Point3I& newValue) {
		modifiableStorage()->setPoint3I(index, newValue);
	}

	/// Sets the value of a Color element at the given index (if this is a point property).
	void setColor(size_t index, const Color& newValue) {
		modifiableStorage()->setColor(index, newValue);
	}

	/// Sets the given element's value of a SymmetricTensor2 property.
	void setSymmetricTensor2(size_t index, const SymmetricTensor2& newValue) {
		modifiableStorage()->setSymmetricTensor2(index, newValue);
	}

	/// Sets the given element's value of a Matrix3 property.
	void setMatrix3(size_t index, const Matrix3& newValue) {
		modifiableStorage()->setMatrix3(index, newValue);
	}

	/// Sets the given element's value of a Quaternion property.
	void setQuaternion(size_t index, const Quaternion& newValue) {
		modifiableStorage()->setQuaternion(index, newValue);
	}

	/// Copies the elements from the given source into this property array using a element mapping.
	void mappedCopyFrom(const PropertyObject* source, const std::vector<size_t>& mapping) {
		modifiableStorage()->mappedCopy(*source->storage(), mapping);
	}

	/// Copies the elements from this property array into the given destination array using an index mapping.
	void mappedCopyTo(PropertyObject* destination, const std::vector<size_t>& mapping) const {
		storage()->mappedCopyTo(*destination->modifiableStorage(), mapping);
	}

	//////////////////////////////// Element types //////////////////////////////

	/// Appends an element type to the list of types.
	void addElementType(const ElementType* type) {
		OVITO_ASSERT(elementTypes().contains(const_cast<ElementType*>(type)) == false);
		_elementTypes.push_back(this, PROPERTY_FIELD(elementTypes), type);
	}

	/// Inserts an element type into the list of types.
	void insertElementType(int index, const ElementType* type) {
		OVITO_ASSERT(elementTypes().contains(const_cast<ElementType*>(type)) == false);
		_elementTypes.insert(this, PROPERTY_FIELD(elementTypes), index, type);
	}

	/// Returns the element type with the given ID, or NULL if no such type exists.
	ElementType* elementType(int id) const {
		for(ElementType* type : elementTypes())
			if(type->numericId() == id)
				return type;
		return nullptr;
	}

	/// Returns the element type with the given human-readable name, or NULL if no such type exists.
	ElementType* elementType(const QString& name) const {
		OVITO_ASSERT(!name.isEmpty());
		for(ElementType* type : elementTypes())
			if(type->name() == name)
				return type;
		return nullptr;
	}

	/// Removes a single element type from this object.
	void removeElementType(int index) {
		_elementTypes.remove(this, PROPERTY_FIELD(elementTypes), index);
	}

	/// Removes all elements types from this object.
	void clearElementTypes() {
		_elementTypes.clear(this, PROPERTY_FIELD(elementTypes));
	}

	/// Builds a mapping from numeric IDs to type colors.
	std::map<int,Color> typeColorMap() const {
		std::map<int,Color> m;
		for(ElementType* type : elementTypes())
			m.insert({type->numericId(), type->color()});
		return m;
	}

	/// Returns an numeric type ID that is not yet used by any of the existing element types.
	int generateUniqueElementTypeId(int startAt = 1) const {
		int maxId = startAt;
		for(ElementType* type : elementTypes())
			maxId = std::max(maxId, type->numericId() + 1);
		return maxId;
	}

	/// Helper method that remaps the existing type IDs to a contiguous range starting at the given
	/// base ID. This method is mainly used for file output, because some file formats
	/// work with numeric particle types only, which must form a contiguous range.
	/// The method returns the mapping of output type IDs to original type IDs
	/// and a copy of the property array in which the original type ID values have
	/// been remapped to the output IDs.
	std::tuple<std::map<int,int>, ConstPropertyPtr> generateContiguousTypeIdMapping(int baseId = 1) const;

	////////////////////////////// Support functions for the Python bindings //////////////////////////////

	/// Indicates to the Python binding layer that this property object has been temporarily put into a
	/// writable state. In this state, the binding layer will allow write access to the property's internal data.
	bool isWritableFromPython() const { return _isWritableFromPython != 0; }

	/// Puts the property array into a writable state.
	/// In the writable state, the Python binding layer will allow write access to the property's internal data.
	void makeWritableFromPython();

	/// Puts the property array back into the default read-only state.
	/// In the read-only state, the Python binding layer will not permit write access to the property's internal data.
	void makeReadOnlyFromPython() {
		OVITO_ASSERT(_isWritableFromPython > 0);
		_isWritableFromPython--;
	}

	/// Returns whether this data object wants to be shown in the pipeline editor
	/// under the data source section.
	/// This implementation returns true only it this is a typed property, i.e. if the 'elementTypes' list contains
	/// some elements. In this case we want the property to appear in the pipeline editor so that the user can
	/// edit the individual types.
	virtual bool showInPipelineEditor() const override {
		return !elementTypes().empty();
	}

	/// Returns the display title of this property object in the user interface.
	virtual QString objectTitle() const override;

protected:

	/// Saves the class' contents to the given stream.
	virtual void saveToStream(ObjectSaveStream& stream, bool excludeRecomputableData) override;

	/// Loads the class' contents from the given stream.
	virtual void loadFromStream(ObjectLoadStream& stream) override;

private:

	/// The internal per-element data.
	DECLARE_RUNTIME_PROPERTY_FIELD(PropertyPtr, storage, setStorage);

	/// Contains the list of defined "types" if this is a typed property.
	DECLARE_MODIFIABLE_VECTOR_REFERENCE_FIELD(ElementType, elementTypes, setElementTypes);

	/// The user-interface title of this property.
	DECLARE_MODIFIABLE_PROPERTY_FIELD(QString, title, setTitle);

	/// This is a special flag used by the Python bindings to indicate that
	/// this property object has been temporarily put into a writable state.
	int _isWritableFromPython = 0;
};

}	// End of namespace
}	// End of namespace
