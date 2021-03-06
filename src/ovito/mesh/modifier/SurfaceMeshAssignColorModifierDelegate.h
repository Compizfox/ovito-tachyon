////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2019 Alexander Stukowski
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


#include <ovito/mesh/Mesh.h>
#include <ovito/mesh/surface/SurfaceMesh.h>
#include <ovito/stdmod/modifiers/AssignColorModifier.h>

namespace Ovito { namespace Mesh {

using namespace Ovito::StdMod;

/**
 * \brief Delegate function for the AssignColorModifier that operates on surface mesh vertices.
 */
class OVITO_MESHMOD_EXPORT SurfaceMeshVerticesAssignColorModifierDelegate : public AssignColorModifierDelegate
{
	/// Give the modifier delegate its own metaclass.
	class OOMetaClass : public AssignColorModifierDelegate::OOMetaClass
	{
	public:

		/// Inherit constructor from base class.
		using AssignColorModifierDelegate::OOMetaClass::OOMetaClass;

		/// Indicates which data objects in the given input data collection the modifier delegate is able to operate on.
		virtual QVector<DataObjectReference> getApplicableObjects(const DataCollection& input) const override;

		/// Indicates which class of data objects the modifier delegate is able to operate on.
		virtual const DataObject::OOMetaClass& getApplicableObjectClass() const override { return SurfaceMeshVertices::OOClass(); }

		/// The name by which Python scripts can refer to this modifier delegate.
		virtual QString pythonDataName() const override { return QStringLiteral("surface_vertices"); }
	};

	Q_OBJECT
	OVITO_CLASS_META(SurfaceMeshVerticesAssignColorModifierDelegate, OOMetaClass)

	Q_CLASSINFO("DisplayName", "Mesh Vertices");

public:

	/// Constructor.
	Q_INVOKABLE SurfaceMeshVerticesAssignColorModifierDelegate(DataSet* dataset) : AssignColorModifierDelegate(dataset) {}

protected:

	/// \brief returns the ID of the standard property that will receive the computed colors.
	virtual int outputColorPropertyId() const override { return SurfaceMeshVertices::ColorProperty; }
};

/**
 * \brief Delegate function for the AssignColorModifier that operates on surface mesh faces.
 */
class OVITO_MESHMOD_EXPORT SurfaceMeshFacesAssignColorModifierDelegate : public AssignColorModifierDelegate
{
	/// Give the modifier delegate its own metaclass.
	class OOMetaClass : public AssignColorModifierDelegate::OOMetaClass
	{
	public:

		/// Inherit constructor from base class.
		using AssignColorModifierDelegate::OOMetaClass::OOMetaClass;

		/// Indicates which data objects in the given input data collection the modifier delegate is able to operate on.
		virtual QVector<DataObjectReference> getApplicableObjects(const DataCollection& input) const override;

		/// Indicates which class of data objects the modifier delegate is able to operate on.
		virtual const DataObject::OOMetaClass& getApplicableObjectClass() const override { return SurfaceMeshFaces::OOClass(); }

		/// The name by which Python scripts can refer to this modifier delegate.
		virtual QString pythonDataName() const override { return QStringLiteral("surface_faces"); }
	};

	Q_OBJECT
	OVITO_CLASS_META(SurfaceMeshFacesAssignColorModifierDelegate, OOMetaClass)

	Q_CLASSINFO("DisplayName", "Mesh Faces");

public:

	/// Constructor.
	Q_INVOKABLE SurfaceMeshFacesAssignColorModifierDelegate(DataSet* dataset) : AssignColorModifierDelegate(dataset) {}

protected:

	/// \brief returns the ID of the standard property that will receive the computed colors.
	virtual int outputColorPropertyId() const override { return SurfaceMeshFaces::ColorProperty; }
};

/**
 * \brief Delegate function for the AssignColorModifier that operates on surface mesh regions.
 */
class OVITO_MESHMOD_EXPORT SurfaceMeshRegionsAssignColorModifierDelegate : public AssignColorModifierDelegate
{
	/// Give the modifier delegate its own metaclass.
	class OOMetaClass : public AssignColorModifierDelegate::OOMetaClass
	{
	public:

		/// Inherit constructor from base class.
		using AssignColorModifierDelegate::OOMetaClass::OOMetaClass;

		/// Indicates which data objects in the given input data collection the modifier delegate is able to operate on.
		virtual QVector<DataObjectReference> getApplicableObjects(const DataCollection& input) const override;

		/// Indicates which class of data objects the modifier delegate is able to operate on.
		virtual const DataObject::OOMetaClass& getApplicableObjectClass() const override { return SurfaceMeshRegions::OOClass(); }

		/// The name by which Python scripts can refer to this modifier delegate.
		virtual QString pythonDataName() const override { return QStringLiteral("surface_regions"); }
	};

	Q_OBJECT
	OVITO_CLASS_META(SurfaceMeshRegionsAssignColorModifierDelegate, OOMetaClass)

	Q_CLASSINFO("DisplayName", "Mesh Regions");

public:

	/// Constructor.
	Q_INVOKABLE SurfaceMeshRegionsAssignColorModifierDelegate(DataSet* dataset) : AssignColorModifierDelegate(dataset) {}

protected:

	/// \brief returns the ID of the standard property that will receive the computed colors.
	virtual int outputColorPropertyId() const override { return SurfaceMeshRegions::ColorProperty; }
};
}	// End of namespace
}	// End of namespace
