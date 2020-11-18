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


#include <ovito/particles/Particles.h>
#include <ovito/grid/objects/VoxelGrid.h>
#include <ovito/core/dataset/io/FileSourceImporter.h>
#include <ovito/stdobj/simcell/SimulationCell.h>
#include <ovito/stdobj/properties/PropertyContainerImportData.h>

namespace Ovito { namespace Particles {

/**
 * Holds the data of a single frame loaded by a ParticleImporter.
 */
class OVITO_PARTICLES_EXPORT ParticleFrameData : public FileSourceImporter::FrameData
{
public:

	/// Constructor.
	ParticleFrameData();

	/// Inserts the loaded data into the provided pipeline state structure. This function is
	/// called by the system from the main thread after the asynchronous loading task has finished.
	virtual OORef<DataCollection> handOver(const DataCollection* existing, bool isNewFile, CloneHelper& cloneHelper, FileSource* fileSource) override;

	/// Returns the current simulation cell matrix.
	const SimulationCell& simulationCell() const { return _simulationCell; }

	/// Returns a reference to the simulation cell.
	SimulationCell& simulationCell() { return _simulationCell; }

	/// Returns the per-particle data.
	PropertyContainerImportData& particles() { return _particleData; }

	/// Returns the per-bond data.
	PropertyContainerImportData& bonds() { return _bondData; }

	/// Determines the PBC shift vectors for bonds using the minimum image convention.
	void generateBondPeriodicImageProperty();

	/// Returns the per-angle data.
	PropertyContainerImportData& angles() { return _angleData; }

	/// Returns the per-dihedral data.
	PropertyContainerImportData& dihedrals() { return _dihedralData; }

	/// Returns the per-improper data.
	PropertyContainerImportData& impropers() { return _improperData; }

	/// Returns the per-voxel data.
	PropertyContainerImportData& voxels() { return _voxelData; }

	/// Returns the shape of the voxel grid.
	const VoxelGrid::GridDimensions& voxelGridShape() const { return _voxelGridShape; }

	/// Sets the shape of the voxel grid.
	void setVoxelGridShape(const VoxelGrid::GridDimensions& shape) { _voxelGridShape = shape; }

	/// Returns the human-readable name being assigned to the loaded voxel grid.
	const QString& voxelGridTitle() const { return _voxelGridTitle; }

	/// Sets the human-readable name that will be assigned to the voxel grid.
	void setVoxelGridTitle(const QString& title) { _voxelGridTitle = title; }

	/// Returns the unique data object ID being assigned to the loaded voxel grid.
	const QString& voxelGridId() const { return _voxelGridId; }

	/// Sets the unique data object ID that will be assigned to the voxel grid.
	void setVoxelGridId(const QString& id) { _voxelGridId = id; }

	/// Returns the metadata read from the file header.
	QVariantMap& attributes() { return _attributes; }

	/// Parsers call this method to indicate that the input file contains
	/// additional frames stored back to back with the currently loaded one.
	void signalAdditionalFrames() { _detectedAdditionalFrames = true; }

	/// Sorts the particles list with respect to particle IDs.
	/// Does nothing if particles do not have IDs.
	void sortParticlesById();

private:

	/// The simulation cell.
	SimulationCell _simulationCell;

	/// Particle properties.
	PropertyContainerImportData _particleData;

	/// Bond properties.
	PropertyContainerImportData _bondData;

	/// Angle properties.
	PropertyContainerImportData _angleData;

	/// Dihedral properties.
	PropertyContainerImportData _dihedralData;

	/// Improper properties.
	PropertyContainerImportData _improperData;

	/// Voxel grid properties.
	PropertyContainerImportData _voxelData;

	/// The dimensions of the voxel grid.
	VoxelGrid::GridDimensions _voxelGridShape{{0,0,0}};

	/// The human-readable name to assign to the loaded voxel grid.
	QString _voxelGridTitle;

	/// The ID to assign to the voxel grid data object.
	QString _voxelGridId;
	
	/// The metadata read from the file header.
	QVariantMap _attributes;

	/// Flag that is set by the parser to indicate that the input file contains more than one frame.
	bool _detectedAdditionalFrames = false;
};

}	// End of namespace
}	// End of namespace
