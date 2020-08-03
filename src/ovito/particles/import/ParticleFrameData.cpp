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
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/particles/objects/ParticlesVis.h>
#include <ovito/particles/objects/BondsObject.h>
#include <ovito/particles/objects/BondsVis.h>
#include <ovito/grid/objects/VoxelGrid.h>
#include <ovito/grid/objects/VoxelGridVis.h>
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include <ovito/stdobj/simcell/SimulationCellVis.h>
#include <ovito/stdobj/properties/PropertyStorage.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include <ovito/core/app/Application.h>
#include <ovito/core/dataset/io/FileSource.h>
#include "ParticleFrameData.h"
#include "ParticleImporter.h"

namespace Ovito { namespace Particles {

/******************************************************************************
* Constructor.
******************************************************************************/
ParticleFrameData::ParticleFrameData() 
{
	// Assume periodic boundary conditions by default.
	_simulationCell.setPbcFlags(true, true, true);

	// Set the default visualization element types.
	particles().setVisElementClass(&ParticlesVis::OOClass());
	bonds().setVisElementClass(&BondsVis::OOClass());
	voxels().setVisElementClass(&VoxelGridVis::OOClass());
}

/******************************************************************************
* Determines the PBC shift vectors for bonds using the minimum image convention.
******************************************************************************/
void ParticleFrameData::generateBondPeriodicImageProperty()
{
	ConstPropertyAccess<Point3> posProperty = particles().findStandardProperty(ParticlesObject::PositionProperty);
	if(!posProperty) return;

	ConstPropertyAccess<ParticleIndexPair> bondTopologyProperty = bonds().findStandardProperty(BondsObject::TopologyProperty);
	if(!bondTopologyProperty) return;

	OVITO_ASSERT(!bonds().findStandardProperty(BondsObject::PeriodicImageProperty));
	PropertyAccess<Vector3I> bondPeriodicImageProperty = bonds().addProperty(BondsObject::OOClass().createStandardStorage(bondTopologyProperty.size(), BondsObject::PeriodicImageProperty, true));

	if(!simulationCell().hasPbc())
		return;

	for(size_t bondIndex = 0; bondIndex < bondTopologyProperty.size(); bondIndex++) {
		size_t index1 = bondTopologyProperty[bondIndex][0];
		size_t index2 = bondTopologyProperty[bondIndex][1];
		OVITO_ASSERT(index1 < posProperty.size() && index2 < posProperty.size());
		Vector3 delta = simulationCell().absoluteToReduced(posProperty[index2] - posProperty[index1]);
		for(size_t dim = 0; dim < 3; dim++) {
			if(simulationCell().hasPbc(dim))
				bondPeriodicImageProperty[bondIndex][dim] = -(int)std::floor(delta[dim] + FloatType(0.5));
		}
	}
}

/******************************************************************************
* Inserts the loaded data into the provided pipeline state structure.
* This function is called by the system from the main thread after the
* asynchronous loading task has finished.
******************************************************************************/
OORef<DataCollection> ParticleFrameData::handOver(const DataCollection* existing, bool isNewFile, CloneHelper& cloneHelper, FileSource* fileSource)
{
	// Start with a fresh data collection that will be populated.
	OORef<DataCollection> output = new DataCollection(fileSource->dataset());

	// Create the simulation cell.
	const SimulationCellObject* existingCell = existing ? existing->getObject<SimulationCellObject>() : nullptr;
	if(!existingCell) {
		// Create a new SimulationCellObject.
		SimulationCellObject* cell = output->createObject<SimulationCellObject>(fileSource, simulationCell());

		// Initialize the simulation cell and its vis element with default values.
		if(Application::instance()->executionContext() == Application::ExecutionContext::Interactive)
			cell->loadUserDefaults();

		// Set up the vis element for the simulation cell.
		if(SimulationCellVis* cellVis = dynamic_object_cast<SimulationCellVis>(cell->visElement())) {

			// Choose an appropriate line width depending on the cell's size.
			FloatType cellDiameter = (
					simulationCell().matrix().column(0) +
					simulationCell().matrix().column(1) +
					simulationCell().matrix().column(2)).length();
			cellVis->setDefaultCellLineWidth(std::max(cellDiameter * FloatType(1.4e-3), FloatType(1e-8)));
			cellVis->setCellLineWidth(cellVis->defaultCellLineWidth());
		}
	}
	else {
		// Adopt pbc flags from input file only if it is a new file.
		// This gives the user the option to change the pbc flags without them
		// being overwritten when a new frame from a simulation sequence is loaded.
		SimulationCellObject* cell = cloneHelper.cloneObject(existingCell, false); 
		output->addObject(cell);
		cell->setData(simulationCell(), isNewFile);
	}

	if(!particles().properties().empty()) {

		// Hand over particles.
		const ParticlesObject* existingParticles = existing ? existing->getObject<ParticlesObject>() : nullptr;
		ParticlesObject* targetParticles = output->createObject<ParticlesObject>(fileSource);
		particles().transferToContainer(existingParticles, targetParticles, isNewFile, cloneHelper);

		// Auto-adjust particle display radius.
		if(isNewFile) {
			if(ParticlesVis* particleVis = targetParticles->visElement<ParticlesVis>()) {
				FloatType cellDiameter = (
						simulationCell().matrix().column(0) +
						simulationCell().matrix().column(1) +
						simulationCell().matrix().column(2)).length();
				// Limit particle radius to a fraction of the cell diameter.
				// This is to avoid extremely large particles when the length scale of the simulation is <<1.
				cellDiameter /= 2;
				if(particleVis->defaultParticleRadius() > cellDiameter && cellDiameter != 0)
					particleVis->setDefaultParticleRadius(cellDiameter);
			}
		}

		// Hand over bonds.
		if(!bonds().properties().empty()) {
			const BondsObject* existingBonds = existingParticles ? existingParticles->bonds() : nullptr;
			OORef<BondsObject> targetBonds = new BondsObject(fileSource->dataset());
			targetParticles->setBonds(targetBonds);
			targetBonds->setDataSource(fileSource);
			bonds().transferToContainer(existingBonds, targetBonds, isNewFile, cloneHelper);
		}

		// Hand over list of angles.
		if(!angles().properties().empty()) {
			const AnglesObject* existingAngles = existingParticles ? existingParticles->angles() : nullptr;
			OORef<AnglesObject> targetAngles = new AnglesObject(fileSource->dataset());
			targetParticles->setAngles(targetAngles);
			targetAngles->setDataSource(fileSource);
			angles().transferToContainer(existingAngles, targetAngles, isNewFile, cloneHelper);
		}

		// Hand over list of dihredrals.
		if(!dihedrals().properties().empty()) {
			const DihedralsObject* existingDihedrals = existingParticles ? existingParticles->dihedrals() : nullptr;
			OORef<DihedralsObject> targetDihedrals = new DihedralsObject(fileSource->dataset());
			targetParticles->setDihedrals(targetDihedrals);
			targetDihedrals->setDataSource(fileSource);
			dihedrals().transferToContainer(existingDihedrals, targetDihedrals, isNewFile, cloneHelper);
		}

		// Hand over list of impropers.
		if(!impropers().properties().empty()) {
			const ImpropersObject* existingImpropers = existingParticles ? existingParticles->impropers() : nullptr;
			OORef<ImpropersObject> targetImpropers = new ImpropersObject(fileSource->dataset());
			targetParticles->setImpropers(targetImpropers);
			targetImpropers->setDataSource(fileSource);
			impropers().transferToContainer(existingImpropers, targetImpropers, isNewFile, cloneHelper);
		}

		targetParticles->verifyIntegrity();
	}

	// Transfer voxel data.
	if(voxelGridShape() != VoxelGrid::GridDimensions{0,0,0}) {

		// Look for an existing VoxelGrid object in the old data collection.
		const VoxelGrid* existingVoxelGrid = nullptr;
		if(existing) {
			if(!voxelGridId().isEmpty()) {
				ConstDataObjectPath path = existing->getObject<VoxelGrid>(voxelGridId());
				if(!path.empty()) existingVoxelGrid = static_object_cast<VoxelGrid>(path.back());
			}
			else existingVoxelGrid = existing->getObject<VoxelGrid>();
		}

		// Create the new VoxelGrid object.
		VoxelGrid* voxelGrid = output->createObject<VoxelGrid>(voxelGridId().isEmpty() ? QStringLiteral("imported") : voxelGridId(), fileSource, voxelGridTitle());
		voxelGrid->setShape(voxelGridShape());
		voxelGrid->setDomain(output->getObject<SimulationCellObject>());

		voxels().transferToContainer(existingVoxelGrid, voxelGrid, isNewFile, cloneHelper);

		// Turn off display of voxel grid by default.
		if(isNewFile && voxelGrid->visElement())
			voxelGrid->visElement()->setEnabled(false);

		// Give the vis element of the grid a more specific title.
		if(VoxelGridVis* gridVis = voxelGrid->visElement<VoxelGridVis>())
			gridVis->setTitle(voxelGridTitle());
	}

	// Hand over timestep information and other metadata as global attributes.
	for(auto a = _attributes.cbegin(); a != _attributes.cend(); ++a) {
		output->addAttribute(a.key(), a.value(), fileSource);
	}

	// If the file parser has detected that the input file contains additional frame data following the
	// current frame, active the 'contains multiple frames' option for the importer. This will trigger
	// a scanning process for the entire file to discover all contained frames.
	if(_detectedAdditionalFrames && isNewFile) {
		if(ParticleImporter* importer = dynamic_object_cast<ParticleImporter>(fileSource->importer())) {
			importer->setMultiTimestepFile(true);
		}
	}

	return output;
}

/******************************************************************************
* Sorts the particles list with respect to particle IDs.
* Does nothing if particles do not have IDs.
******************************************************************************/
void ParticleFrameData::sortParticlesById()
{
	std::vector<size_t> invertedPermutation = particles().sortElementsById();
	if(!invertedPermutation.empty()) {

		// Update bond topology data to match new particle ordering.
		if(PropertyAccess<ParticleIndexPair> bondTopology = bonds().findStandardProperty(BondsObject::TopologyProperty)) {
			for(ParticleIndexPair& bond : bondTopology) {
				for(qlonglong& idx : bond) {
					if(idx >= 0 && idx < invertedPermutation.size())
						idx = invertedPermutation[idx];
				}
			}
		}

		// Update angle topology data to match new particle ordering.
		if(PropertyAccess<ParticleIndexTriplet> angleTopology = angles().findStandardProperty(AnglesObject::TopologyProperty)) {
			for(ParticleIndexTriplet& angle : angleTopology) {
				for(qlonglong& idx : angle) {
					if(idx >= 0 && idx < invertedPermutation.size())
						idx = invertedPermutation[idx];
				}
			}
		}

		// Update dihedral topology data to match new particle ordering.
		if(PropertyAccess<ParticleIndexQuadruplet> dihedralTopology = dihedrals().findStandardProperty(DihedralsObject::TopologyProperty)) {
			for(ParticleIndexQuadruplet& dihedral : dihedralTopology) {
				for(qlonglong& idx : dihedral) {
					if(idx >= 0 && idx < invertedPermutation.size())
						idx = invertedPermutation[idx];
				}
			}
		}

		// Update improper topology data to match new particle ordering.
		if(PropertyAccess<ParticleIndexQuadruplet> improperTopology = impropers().findStandardProperty(ImpropersObject::TopologyProperty)) {
			for(ParticleIndexQuadruplet& improper : improperTopology) {
				for(qlonglong& idx : improper) {
					if(idx >= 0 && idx < invertedPermutation.size())
						idx = invertedPermutation[idx];
				}
			}
		}
	}
}

}	// End of namespace
}	// End of namespace
