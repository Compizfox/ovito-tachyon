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
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include "ParticleExpressionEvaluator.h"

namespace Ovito { namespace Particles {

/******************************************************************************
* Initializes the list of input variables from the given input state.
******************************************************************************/
void ParticleExpressionEvaluator::createInputVariables(const std::vector<ConstPropertyPtr>& inputProperties, const SimulationCell* simCell, const QVariantMap& attributes, int animationFrame)
{
	PropertyExpressionEvaluator::createInputVariables(inputProperties, simCell, attributes, animationFrame);

	// Create computed variables for reduced particle coordinates.
	if(simCell) {
		// Look for the 'Position' particle property in the inputs.
		auto iter = boost::find_if(inputProperties, [](const ConstPropertyPtr& property) {
			return property->type() == ParticlesObject::PositionProperty;
		});
		if(iter != inputProperties.end()) {
			SimulationCell cellData = *simCell;
			ConstPropertyAccess<Point3> posProperty = *iter;
			registerComputedVariable("ReducedPosition.X", [posProperty,cellData](size_t particleIndex) -> double {
				return cellData.inverseMatrix().prodrow(posProperty[particleIndex], 0);
			});
			registerComputedVariable("ReducedPosition.Y", [posProperty,cellData](size_t particleIndex) -> double {
				return cellData.inverseMatrix().prodrow(posProperty[particleIndex], 1);
			});
			registerComputedVariable("ReducedPosition.Z", [posProperty,cellData](size_t particleIndex) -> double {
				return cellData.inverseMatrix().prodrow(posProperty[particleIndex], 2);
			});
		}
	}
}

/******************************************************************************
* Specifies the expressions to be evaluated for each bond and creates the input variables.
******************************************************************************/
void BondExpressionEvaluator::initialize(const QStringList& expressions, const PipelineFlowState& state, const ConstDataObjectPath& containerPath, int animationFrame) 
{
	PropertyExpressionEvaluator::initialize(expressions, state, containerPath, animationFrame);

	// Look for the particles object, which is the parent of the bonds object.
	if(containerPath.size() >= 2) { 
		if(const ParticlesObject* particles = dynamic_object_cast<ParticlesObject>(containerPath[containerPath.size() - 2])) {
			const BondsObject* bonds = static_object_cast<BondsObject>(containerPath.back());
			_topologyArray = bonds->getPropertyStorage(BondsObject::TopologyProperty);

			// Define computed variable 'BondLength', which yields the length of the bonds.
			if(ConstPropertyAccessAndRef<Point3> positions = particles->getProperty(ParticlesObject::PositionProperty)) {
				if(ConstPropertyAccessAndRef<ParticleIndexPair> topology = bonds->getProperty(BondsObject::TopologyProperty)) {
					ConstPropertyAccessAndRef<Vector3I> periodicImages = bonds->getProperty(BondsObject::PeriodicImageProperty);
					SimulationCell simCell;
					if(const SimulationCellObject* simCellObj = state.getObject<SimulationCellObject>())
						simCell = simCellObj->data();
					else
						periodicImages.reset();

					registerComputedVariable("BondLength", [positions=std::move(positions),topology=std::move(topology),periodicImages=std::move(periodicImages),simCell](size_t bondIndex) -> double {
						size_t index1 = topology[bondIndex][0];
						size_t index2 = topology[bondIndex][1];
						if(positions.size() > index1 && positions.size() > index2) {
							const Point3& p1 = positions[index1];
							const Point3& p2 = positions[index2];
							Vector3 delta = p2 - p1;
							if(periodicImages) {
								if(int dx = periodicImages[bondIndex][0]) delta += simCell.matrix().column(0) * (FloatType)dx;
								if(int dy = periodicImages[bondIndex][1]) delta += simCell.matrix().column(1) * (FloatType)dy;
								if(int dz = periodicImages[bondIndex][2]) delta += simCell.matrix().column(2) * (FloatType)dz;
							}
							return delta.length();
						}
						else return 0;
					},
					tr("dynamically calculated"));
				}
			}

			// Build list of particle properties that will be made available as expression variables.
			std::vector<ConstPropertyPtr> inputParticleProperties;
			for(const PropertyObject* prop : particles->properties()) {
				inputParticleProperties.push_back(prop->storage());
			}
			registerPropertyVariables(inputParticleProperties, 1, "@1.");
			registerPropertyVariables(inputParticleProperties, 2, "@2.");
		}
	}
}

/******************************************************************************
* Updates the stored value of variables that depends on the current element index.
******************************************************************************/
void BondExpressionEvaluator::updateVariables(Worker& worker, size_t elementIndex)
{
	PropertyExpressionEvaluator::updateVariables(worker, elementIndex);

	// Update values of variables referring to particle properties.
	if(_topologyArray) {
		size_t particleIndex1 = _topologyArray[elementIndex][0];
		size_t particleIndex2 = _topologyArray[elementIndex][1];
		worker.updateVariables(1, particleIndex1);
		worker.updateVariables(2, particleIndex2);
	}
}

/********************************§**********************************************
* Returns a human-readable text listing the input variables.
******************************************************************************/
QString BondExpressionEvaluator::inputVariableTable() const
{
	QString table = PropertyExpressionEvaluator::inputVariableTable();
	table.append(QStringLiteral("<p><b>Particle properties:</b><ul>"));
	table.append(QStringLiteral("<li>@1... (<i style=\"color: #555;\">property of first particle</i>)</li>"));
	table.append(QStringLiteral("<li>@2... (<i style=\"color: #555;\">property of second particle</i>)</li>"));
	table.append(QStringLiteral("</ul></p>"));
	return table;
}

}	// End of namespace
}	// End of namespace
