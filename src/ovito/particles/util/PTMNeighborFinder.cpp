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

#include <ovito/particles/Particles.h>
#include <ovito/core/utilities/concurrent/Task.h>
#include <ovito/particles/modifier/analysis/ptm/PTMAlgorithm.h>
#include "PTMNeighborFinder.h"

namespace Ovito { namespace Particles {

/******************************************************************************
* Prepares the neighbor list builder.
******************************************************************************/
bool PTMNeighborFinder::prepare(NearestNeighborFinder::Query<PTMAlgorithm::MAX_INPUT_NEIGHBORS> *neighQuery,
								ConstPropertyAccess<PTMAlgorithm::StructureType> structuresArray,
								ConstPropertyAccess<Quaternion> orientationsArray,
								ConstPropertyAccess<qlonglong> correspondencesArray,
								Task* promise)
{
	if(promise) promise->setProgressMaximum(0);
	OVITO_ASSERT(neighQuery);
	OVITO_ASSERT(structuresArray);
	OVITO_ASSERT(orientationsArray);
	OVITO_ASSERT(correspondencesArray);

	_neighQuery = neighQuery;
	_structuresArray = structuresArray;
	_orientationsArray = orientationsArray;
	_correspondencesArray = correspondencesArray;
	structureType = PTMAlgorithm::OTHER;
	rmsd = std::numeric_limits<FloatType>::infinity();
	_numNeighbors = 0;
	_templateIndex = 0;
	return true;
}

}	// End of namespace
}	// End of namespace

