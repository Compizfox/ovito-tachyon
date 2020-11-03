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
								ConstPropertyAccess<qlonglong> correspondenceArray,
								ConstPropertyAccess<PTMAlgorithm::StructureType> structuresArray,
								Task* promise)
{
	if(promise) promise->setProgressMaximum(0);

	_neighQuery = neighQuery;
	_correspondenceArray = correspondenceArray;
	_structuresArray = structuresArray;
    structureType = PTMAlgorithm::OTHER;
	return true;
}

}	// End of namespace
}	// End of namespace


#if 0
/******************************************************************************
* Returns the number of neighbors for the PTM structure found for the current particle.
******************************************************************************/
int PTMAlgorithm::Kernel::numTemplateNeighbors() const
{
	return ptm_num_nbrs[_structureType];
}

/******************************************************************************
* Returns the neighbor information corresponding to the i-th neighbor in the
* PTM template identified for the current particle.
******************************************************************************/
const NearestNeighborFinder::Neighbor& PTMAlgorithm::Kernel::getTemplateNeighbor(int index) const
{
	OVITO_ASSERT(_structureType != OTHER);
	OVITO_ASSERT(index >= 0 && index < numTemplateNeighbors());
	int mappedIndex = _env.correspondences[index + 1] - 1;
	return getNearestNeighbor(mappedIndex);
}

/******************************************************************************
* Returns the ideal vector corresponding to the i-th neighbor in the PTM template
* identified for the current particle.
******************************************************************************/
const Vector_3<double>& PTMAlgorithm::Kernel::getIdealNeighborVector(int index) const
{
	OVITO_ASSERT(_structureType != OTHER);
	OVITO_ASSERT(index >= 0 && index < numTemplateNeighbors());
	OVITO_ASSERT(_bestTemplate != nullptr);
	return *reinterpret_cast<const Vector_3<double>*>(_bestTemplate[index + 1]);
}
#endif

