///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2017) Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  OVITO is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include <plugins/particles/Particles.h>
#include "ParticleBondMap.h"

namespace Ovito { namespace Particles {

#if 0
/******************************************************************************
* Reduces the size of the storage array, removing elements for which 
* the corresponding bits in the bit array are set.
******************************************************************************/
void BondsStorage::filterResize(const boost::dynamic_bitset<>& mask)
{
	OVITO_ASSERT(size() == mask.size());
	// Find the first bond to delete:
	size_t i;
	for(i = 0; i != mask.size(); i++)
		if(mask.test(i)) break;
	if(i == mask.size()) return;
	// Continue from here:
	auto b = begin() + i;
	for(; i != mask.size(); ++i) {
		if(!mask.test(i))
			*b++ = std::move((*this)[i]);
	}
	erase(b, end());
}
#endif

/******************************************************************************
* Initializes the helper class.
******************************************************************************/
ParticleBondMap::ParticleBondMap(ConstPropertyPtr bondTopology, ConstPropertyPtr bondPeriodicImages) :
	_bondTopology(std::move(bondTopology)),
	_bondPeriodicImages(std::move(bondPeriodicImages)),
	_nextBond(_bondTopology->size()*2, _bondTopology->size()*2)
{
	for(size_t bondIndex = _bondTopology->size(); bondIndex-- != 0; ) {
		size_t index1 = _bondTopology->getInt64Component(bondIndex, 0);
		size_t index2 = _bondTopology->getInt64Component(bondIndex, 1);
		if(index1 >= _startIndices.size())
			_startIndices.resize(index1 + 1, endOfListValue());
		if(index2 >= _startIndices.size())
			_startIndices.resize(index2 + 1, endOfListValue());

		size_t evenIndex = bondIndex * 2;
		size_t oddIndex  = evenIndex + 1;
		_nextBond[evenIndex] = _startIndices[index1];
		_nextBond[oddIndex]  = _startIndices[index2];
		_startIndices[index1] = evenIndex;
		_startIndices[index2] = oddIndex;
	}
}


}	// End of namespace
}	// End of namespace