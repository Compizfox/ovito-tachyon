///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2018) Alexander Stukowski
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
#include "ParticlesObject.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(ParticlesObject);
DEFINE_REFERENCE_FIELD(ParticlesObject, bonds);
SET_PROPERTY_FIELD_LABEL(ParticlesObject, bonds, "Bonds");

/******************************************************************************
* Constructor.
******************************************************************************/
ParticlesObject::ParticlesObject(DataSet* dataset) : PropertyContainer(dataset)
{
}

}	// End of namespace
}	// End of namespace
