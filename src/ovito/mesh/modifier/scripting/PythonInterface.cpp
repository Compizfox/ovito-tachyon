///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2019) Alexander Stukowski
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

#include <ovito/mesh/Mesh.h>
#include <ovito/pyscript/binding/PythonBinding.h>
#include <ovito/core/app/PluginManager.h>
#include <ovito/mesh/modifier/SurfaceMeshAffineTransformationModifierDelegate.h>
#include <ovito/mesh/modifier/SurfaceMeshReplicateModifierDelegate.h>
#include <ovito/mesh/modifier/SurfaceMeshSliceModifierDelegate.h>

namespace Ovito { namespace Mesh {

using namespace PyScript;

PYBIND11_MODULE(MeshModPython, m)
{
	// Register the classes of this plugin with the global PluginManager.
	PluginManager::instance().registerLoadedPluginClasses();

	py::options options;
	options.disable_function_signatures();

	ovito_class<SurfaceMeshAffineTransformationModifierDelegate, AffineTransformationModifierDelegate>{m};
	ovito_class<SurfaceMeshReplicateModifierDelegate, ReplicateModifierDelegate>{m};
	ovito_class<SurfaceMeshSliceModifierDelegate, SliceModifierDelegate>{m};
}

OVITO_REGISTER_PLUGIN_PYTHON_INTERFACE(MeshModPython);

}	// End of namespace
}	// End of namespace