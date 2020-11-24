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
#include <ovito/particles/modifier/modify/LoadTrajectoryModifier.h>
#include <ovito/core/dataset/scene/PipelineSceneNode.h>
#include <ovito/core/app/Application.h>
#include <ovito/core/dataset/io/FileSource.h>
#include "ParticleImporter.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(ParticleImporter);
DEFINE_PROPERTY_FIELD(ParticleImporter, sortParticles);
SET_PROPERTY_FIELD_LABEL(ParticleImporter, sortParticles, "Sort particles by ID");

/******************************************************************************
* Is called when the value of a property of this object has changed.
******************************************************************************/
void ParticleImporter::propertyChanged(const PropertyFieldDescriptor& field)
{
	FileSourceImporter::propertyChanged(field);

	if(field == PROPERTY_FIELD(sortParticles)) {
		// Reload input file(s) when this option has been changed.
		// But no need to refetch the files from the remote location. Reparsing the cached files is sufficient.
		requestReload();
	}
}

/******************************************************************************
* Is called when importing multiple files of different formats.
******************************************************************************/
bool ParticleImporter::importFurtherFiles(std::vector<std::pair<QUrl, OORef<FileImporter>>> sourceUrlsAndImporters, ImportMode importMode, bool autodetectFileSequences, PipelineSceneNode* pipeline)
{
	OVITO_ASSERT(!sourceUrlsAndImporters.empty());
	OORef<ParticleImporter> nextImporter = dynamic_object_cast<ParticleImporter>(sourceUrlsAndImporters.front().second);
	if(this->isTrajectoryFormat() == false && nextImporter && nextImporter->isTrajectoryFormat() == true) {

		// Create a new file source for loading the trajectory.
		OORef<FileSource> fileSource = new FileSource(dataset());
		// Load user-defined default settings.
		if(Application::instance()->executionContext() == Application::ExecutionContext::Interactive)
			fileSource->loadUserDefaults();

		// Concatenate all files from the input list having the same file format into one sequence,
		// which gets handled by the trajectory importer.
		std::vector<QUrl> sourceUrls;
		sourceUrls.push_back(std::move(sourceUrlsAndImporters.front().first));
		auto iter = std::next(sourceUrlsAndImporters.begin());
		for(; iter != sourceUrlsAndImporters.end(); ++iter) {
			if(iter->second->getOOClass() != nextImporter->getOOClass())
				break;
			sourceUrls.push_back(std::move(iter->first));		
		}
		sourceUrlsAndImporters.erase(sourceUrlsAndImporters.begin(), iter);

		// Set the input file location(s) and importer.
		if(!fileSource->setSource(std::move(sourceUrls), nextImporter, autodetectFileSequences))
			return {};

		// Create a modifier for injecting the trajectory data into the existing pipeline.
		OORef<LoadTrajectoryModifier> loadTrjMod = new LoadTrajectoryModifier(dataset());
		loadTrjMod->setTrajectorySource(std::move(fileSource));
		pipeline->applyModifier(loadTrjMod);

		if(sourceUrlsAndImporters.empty())
			return true;
	}
	return FileSourceImporter::importFurtherFiles(std::move(sourceUrlsAndImporters), importMode, autodetectFileSequences, pipeline);
}

}	// End of namespace
}	// End of namespace
