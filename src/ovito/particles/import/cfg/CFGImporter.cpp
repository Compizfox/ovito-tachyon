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
#include <ovito/particles/import/ParticleFrameData.h>
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/particles/objects/ParticleType.h>
#include <ovito/stdobj/properties/InputColumnMapping.h>
#include <ovito/core/utilities/io/CompressedTextReader.h>
#include "CFGImporter.h"

#include <boost/algorithm/string.hpp>

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(CFGImporter);

struct CFGHeader {

	qlonglong numParticles;
	FloatType unitMultiplier;
	Matrix3 H0;
	Matrix3 transform;
	FloatType rateScale;
	bool isExtendedFormat;
	bool containsVelocities;
	QStringList auxiliaryFields;

	void parse(CompressedTextReader& stream);
};

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool CFGImporter::OOMetaClass::checkFileFormat(const FileHandle& file) const
{
	// Open input file.
	CompressedTextReader stream(file);

	// Look for the magic string 'Number of particles'.
	// It must appear within the first 20 lines of the CFG file.
	for(int i = 0; i < 20 && !stream.eof(); i++) {
		const char* line = stream.readLineTrimLeft(1024);

		// CFG files start with the string "Number of particles".
		if(boost::algorithm::starts_with(line, "Number of particles"))
			return true;

		// Terminate early if line is non-empty and contains anything else than a comment (#).
		if(line[0] > ' ' && line[0] != '#')
			return false;
	}

	return false;
}

/******************************************************************************
* Parses the header of a CFG file.
******************************************************************************/
void CFGHeader::parse(CompressedTextReader& stream)
{
	using namespace std;

	numParticles = -1;
	unitMultiplier = 1;
	H0.setIdentity();
	transform.setIdentity();
	rateScale = 1;
	isExtendedFormat = false;
	containsVelocities = true;
	int entry_count = 0;

	while(!stream.eof()) {

		string line(stream.readLineTrimLeft());

		// Ignore comments
		size_t commentChar = line.find('#');
		if(commentChar != string::npos) line.resize(commentChar);

		// Skip empty lines.
		size_t trimmedLine = line.find_first_not_of(" \t\n\r");
		if(trimmedLine == string::npos) continue;
		if(trimmedLine != 0) line = line.substr(trimmedLine);

		size_t splitChar = line.find('=');
		if(splitChar == string::npos) {
			if(stream.lineStartsWithToken(".NO_VELOCITY.")) {
				containsVelocities = false;
				continue;
			}
			break;
		}

		string key = line.substr(0, line.find_last_not_of(" \t\n\r", splitChar - 1) + 1);
		size_t valuestart = line.find_first_not_of(" \t\n\r", splitChar + 1);
		if(valuestart == string::npos) valuestart = splitChar+1;
		string value = line.substr(valuestart);

		if(key == "Number of particles") {
			if(sscanf(value.c_str(), "%lld", &numParticles) != 1 || numParticles < 0 || numParticles > 100'000'000'000ll)
				throw Exception(CFGImporter::tr("CFG file parsing error. Invalid number of atoms (line %1): %2").arg(stream.lineNumber()).arg(QString::fromStdString(value)));
		}
		else if(key == "A") unitMultiplier = atof(value.c_str());
		else if(key == "H0(1,1)") H0(0,0) = atof(value.c_str()) * unitMultiplier;
		else if(key == "H0(1,2)") H0(0,1) = atof(value.c_str()) * unitMultiplier;
		else if(key == "H0(1,3)") H0(0,2) = atof(value.c_str()) * unitMultiplier;
		else if(key == "H0(2,1)") H0(1,0) = atof(value.c_str()) * unitMultiplier;
		else if(key == "H0(2,2)") H0(1,1) = atof(value.c_str()) * unitMultiplier;
		else if(key == "H0(2,3)") H0(1,2) = atof(value.c_str()) * unitMultiplier;
		else if(key == "H0(3,1)") H0(2,0) = atof(value.c_str()) * unitMultiplier;
		else if(key == "H0(3,2)") H0(2,1) = atof(value.c_str()) * unitMultiplier;
		else if(key == "H0(3,3)") H0(2,2) = atof(value.c_str()) * unitMultiplier;
		else if(key == "Transform(1,1)") transform(0,0) = atof(value.c_str());
		else if(key == "Transform(1,2)") transform(0,1) = atof(value.c_str());
		else if(key == "Transform(1,3)") transform(0,2) = atof(value.c_str());
		else if(key == "Transform(2,1)") transform(1,0) = atof(value.c_str());
		else if(key == "Transform(2,2)") transform(1,1) = atof(value.c_str());
		else if(key == "Transform(2,3)") transform(1,2) = atof(value.c_str());
		else if(key == "Transform(3,1)") transform(2,0) = atof(value.c_str());
		else if(key == "Transform(3,2)") transform(2,1) = atof(value.c_str());
		else if(key == "Transform(3,3)") transform(2,2) = atof(value.c_str());
		else if(key == "eta(1,1)") {}
		else if(key == "eta(1,2)") {}
		else if(key == "eta(1,3)") {}
		else if(key == "eta(2,2)") {}
		else if(key == "eta(2,3)") {}
		else if(key == "eta(3,3)") {}
		else if(key == "R") rateScale = atof(value.c_str());
		else if(key == "entry_count") {
			entry_count = atoi(value.c_str());
			isExtendedFormat = true;
		}
		else if(key.compare(0, 10, "auxiliary[") == 0) {
			isExtendedFormat = true;
			size_t endOfName = value.find_first_of(" \t");
			auxiliaryFields.push_back(QString::fromStdString(value.substr(0, endOfName)).trimmed());
		}
		else {
			throw Exception(CFGImporter::tr("Unknown key in CFG file header at line %1: %2").arg(stream.lineNumber()).arg(QString::fromStdString(line)));
		}
	}
	if(numParticles < 0)
		throw Exception(CFGImporter::tr("Invalid file header. This is not a valid CFG file."));
}

/******************************************************************************
* Parses the given input file.
******************************************************************************/
FileSourceImporter::FrameDataPtr CFGImporter::FrameLoader::loadFile()
{
	// Open file for reading.
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Reading CFG file %1").arg(fileHandle().toString()));

	// Jump to byte offset.
	if(frame().byteOffset != 0)
		stream.seek(frame().byteOffset, frame().lineNumber);

	// Parse header of CFG file.
	CFGHeader header;
	header.parse(stream);

	// Create the destination container for loaded data.
	auto frameData = std::make_shared<ParticleFrameData>();

	ParticleInputColumnMapping cfgMapping;
	if(header.isExtendedFormat == false) {
		cfgMapping.resize(8);
		cfgMapping.mapStandardColumn(0, ParticlesObject::MassProperty);
		cfgMapping.mapStandardColumn(1, ParticlesObject::TypeProperty);
		cfgMapping.mapStandardColumn(2, ParticlesObject::PositionProperty, 0);
		cfgMapping.mapStandardColumn(3, ParticlesObject::PositionProperty, 1);
		cfgMapping.mapStandardColumn(4, ParticlesObject::PositionProperty, 2);
		cfgMapping.mapStandardColumn(5, ParticlesObject::VelocityProperty, 0);
		cfgMapping.mapStandardColumn(6, ParticlesObject::VelocityProperty, 1);
		cfgMapping.mapStandardColumn(7, ParticlesObject::VelocityProperty, 2);
	}
	else {
		cfgMapping.resize(header.containsVelocities ? 6 : 3);
		cfgMapping.mapStandardColumn(0, ParticlesObject::PositionProperty, 0);
		cfgMapping.mapStandardColumn(1, ParticlesObject::PositionProperty, 1);
		cfgMapping.mapStandardColumn(2, ParticlesObject::PositionProperty, 2);
		if(header.containsVelocities) {
			cfgMapping.mapStandardColumn(3, ParticlesObject::VelocityProperty, 0);
			cfgMapping.mapStandardColumn(4, ParticlesObject::VelocityProperty, 1);
			cfgMapping.mapStandardColumn(5, ParticlesObject::VelocityProperty, 2);
		}
		generateAutomaticColumnMapping(cfgMapping, header.auxiliaryFields);
	}

	setProgressMaximum(header.numParticles);

	// Prepare the mapping between input file columns and particle properties.
	InputColumnReader columnParser(cfgMapping, frameData->particles(), header.numParticles);

	// Create particle mass and type properties.
	int currentAtomType = 0;
	FloatType currentMass = 0;
	PropertyAccess<int> typeProperty;
	PropertyAccess<FloatType> massProperty;
	PropertyContainerImportData::TypeList* typeList = nullptr;
	if(header.isExtendedFormat) {
		typeProperty = frameData->particles().createStandardProperty<ParticlesObject>(header.numParticles, ParticlesObject::TypeProperty, false);
		typeList = frameData->particles().createPropertyTypesList(typeProperty, ParticleType::OOClass());
		massProperty = frameData->particles().createStandardProperty<ParticlesObject>(header.numParticles, ParticlesObject::MassProperty, false);
	}

	// Read per-particle data.
	bool isFirstLine = true;
	for(qlonglong particleIndex = 0; particleIndex < header.numParticles; ) {

		// Update progress indicator.
		if(!setProgressValueIntermittent(particleIndex))
			return {};

		if(!isFirstLine)
			stream.readLine();
		else
			isFirstLine = false;

		if(header.isExtendedFormat) {
			bool isNewType = true;
			const char* line = stream.line();
			while(*line != '\0' && *line <= ' ') ++line;
			for(; *line != '\0'; ++line) {
				if(*line <= ' ') {
					for(; *line != '\0'; ++line) {
						if(*line > ' ') {
							isNewType = false;
							break;
						}
					}
					break;
				}
			}
			if(isNewType) {
				// Parse mass and atom type name.
				currentMass = atof(stream.line());
				const char* line = stream.readLineTrimLeft();
				const char* line_end = line;
				while(*line_end != '\0' && *line_end > ' ') ++line_end;
				currentAtomType = typeList->addTypeName(line, line_end);

				continue;
			}

			typeProperty[particleIndex] = currentAtomType;
			massProperty[particleIndex] = currentMass;
		}

		try {
			columnParser.readElement(particleIndex, stream.line());
			particleIndex++;
		}
		catch(Exception& ex) {
			throw ex.prependGeneralMessage(tr("Parsing error in line %1 of CFG file.").arg(stream.lineNumber()));
		}
	}

	// Since we created particle types on the go while reading the particles, the assigned particle type IDs
	// depend on the storage order of particles in the file. We rather want a well-defined particle type ordering, that's
	// why we sort them now.
	columnParser.sortElementTypes();
	if(header.isExtendedFormat)
		typeList->sortTypesByName(typeProperty);

	AffineTransformation H((header.transform * header.H0).transposed());
	H.translation() = H * Vector3(-0.5, -0.5, -0.5);
	frameData->simulationCell().setMatrix(H);

	// The CFG file stores particle positions in reduced coordinates.
	// Rescale them now to absolute (Cartesian) coordinates.
	// However, do this only if no absolute coordinates have been read from the extra data columns in the CFG file.
	if(PropertyAccess<Point3> posProperty = frameData->particles().findStandardProperty(ParticlesObject::PositionProperty)) {
		for(Point3& p : posProperty)
			p = H * p;
	}

	// Sort particles by ID if requested.
	if(_sortParticles)
		frameData->sortParticlesById();

	frameData->setStatus(tr("Number of particles: %1").arg(header.numParticles));
	return frameData;
}

/******************************************************************************
 * Guesses the mapping of input file columns to internal particle properties.
 *****************************************************************************/
void CFGImporter::generateAutomaticColumnMapping(ParticleInputColumnMapping& columnMapping, const QStringList& columnNames)
{
	int i = columnMapping.size();
	columnMapping.resize(columnMapping.size() + columnNames.size());
	for(int j = 0; j < columnNames.size(); j++, i++) {
		columnMapping[i].columnName = columnNames[j];
		QString name = columnNames[j].toLower();
		if(name == "vx" || name == "velocities") columnMapping.mapStandardColumn(i, ParticlesObject::VelocityProperty, 0);
		else if(name == "vy") columnMapping.mapStandardColumn(i, ParticlesObject::VelocityProperty, 1);
		else if(name == "vz") columnMapping.mapStandardColumn(i, ParticlesObject::VelocityProperty, 2);
		else if(name == "v") columnMapping.mapStandardColumn(i, ParticlesObject::VelocityMagnitudeProperty);
		else if(name == "id") columnMapping.mapStandardColumn(i, ParticlesObject::IdentifierProperty);
		else if(name == "radius") columnMapping.mapStandardColumn(i, ParticlesObject::RadiusProperty);
		else if(name == "q") columnMapping.mapStandardColumn(i, ParticlesObject::ChargeProperty);
		else if(name == "ix") columnMapping.mapStandardColumn(i, ParticlesObject::PeriodicImageProperty, 0);
		else if(name == "iy") columnMapping.mapStandardColumn(i, ParticlesObject::PeriodicImageProperty, 1);
		else if(name == "iz") columnMapping.mapStandardColumn(i, ParticlesObject::PeriodicImageProperty, 2);
		else if(name == "fx") columnMapping.mapStandardColumn(i, ParticlesObject::ForceProperty, 0);
		else if(name == "fy") columnMapping.mapStandardColumn(i, ParticlesObject::ForceProperty, 1);
		else if(name == "fz") columnMapping.mapStandardColumn(i, ParticlesObject::ForceProperty, 2);
		else if(name == "mux") columnMapping.mapStandardColumn(i, ParticlesObject::DipoleOrientationProperty, 0);
		else if(name == "muy") columnMapping.mapStandardColumn(i, ParticlesObject::DipoleOrientationProperty, 1);
		else if(name == "muz") columnMapping.mapStandardColumn(i, ParticlesObject::DipoleOrientationProperty, 2);
		else if(name == "mu") columnMapping.mapStandardColumn(i, ParticlesObject::DipoleMagnitudeProperty);
		else if(name == "omegax") columnMapping.mapStandardColumn(i, ParticlesObject::AngularVelocityProperty, 0);
		else if(name == "omegay") columnMapping.mapStandardColumn(i, ParticlesObject::AngularVelocityProperty, 1);
		else if(name == "omegaz") columnMapping.mapStandardColumn(i, ParticlesObject::AngularVelocityProperty, 2);
		else if(name == "angmomx") columnMapping.mapStandardColumn(i, ParticlesObject::AngularMomentumProperty, 0);
		else if(name == "angmomy") columnMapping.mapStandardColumn(i, ParticlesObject::AngularMomentumProperty, 1);
		else if(name == "angmomz") columnMapping.mapStandardColumn(i, ParticlesObject::AngularMomentumProperty, 2);
		else if(name == "tqx") columnMapping.mapStandardColumn(i, ParticlesObject::TorqueProperty, 0);
		else if(name == "tqy") columnMapping.mapStandardColumn(i, ParticlesObject::TorqueProperty, 1);
		else if(name == "tqz") columnMapping.mapStandardColumn(i, ParticlesObject::TorqueProperty, 2);
		else if(name == "spin") columnMapping.mapStandardColumn(i, ParticlesObject::SpinProperty);
		else {
			columnMapping.mapCustomColumn(i, columnNames[j], PropertyStorage::Float);
		}
	}
}

}	// End of namespace
}	// End of namespace
