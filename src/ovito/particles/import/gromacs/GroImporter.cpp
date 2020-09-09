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
#include <ovito/core/utilities/io/CompressedTextReader.h>
#include <ovito/core/utilities/io/NumberParsing.h>
#include "GroImporter.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(GroImporter);

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool GroImporter::OOMetaClass::checkFileFormat(const FileHandle& file) const
{
	// Open input file.
	CompressedTextReader stream(file);

	// Skip first comment line.
	stream.readLine(1024);

	// Read second line.
	const char* p = stream.readLineTrimLeft(128);
	if(*p == '\0')
		return false;

	// Expect a digit (number of atoms).
	if(*p < '0' || *p > '9') return false;
	// Skip digits.
	while(*p >= '0' && *p <= '9') {
		if(*p == '\0') break;
		++p;
	}
	// The number of atoms should be followed by nothing but whitespace.
	bool foundNewline = false;
	while(*p != '\0') {
		if(*p > ' ') return false;
		if(*p == '\n' || *p == '\r')
			foundNewline = true;
		++p;
	}

	if(!foundNewline) return false;

	// Read third line.
	int i1, i2;
	char s1[6], s2[6];
	if(sscanf(stream.readLine(), "%5i%5s%5s%5i", &i1, s1, s2, &i2) != 4 || i1 < 1 || i2 < 1) {
		return false;
	}

	return true;
}

/******************************************************************************
* Scans the data file and builds a list of source frames.
******************************************************************************/
void GroImporter::FrameFinder::discoverFramesInFile(QVector<FileSourceImporter::Frame>& frames)
{
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Scanning file %1").arg(fileHandle().toString()));
	setProgressMaximum(stream.underlyingSize());

	int frameNumber = 0;
	QString filename = fileHandle().sourceUrl().fileName();
	Frame frame(fileHandle());

	while(!stream.eof() && !isCanceled()) {
		frame.byteOffset = stream.byteOffset();
		frame.lineNumber = stream.lineNumber();

		// Skip comment line.
		stream.readLine();

		// Parse number of atoms.
		const char* line = stream.readLineTrimLeft();

		if(line[0] == '\0') break;

		unsigned long long numParticlesLong;
		int charCount;
		if(sscanf(line, "%llu%n", &numParticlesLong, &charCount) != 1)
			throw Exception(tr("Invalid number of atoms in line %1 of Gromacs file: %2").arg(stream.lineNumber()).arg(stream.lineString().trimmed()));

		// Check trailing whitespace. There should be nothing else but the number of atoms on the second line.
		for(const char* p = line + charCount; *p != '\0'; ++p) {
			if(*p > ' ')
				throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected token following number of atoms:\n\n\"%2\"").arg(stream.lineNumber()).arg(stream.lineString().trimmed()));
		}

		// Create a new record for the time step.
		frame.label = QStringLiteral("%1 (Frame %2)").arg(filename).arg(frameNumber++);
		frames.push_back(frame);

		// Skip atom lines.
		for(unsigned long long i = 0; i < numParticlesLong; i++) {
			stream.readLine();
			if(!setProgressValueIntermittent(stream.underlyingByteOffset()))
				return;
		}

		// Skip cell geometry line.
		stream.readLine();
	}
}

/******************************************************************************
* Parses the given input file.
******************************************************************************/
FileSourceImporter::FrameDataPtr GroImporter::FrameLoader::loadFile()
{
	// Open file for reading.
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Reading Gromacs file %1").arg(fileHandle().toString()));

	// Jump to byte offset.
	if(frame().byteOffset != 0)
		stream.seek(frame().byteOffset, frame().lineNumber);

	// Create the destination container for loaded data.
	auto frameData = std::make_shared<ParticleFrameData>();

	// Read comment line.
	stream.readLine();
	QString commentLine = stream.lineString().trimmed();

	// Parse number of atoms.
	unsigned long long numParticles;
	int charCount;
	if(sscanf(stream.readLine(), "%llu%n", &numParticles, &charCount) != 1)
		throw Exception(tr("Invalid number of particles in line %1 of Gromacs file: %2").arg(stream.lineNumber()).arg(stream.lineString().trimmed()));

	// Check trailing whitespace. There should be nothing else but the number of atoms on the first line.
	for(const char* p = stream.line() + charCount; *p != '\0'; ++p) {
		if(*p > ' ')
			throw Exception(tr("Parsing error in line %1 of Gromacs file. The second line of a .gro file should contain just the number of atoms and nothing else.").arg(stream.lineNumber()));
	}
	if(numParticles > (unsigned long long)std::numeric_limits<int>::max())
		throw Exception(tr("Too many atoms in Gromacs file. This program version can read files with up to %1 atoms only.").arg(std::numeric_limits<int>::max()));
	setProgressMaximum(numParticles);

	// Create particle properties.
	PropertyAccess<Point3> posProperty = frameData->particles().createStandardProperty<ParticlesObject>(numParticles, ParticlesObject::PositionProperty, true);
	PropertyAccess<int> residueTypeProperty = frameData->particles().addProperty(std::make_shared<PropertyStorage>(numParticles, PropertyStorage::Int, 1, 0, QStringLiteral("Residue Type"), false));
	PropertyAccess<qlonglong> residueNumberProperty = frameData->particles().addProperty(std::make_shared<PropertyStorage>(numParticles, PropertyStorage::Int64, 1, 0, QStringLiteral("Residue Identifier"), false));
	PropertyAccess<int> typeProperty = frameData->particles().createStandardProperty<ParticlesObject>(numParticles, ParticlesObject::TypeProperty, false);
	PropertyAccess<qlonglong> identifierProperty = frameData->particles().createStandardProperty<ParticlesObject>(numParticles, ParticlesObject::IdentifierProperty, true);
	PropertyAccess<Vector3> velocityProperty;

	// Create particle and residue type lists, because we need to populate them while parsing.
	std::unique_ptr<PropertyContainerImportData::TypeList> typeList = std::make_unique<PropertyContainerImportData::TypeList>(ParticleType::OOClass());
	std::unique_ptr<PropertyContainerImportData::TypeList> residueTypeList = std::make_unique<PropertyContainerImportData::TypeList>(ElementType::OOClass());

	// Parse list of atoms.
	for(size_t i = 0; i < numParticles; i++) {
		if(!setProgressValueIntermittent(i)) return {};
		const char* token = stream.readLine();

		// Parse residue number (5 characters).
		const char* token_end = token + 5;
		while(token < token_end && *token <= ' ') {
			if(*token == '\0')
				throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
			++token;
		}
		int residueNumber;
		bool ok = parseInt(token, token_end, residueNumber);
		if(!ok)
			throw Exception(tr("Parsing error in line %1 of Gromacs file. Invalid residue number.").arg(stream.lineNumber()));
		token = token_end;

		// Parse residue name (5 characters).
		char residueName[6];
		char* residueNameStart = residueName;
		token_end = token + 5;
		for(; token != token_end; token++, residueNameStart++) {
			if(*token == '\0')
				throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
			if(*token > ' ')
				break;
		}
		char* residueNameEnd = residueNameStart;
		for(; token != token_end; token++, residueNameEnd++) {
			if(*token == '\0')
				throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
			if(*token <= ' ')
				break;
			*residueNameEnd = *token;
		}
		token = token_end;

		// Parse atom name (5 characters).
		char atomName[6];
		char* atomNameStart = atomName;
		token_end = token + 5;
		for(; token != token_end; token++, atomNameStart++) {
			if(*token == '\0')
				throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
			if(*token > ' ')
				break;
		}
		char* atomNameEnd = atomNameStart;
		for(; token != token_end; token++, atomNameEnd++) {
			if(*token == '\0')
				throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
			if(*token <= ' ')
				break;
			*atomNameEnd = *token;
		}
		token = token_end;

		// Parse atom number (5 characters).
		token_end = token + 5;
		while(token < token_end && *token <= ' ') {
			if(*token == '\0')
				throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
			++token;
		}
		int atomNumber;
		ok = parseInt(token, token_end, atomNumber);
		if(!ok || atomNumber <= 0 || atomNumber > numParticles)
			throw Exception(tr("Parsing error in line %1 of Gromacs file. Invalid atom number.").arg(stream.lineNumber()));
		token = token_end;
		size_t atomIndex = atomNumber - 1;

		// Store parsed value in property arrays.
		identifierProperty.set(atomIndex, atomNumber);
		typeProperty.set(atomIndex, typeList->addTypeName(atomNameStart, atomNameEnd));
		residueTypeProperty.set(atomIndex, residueTypeList->addTypeName(residueNameStart, residueNameEnd));
		residueNumberProperty.set(atomIndex, residueNumber);

		// Parse atomic xyz coordinates.
		// First, determine column width by counting distance between decimal points.
		const char* c = token;
		while(*c != '\0' && *c != '.')
			c++;
		if(*c == '\0')
			throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
		int columnWidth = 1;
		for(const char* c2 = c + 1; *c2 != '\0' && *c2 != '.'; ++c2)
			columnWidth++;
		Point3& pos = posProperty[atomIndex];
		for(size_t dim = 0; dim < 3; dim++) {
			token_end = token + columnWidth;
			while(token < token_end && *token <= ' ') {
				if(*token == '\0')
					throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
				++token;
			}
			ok = parseFloatType(token, token_end, pos[dim]);
			if(!ok)
				throw Exception(tr("Parsing error in line %1 of Gromacs file. Invalid atomic coordinate (col width=%2).").arg(stream.lineNumber()).arg(columnWidth));
			token = token_end;
		}

		// Convert coordinates from nanometers to angstroms.
		pos *= FloatType(10);

		// Parse atomic velocity vectors (optional). Gromacs files use velocity units nm/ps (or km/s).
		// First, determine column width by counting distance between decimal points.
		c = token;
		while(*c != '\0' && *c != '.')
			c++;
		if(*c != '\0') {
			columnWidth = 1;
			for(const char* c2 = c + 1; *c2 != '\0' && *c2 != '.'; ++c2)
				columnWidth++;
			if(!velocityProperty)
				velocityProperty = frameData->particles().createStandardProperty<ParticlesObject>(numParticles, ParticlesObject::VelocityProperty, false);
			Vector3& v = velocityProperty[atomIndex];
			for(size_t dim = 0; dim < 3; dim++) {
				token_end = token + columnWidth;
				while(token < token_end && *token <= ' ') {
					if(*token == '\0')
						throw Exception(tr("Parsing error in line %1 of Gromacs file. Unexpected end of line.").arg(stream.lineNumber()));
					++token;
				}
				ok = parseFloatType(token, token_end, v[dim]);
				if(!ok)
					throw Exception(tr("Parsing error in line %1 of Gromacs file. Invalid atomic velocity vector (col width=%2).").arg(stream.lineNumber()).arg(columnWidth));
				token = token_end;
			}
		}
	}

	// Since we created particle types on the go while reading the particles, the assigned particle type IDs
	// depend on the storage order of particles in the file. We rather want a well-defined particle type ordering, that's
	// why we sort them now.
	typeList->sortTypesByName(typeProperty);
	residueTypeList->sortTypesByName(residueTypeProperty);
	frameData->particles().setPropertyTypesList(typeProperty, std::move(typeList));
	frameData->particles().setPropertyTypesList(residueTypeProperty, std::move(residueTypeList));

	// Parse simulation cell definition.
	AffineTransformation cell = AffineTransformation::Identity();
	if(sscanf(stream.readLine(), 
		FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " "
		FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " "
		FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING,
		&cell(0,0), &cell(1,1), &cell(2,2),
		&cell(1,0), &cell(2,0), &cell(0,1),
		&cell(2,1), &cell(0,2), &cell(1,2)
	  ) < 3)
		throw Exception(tr("Parsing error in line %1 of Gromacs file. Invalid simulation cell definition: %2").arg(stream.lineNumber()).arg(stream.lineString()));
	// Convert cell size from nanometers to angstroms.
	frameData->simulationCell().setMatrix(cell * FloatType(10));

	// Detect if there are more simulation frames following in the file.
	if(!stream.eof())
		frameData->signalAdditionalFrames();

	if(commentLine.isEmpty())
		frameData->setStatus(tr("%1 atoms").arg(numParticles));
	else
		frameData->setStatus(tr("%1 atoms\n%2").arg(numParticles).arg(commentLine));

	return frameData;
}

}	// End of namespace
}	// End of namespace
