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
#include <ovito/core/app/Application.h>
#include <ovito/core/utilities/io/CompressedTextReader.h>
#include <ovito/core/utilities/io/FileManager.h>
#include "LAMMPSDumpLocalImporter.h"

#include <QRegularExpression>

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(LAMMPSDumpLocalImporter);
DEFINE_PROPERTY_FIELD(LAMMPSDumpLocalImporter, columnMapping);
SET_PROPERTY_FIELD_LABEL(LAMMPSDumpLocalImporter, columnMapping, "File column mapping");

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool LAMMPSDumpLocalImporter::OOMetaClass::checkFileFormat(const FileHandle& file) const
{
	// Open input file.
	CompressedTextReader stream(file);

	// Read first line.
	stream.readLine(15);

	// Dump files written by LAMMPS start with one of the following keywords: TIMESTEP, UNITS or TIME.  
	if(!stream.lineStartsWith("ITEM: TIMESTEP") && !stream.lineStartsWith("ITEM: UNITS") && !stream.lineStartsWith("ITEM: TIME"))
		return false;

	// Continue reading until "ITEM: NUMBER OF ENTRIES" line is encountered.
	for(int i = 0; i < 20; i++) {
		if(stream.eof())
			return false;
		stream.readLine();
		if(stream.lineStartsWith("ITEM: NUMBER OF ENTRIES"))
			return true;
	}

	return false;
}

/******************************************************************************
* Inspects the header of the given file and returns the number of file columns.
******************************************************************************/
Future<InputColumnMapping> LAMMPSDumpLocalImporter::inspectFileHeader(const Frame& frame)
{
	// Retrieve file.
	return Application::instance()->fileManager()->fetchUrl(dataset()->taskManager(), frame.sourceFile)
		.then(executor(), [this, frame](const FileHandle& fileHandle) {

			// Start task that inspects the file header to determine the contained data columns.
			activateCLocale();
			FrameLoaderPtr inspectionTask = std::make_shared<FrameLoader>(frame, fileHandle);
			return dataset()->taskManager().runTaskAsync(inspectionTask)
				.then([](const FileSourceImporter::FrameDataPtr& frameData) {
					return static_cast<LAMMPSFrameData*>(frameData.get())->detectedColumnMapping();
				});
		});
}

/******************************************************************************
* Scans the data file and builds a list of source frames.
******************************************************************************/
void LAMMPSDumpLocalImporter::FrameFinder::discoverFramesInFile(QVector<FileSourceImporter::Frame>& frames)
{
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Scanning LAMMPS dump local file %1").arg(fileHandle().toString()));
	setProgressMaximum(stream.underlyingSize());

	// Regular expression for whitespace characters.
	QRegularExpression ws_re(QStringLiteral("\\s+"));

	unsigned long long timestep = 0;
	size_t numElements = 0;
	Frame frame(fileHandle());

	while(!stream.eof() && !isCanceled()) {
		qint64 byteOffset = stream.byteOffset();
		int lineNumber = stream.lineNumber();

		// Parse next line.
		stream.readLine();

		do {
			if(stream.lineStartsWith("ITEM: TIMESTEP")) {
				if(sscanf(stream.readLine(), "%llu", &timestep) != 1)
					throw Exception(tr("LAMMPS dump local file parsing error. Invalid timestep number (line %1):\n%2").arg(stream.lineNumber()).arg(stream.lineString()));
				frame.byteOffset = byteOffset;
				frame.lineNumber = lineNumber;
				frame.label = QString("Timestep %1").arg(timestep);
				frames.push_back(frame);
				break;
			}
			else if(stream.lineStartsWithToken("ITEM: TIME")) {
				stream.readLine();
				stream.readLine();
			}
			else if(stream.lineStartsWith("ITEM: NUMBER OF ENTRIES")) {
				// Parse number of entries.
				unsigned long long u;
				if(sscanf(stream.readLine(), "%llu", &u) != 1)
					throw Exception(tr("LAMMPS dump local file parsing error. Invalid number of entries in line %1:\n%2").arg(stream.lineNumber()).arg(stream.lineString()));
				if(u > 100'000'000'000ll)
					throw Exception(tr("LAMMPS dump local file parsing error. Number of entries in line %1 is too large. The LAMMPS dump local file reader doesn't accept files with more than 100 entries.").arg(stream.lineNumber()));
				numElements = (size_t)u;
				break;
			}
			else if(stream.lineStartsWith("ITEM: ENTRIES")) {
				for(size_t i = 0; i < numElements; i++) {
					stream.readLine();
					if(!setProgressValueIntermittent(stream.underlyingByteOffset()))
						return;
				}
				break;
			}
			else if(stream.lineStartsWith("ITEM:")) {
				// Skip lines up to next ITEM:
				while(!stream.eof()) {
					byteOffset = stream.byteOffset();
					stream.readLine();
					if(stream.lineStartsWith("ITEM:"))
						break;
				}
			}
			else {
				throw Exception(tr("LAMMPS dump local file parsing error. Line %1 of file %2 is invalid.").arg(stream.lineNumber()).arg(stream.filename()));
			}
		}
		while(!stream.eof());
	}
}

/******************************************************************************
* Parses the given input file.
******************************************************************************/
FileSourceImporter::FrameDataPtr LAMMPSDumpLocalImporter::FrameLoader::loadFile()
{
	// Open file for reading.
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Reading LAMMPS dump local file %1").arg(fileHandle().toString()));

	// Jump to byte offset.
	if(frame().byteOffset != 0)
		stream.seek(frame().byteOffset, frame().lineNumber);

	// Create the destination container for loaded data.
	auto frameData = std::make_shared<LAMMPSFrameData>();

	// Hide particles, because this importer loads non-particle data.
	frameData->particles().setVisElementClass(nullptr);

	// Regular expression for whitespace characters.
	QRegularExpression ws_re(QStringLiteral("\\s+"));

	unsigned long long timestep;
	size_t numElements = 0;

	while(!stream.eof()) {

		// Parse next line.
		stream.readLine();

		do {
			if(stream.lineStartsWith("ITEM: TIMESTEP")) {
				if(sscanf(stream.readLine(), "%llu", &timestep) != 1)
					throw Exception(tr("LAMMPS dump local file parsing error. Invalid timestep number (line %1):\n%2").arg(stream.lineNumber()).arg(stream.lineString()));
				frameData->attributes().insert(QStringLiteral("Timestep"), QVariant::fromValue(timestep));
				break;
			}
			else if(stream.lineStartsWithToken("ITEM: TIME")) {
				FloatType simulationTime;
				if(sscanf(stream.readLine(), FLOATTYPE_SCANF_STRING, &simulationTime) != 1)
					throw Exception(tr("LAMMPS dump local file parsing error. Invalid time value (line %1):\n%2").arg(stream.lineNumber()).arg(stream.lineString()));
				frameData->attributes().insert(QStringLiteral("Time"), QVariant::fromValue(simulationTime));
				break;
			}
			else if(stream.lineStartsWith("ITEM: NUMBER OF ENTRIES")) {
				// Parse number of entries.
				unsigned long long u;
				if(sscanf(stream.readLine(), "%llu", &u) != 1)
					throw Exception(tr("LAMMPS dump local file parsing error. Invalid number of entries in line %1:\n%2").arg(stream.lineNumber()).arg(stream.lineString()));

				numElements = (size_t)u;
				setProgressMaximum(u);
				break;
			}
			else if(stream.lineStartsWith("ITEM: BOX BOUNDS xy xz yz")) {

				// Parse optional boundary condition flags.
				QStringList tokens = stream.lineString().mid(qstrlen("ITEM: BOX BOUNDS xy xz yz")).split(ws_re, QString::SkipEmptyParts);
				if(tokens.size() >= 3)
					frameData->simulationCell().setPbcFlags(tokens[0] == "pp", tokens[1] == "pp", tokens[2] == "pp");

				// Parse triclinic simulation box.
				FloatType tiltFactors[3];
				Box3 simBox;
				for(int k = 0; k < 3; k++) {
					if(sscanf(stream.readLine(), FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING, &simBox.minc[k], &simBox.maxc[k], &tiltFactors[k]) != 3)
						throw Exception(tr("Invalid box size in line %1 of LAMMPS dump local file: %2").arg(stream.lineNumber()).arg(stream.lineString()));
				}

				// LAMMPS only stores the outer bounding box of the simulation cell in the dump file.
				// We have to determine the size of the actual triclinic cell.
				simBox.minc.x() -= std::min(std::min(std::min(tiltFactors[0], tiltFactors[1]), tiltFactors[0]+tiltFactors[1]), (FloatType)0);
				simBox.maxc.x() -= std::max(std::max(std::max(tiltFactors[0], tiltFactors[1]), tiltFactors[0]+tiltFactors[1]), (FloatType)0);
				simBox.minc.y() -= std::min(tiltFactors[2], (FloatType)0);
				simBox.maxc.y() -= std::max(tiltFactors[2], (FloatType)0);
				frameData->simulationCell().setMatrix(AffineTransformation(
						Vector3(simBox.sizeX(), 0, 0),
						Vector3(tiltFactors[0], simBox.sizeY(), 0),
						Vector3(tiltFactors[1], tiltFactors[2], simBox.sizeZ()),
						simBox.minc - Point3::Origin()));
				break;
			}
			else if(stream.lineStartsWith("ITEM: BOX BOUNDS")) {
				// Parse optional boundary condition flags.
				QStringList tokens = stream.lineString().mid(qstrlen("ITEM: BOX BOUNDS")).split(ws_re, QString::SkipEmptyParts);
				if(tokens.size() >= 3)
					frameData->simulationCell().setPbcFlags(tokens[0] == "pp", tokens[1] == "pp", tokens[2] == "pp");

				// Parse orthogonal simulation box size.
				Box3 simBox;
				for(int k = 0; k < 3; k++) {
					if(sscanf(stream.readLine(), FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING, &simBox.minc[k], &simBox.maxc[k]) != 2)
						throw Exception(tr("Invalid box size in line %1 of LAMMPS dump local file: %2").arg(stream.lineNumber()).arg(stream.lineString()));
				}

				frameData->simulationCell().setMatrix(AffineTransformation(
						Vector3(simBox.sizeX(), 0, 0),
						Vector3(0, simBox.sizeY(), 0),
						Vector3(0, 0, simBox.sizeZ()),
						simBox.minc - Point3::Origin()));
				break;
			}
			else if(stream.lineStartsWith("ITEM: ENTRIES")) {

				// Read the column names list.
				QStringList tokens = stream.lineString().split(ws_re, QString::SkipEmptyParts);
				OVITO_ASSERT(tokens[0] == "ITEM:" && tokens[1] == "ENTRIES");
				QStringList fileColumnNames = tokens.mid(2);

				// Stop here if we are only inspecting the file's header.
				if(_parseFileHeaderOnly) {
					if(fileColumnNames.isEmpty()) {
						// If no file columns names are available, count at least the number
						// of data columns.
						stream.readLine();
						int columnCount = stream.lineString().split(ws_re, QString::SkipEmptyParts).size();
						frameData->detectedColumnMapping().resize(columnCount);
					}
					else {
						frameData->detectedColumnMapping().resize(fileColumnNames.size());
						for(int i = 0; i < fileColumnNames.size(); i++)
							frameData->detectedColumnMapping()[i].columnName = fileColumnNames[i];
					}
					return frameData;
				}

				// Parse data columns.
				InputColumnReader columnParser(_columnMapping, frameData->bonds(), numElements);

				// If possible, use memory-mapped file access for best performance.
				const char* s_start;
				const char* s_end;
				std::tie(s_start, s_end) = stream.mmap();
				auto s = s_start;
				int lineNumber = stream.lineNumber() + 1;
				try {
					for(size_t i = 0; i < numElements; i++, lineNumber++) {
						if(!setProgressValueIntermittent(i)) return {};
						if(!s)
							columnParser.readElement(i, stream.readLine());
						else
							s = columnParser.readElement(i, s, s_end);
					}
				}
				catch(Exception& ex) {
					throw ex.prependGeneralMessage(tr("Parsing error in line %1 of LAMMPS dump local file.").arg(lineNumber));
				}
				if(s) {
					stream.munmap();
					stream.seek(stream.byteOffset() + (s - s_start));
				}

				// Sort the type lists since we created elements on the go and their order depends on the occurrence of types in the file.
				columnParser.sortElementTypes();

				// If the bond "Topology" property was loaded, we need to shift particle indices by 1, because LAMMPS
				// uses 1-based atom IDs and OVITO uses 0-based indices.
				if(PropertyAccess<ParticleIndexPair> topologyProperty = frameData->bonds().findStandardProperty(BondsObject::TopologyProperty)) {
					for(ParticleIndexPair& ab : topologyProperty) {
						ab[0] -= 1;
						ab[1] -= 1;
					}
				}

				// Detect if there are more simulation frames following in the file.
				if(!stream.eof()) {
					stream.readLine();
					if(stream.lineStartsWith("ITEM: TIMESTEP") || stream.lineStartsWith("ITEM: TIME"))
						frameData->signalAdditionalFrames();
				}

				frameData->setStatus(tr("%1 bonds at timestep %2").arg(numElements).arg(timestep));
				return frameData; // Done!
			}
			else if(stream.lineStartsWith("ITEM:")) {
				// For the sake of forward compatibility, we ignore unknown ITEM sections.
				// Skip lines until the next "ITEM:" is reached.
				while(!stream.eof() && !isCanceled()) {
					stream.readLine();
					if(stream.lineStartsWith("ITEM:"))
						break;
				}
			}
			else {
				throw Exception(tr("LAMMPS dump local file parsing error. Line %1 of file %2 is invalid.").arg(stream.lineNumber()).arg(stream.filename()));
			}
		}
		while(!stream.eof());
	}

	throw Exception(tr("LAMMPS dump local file parsing error. Unexpected end of file at line %1 or \"ITEM: ENTRIES\" section is not present in dump file.").arg(stream.lineNumber()));
}

}	// End of namespace
}	// End of namespace
