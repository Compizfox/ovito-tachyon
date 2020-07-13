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
#include "LAMMPSTextDumpImporter.h"

#include <QRegularExpression>

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(LAMMPSTextDumpImporter);
DEFINE_PROPERTY_FIELD(LAMMPSTextDumpImporter, useCustomColumnMapping);
SET_PROPERTY_FIELD_LABEL(LAMMPSTextDumpImporter, useCustomColumnMapping, "Custom file column mapping");

/******************************************************************************
 * Sets the user-defined mapping between data columns in the input file and
 * the internal particle properties.
 *****************************************************************************/
void LAMMPSTextDumpImporter::setCustomColumnMapping(const InputColumnMapping& mapping)
{
	_customColumnMapping = mapping;
	notifyTargetChanged();
}

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool LAMMPSTextDumpImporter::OOMetaClass::checkFileFormat(const FileHandle& file) const
{
	// Open input file.
	CompressedTextReader stream(file);

	// Read first line.
	stream.readLine(15);
	if(stream.lineStartsWith("ITEM: TIMESTEP"))
		return true;

	return false;
}

/******************************************************************************
* Inspects the header of the given file and returns the number of file columns.
******************************************************************************/
Future<InputColumnMapping> LAMMPSTextDumpImporter::inspectFileHeader(const Frame& frame)
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
void LAMMPSTextDumpImporter::FrameFinder::discoverFramesInFile(QVector<FileSourceImporter::Frame>& frames)
{
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Scanning LAMMPS dump file %1").arg(fileHandle().toString()));
	setProgressMaximum(stream.underlyingSize());

	// Regular expression for whitespace characters.
	QRegularExpression ws_re(QStringLiteral("\\s+"));

	int timestep = 0;
	size_t numParticles = 0;
	Frame frame(fileHandle());

	while(!stream.eof() && !isCanceled()) {
		qint64 byteOffset = stream.byteOffset();
		int lineNumber = stream.lineNumber();

		// Parse next line.
		stream.readLine();

		do {
			if(stream.lineStartsWith("ITEM: TIMESTEP")) {
				if(sscanf(stream.readLine(), "%i", &timestep) != 1)
					throw Exception(tr("LAMMPS dump file parsing error. Invalid timestep number (line %1):\n%2").arg(stream.lineNumber()).arg(QString::fromLocal8Bit(stream.line())));
				frame.byteOffset = byteOffset;
				frame.lineNumber = lineNumber;
				frame.label = QString("Timestep %1").arg(timestep);
				frames.push_back(frame);
				break;
			}
			else if(stream.lineStartsWith("ITEM: NUMBER OF ATOMS")) {
				// Parse number of atoms.
				unsigned long long u;
				if(sscanf(stream.readLine(), "%llu", &u) != 1)
					throw Exception(tr("LAMMPS dump file parsing error. Invalid number of atoms in line %1:\n%2").arg(stream.lineNumber()).arg(stream.lineString()));
				if(u > 100'000'000'000ll)
					throw Exception(tr("LAMMPS dump file parsing error. Number of atoms in line %1 is too large. The LAMMPS dump file reader doesn't accept files with more than 100 billion atoms.").arg(stream.lineNumber()));
				numParticles = (size_t)u;
				break;
			}
			else if(stream.lineStartsWith("ITEM: ATOMS")) {
				for(size_t i = 0; i < numParticles; i++) {
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
				throw Exception(tr("LAMMPS dump file parsing error. Line %1 of file %2 is invalid.").arg(stream.lineNumber()).arg(stream.filename()));
			}
		}
		while(!stream.eof());
	}
}

/******************************************************************************
* Parses the given input file.
******************************************************************************/
FileSourceImporter::FrameDataPtr LAMMPSTextDumpImporter::FrameLoader::loadFile()
{
	// Open file for reading.
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Reading LAMMPS dump file %1").arg(fileHandle().toString()));

	// Jump to byte offset.
	if(frame().byteOffset != 0)
		stream.seek(frame().byteOffset, frame().lineNumber);

	// Create the destination container for loaded data.
	auto frameData = std::make_shared<LAMMPSFrameData>();

	// Regular expression for whitespace characters.
	QRegularExpression ws_re(QStringLiteral("\\s+"));

	int timestep;
	size_t numParticles = 0;

	while(!stream.eof()) {

		// Parse next line.
		stream.readLine();

		do {
			if(stream.lineStartsWith("ITEM: TIMESTEP")) {
				if(sscanf(stream.readLine(), "%i", &timestep) != 1)
					throw Exception(tr("LAMMPS dump file parsing error. Invalid timestep number (line %1):\n%2").arg(stream.lineNumber()).arg(stream.lineString()));
				frameData->attributes().insert(QStringLiteral("Timestep"), QVariant::fromValue(timestep));
				break;
			}
			else if(stream.lineStartsWith("ITEM: NUMBER OF ATOMS")) {
				// Parse number of atoms.
				unsigned long long u;
				if(sscanf(stream.readLine(), "%llu", &u) != 1)
					throw Exception(tr("LAMMPS dump file parsing error. Invalid number of atoms in line %1:\n%2").arg(stream.lineNumber()).arg(stream.lineString()));
				if(u >= 2147483648ull)
					throw Exception(tr("LAMMPS dump file parsing error. Number of atoms in line %1 exceeds internal limit of 2^31 atoms:\n%2").arg(stream.lineNumber()).arg(stream.lineString()));

				numParticles = (size_t)u;
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
						throw Exception(tr("Invalid box size in line %1 of LAMMPS dump file: %2").arg(stream.lineNumber()).arg(stream.lineString()));
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
						throw Exception(tr("Invalid box size in line %1 of dump file: %2").arg(stream.lineNumber()).arg(stream.lineString()));
				}

				frameData->simulationCell().setMatrix(AffineTransformation(
						Vector3(simBox.sizeX(), 0, 0),
						Vector3(0, simBox.sizeY(), 0),
						Vector3(0, 0, simBox.sizeZ()),
						simBox.minc - Point3::Origin()));
				break;
			}
			else if(stream.lineStartsWith("ITEM: ATOMS")) {

				// Read the column names list.
				QStringList tokens = stream.lineString().split(ws_re, QString::SkipEmptyParts);
				OVITO_ASSERT(tokens[0] == "ITEM:" && tokens[1] == "ATOMS");
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
						frameData->detectedColumnMapping() = generateAutomaticColumnMapping(fileColumnNames);
					}
					return frameData;
				}

				// Set up column-to-property mapping.
				InputColumnMapping columnMapping;
				if(_useCustomColumnMapping)
					columnMapping = _customColumnMapping;
				else
					columnMapping = generateAutomaticColumnMapping(fileColumnNames);

				// Parse data columns.
				InputColumnReader columnParser(columnMapping, *frameData, numParticles);

				// If possible, use memory-mapped file access for best performance.
				const char* s_start;
				const char* s_end;
				std::tie(s_start, s_end) = stream.mmap();
				auto s = s_start;
				int lineNumber = stream.lineNumber() + 1;
				try {
					for(size_t i = 0; i < numParticles; i++, lineNumber++) {
						if(!setProgressValueIntermittent(i)) return {};
						if(!s)
							columnParser.readParticle(i, stream.readLine());
						else
							s = columnParser.readParticle(i, s, s_end);
					}
				}
				catch(Exception& ex) {
					throw ex.prependGeneralMessage(tr("Parsing error in line %1 of LAMMPS dump file.").arg(lineNumber));
				}
				if(s) {
					stream.munmap();
					stream.seek(stream.byteOffset() + (s - s_start));
				}

				// Sort the particle type list since we created particles on the go and their order depends on the occurrence of types in the file.
				columnParser.sortParticleTypes();

				// Determine if particle coordinates are given in reduced form and need to be rescaled to absolute form.
				bool reducedCoordinates = false;
				if(!fileColumnNames.empty()) {
					// If the dump file contains column names, then we can use them to detect 
					// the type of particle coordinates. Reduced coordinates are found in columns
					// "xs, ys, zs" or "xsu, ysu, zsu".
					for(int i = 0; i < (int)columnMapping.size() && i < fileColumnNames.size(); i++) {
						if(columnMapping[i].property.type() == ParticlesObject::PositionProperty) {
							reducedCoordinates = (
									fileColumnNames[i] == "xs" || fileColumnNames[i] == "xsu" ||
									fileColumnNames[i] == "ys" || fileColumnNames[i] == "ysu" ||
									fileColumnNames[i] == "zs" || fileColumnNames[i] == "zsu");
							// break; Note: Do not stop the loop here, because the 'Position' particle 
							// property may be associated with several file columns, and it's the last column that 
							// ends up getting imported into OVITO. 
						}
					}
				}
				else {
					// If no column names are available, use the following heuristic:
					// Assume reduced coordinates if all particle coordinates are within the [-0.02,1.02] interval.
					// We allow coordinates to be slightly outside the [0,1] interval, because LAMMPS
					// wraps around particles at the periodic boundaries only occasionally.
					if(ConstPropertyAccess<Point3> posProperty = frameData->findStandardParticleProperty(ParticlesObject::PositionProperty)) {
						// Compute bound box of particle positions.
						Box3 boundingBox;
						boundingBox.addPoints(posProperty);
						// Check if bounding box is inside the (slightly extended) unit cube.
						if(Box3(Point3(FloatType(-0.02)), Point3(FloatType(1.02))).containsBox(boundingBox))
							reducedCoordinates = true;
					}
				}

				if(reducedCoordinates) {
					// Convert all atom coordinates from reduced to absolute (Cartesian) format.
					if(PropertyAccess<Point3> posProperty = frameData->findStandardParticleProperty(ParticlesObject::PositionProperty)) {
						const AffineTransformation simCell = frameData->simulationCell().matrix();
						for(Point3& p : posProperty)
							p = simCell * p;
					}
				}

				// If a "diameter" column was loaded and stored in the "Radius" particle property,
				// we need to divide values by two.
				if(!fileColumnNames.empty()) {
					for(int i = 0; i < (int)columnMapping.size() && i < fileColumnNames.size(); i++) {
						if(columnMapping[i].property.type() == ParticlesObject::RadiusProperty && fileColumnNames[i] == "diameter") {
							if(PropertyAccess<FloatType> radiusProperty = frameData->findStandardParticleProperty(ParticlesObject::RadiusProperty)) {
								for(FloatType& r : radiusProperty)
									r /= 2;
							}
							break;
						}
					}
				}

				// Detect dimensionality of system.
				frameData->simulationCell().set2D(!columnMapping.hasZCoordinates());

				// Detect if there are more simulation frames following in the file.
				if(!stream.eof()) {
					stream.readLine();
					if(stream.lineStartsWith("ITEM: TIMESTEP"))
						frameData->signalAdditionalFrames();
				}

				// Sort particles by ID.
				if(_sortParticles)
					frameData->sortParticlesById();

				frameData->setStatus(tr("%1 particles at timestep %2").arg(numParticles).arg(timestep));
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
				throw Exception(tr("LAMMPS dump file parsing error. Line %1 of file %2 is invalid.").arg(stream.lineNumber()).arg(stream.filename()));
			}
		}
		while(!stream.eof());
	}

	throw Exception(tr("LAMMPS dump file parsing error. Unexpected end of file at line %1 or \"ITEM: ATOMS\" section is not present in dump file.").arg(stream.lineNumber()));
}

/******************************************************************************
 * Guesses the mapping of input file columns to internal particle properties.
 *****************************************************************************/
InputColumnMapping LAMMPSTextDumpImporter::generateAutomaticColumnMapping(const QStringList& columnNames)
{
	InputColumnMapping columnMapping;
	columnMapping.resize(columnNames.size());
	for(int i = 0; i < columnNames.size(); i++) {
		QString name = columnNames[i].toLower();
		columnMapping[i].columnName = columnNames[i];
		if(name == "x" || name == "xu" || name == "coordinates") columnMapping.mapStandardColumn(i, ParticlesObject::PositionProperty, 0);
		else if(name == "y" || name == "yu") columnMapping.mapStandardColumn(i, ParticlesObject::PositionProperty, 1);
		else if(name == "z" || name == "zu") columnMapping.mapStandardColumn(i, ParticlesObject::PositionProperty, 2);
		else if(name == "xs" || name == "xsu") { columnMapping.mapStandardColumn(i, ParticlesObject::PositionProperty, 0); }
		else if(name == "ys" || name == "ysu") { columnMapping.mapStandardColumn(i, ParticlesObject::PositionProperty, 1); }
		else if(name == "zs" || name == "zsu") { columnMapping.mapStandardColumn(i, ParticlesObject::PositionProperty, 2); }
		else if(name == "vx" || name == "velocities") columnMapping.mapStandardColumn(i, ParticlesObject::VelocityProperty, 0);
		else if(name == "vy") columnMapping.mapStandardColumn(i, ParticlesObject::VelocityProperty, 1);
		else if(name == "vz") columnMapping.mapStandardColumn(i, ParticlesObject::VelocityProperty, 2);
		else if(name == "id") columnMapping.mapStandardColumn(i, ParticlesObject::IdentifierProperty);
		else if(name == "type" || name == "element" || name == "atom_types") columnMapping.mapStandardColumn(i, ParticlesObject::TypeProperty);
		else if(name == "mass") columnMapping.mapStandardColumn(i, ParticlesObject::MassProperty);
		else if(name == "radius" || name == "diameter") columnMapping.mapStandardColumn(i, ParticlesObject::RadiusProperty);
		else if(name == "mol") columnMapping.mapStandardColumn(i, ParticlesObject::MoleculeProperty);
		else if(name == "q") columnMapping.mapStandardColumn(i, ParticlesObject::ChargeProperty);
		else if(name == "ix") columnMapping.mapStandardColumn(i, ParticlesObject::PeriodicImageProperty, 0);
		else if(name == "iy") columnMapping.mapStandardColumn(i, ParticlesObject::PeriodicImageProperty, 1);
		else if(name == "iz") columnMapping.mapStandardColumn(i, ParticlesObject::PeriodicImageProperty, 2);
		else if(name == "fx" || name == "forces") columnMapping.mapStandardColumn(i, ParticlesObject::ForceProperty, 0);
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
		else if(name == "c_cna" || name == "pattern") columnMapping.mapStandardColumn(i, ParticlesObject::StructureTypeProperty);
		else if(name == "c_epot") columnMapping.mapStandardColumn(i, ParticlesObject::PotentialEnergyProperty);
		else if(name == "c_kpot") columnMapping.mapStandardColumn(i, ParticlesObject::KineticEnergyProperty);
		else if(name == "c_stress[1]") columnMapping.mapStandardColumn(i, ParticlesObject::StressTensorProperty, 0);
		else if(name == "c_stress[2]") columnMapping.mapStandardColumn(i, ParticlesObject::StressTensorProperty, 1);
		else if(name == "c_stress[3]") columnMapping.mapStandardColumn(i, ParticlesObject::StressTensorProperty, 2);
		else if(name == "c_stress[4]") columnMapping.mapStandardColumn(i, ParticlesObject::StressTensorProperty, 3);
		else if(name == "c_stress[5]") columnMapping.mapStandardColumn(i, ParticlesObject::StressTensorProperty, 4);
		else if(name == "c_stress[6]") columnMapping.mapStandardColumn(i, ParticlesObject::StressTensorProperty, 5);
		else if(name == "c_orient[1]") columnMapping.mapStandardColumn(i, ParticlesObject::OrientationProperty, 0);
		else if(name == "c_orient[2]") columnMapping.mapStandardColumn(i, ParticlesObject::OrientationProperty, 1);
		else if(name == "c_orient[3]") columnMapping.mapStandardColumn(i, ParticlesObject::OrientationProperty, 2);
		else if(name == "c_orient[4]") columnMapping.mapStandardColumn(i, ParticlesObject::OrientationProperty, 3);
		else if(name == "c_shape[1]") columnMapping.mapStandardColumn(i, ParticlesObject::AsphericalShapeProperty, 0);
		else if(name == "c_shape[2]") columnMapping.mapStandardColumn(i, ParticlesObject::AsphericalShapeProperty, 1);
		else if(name == "c_shape[3]") columnMapping.mapStandardColumn(i, ParticlesObject::AsphericalShapeProperty, 2);
		else if(name == "selection") columnMapping.mapStandardColumn(i, ParticlesObject::SelectionProperty, 0);
		else {
			columnMapping.mapCustomColumn(i, name, PropertyStorage::Float);
		}
	}
	return columnMapping;
}

/******************************************************************************
 * Saves the class' contents to the given stream.
 *****************************************************************************/
void LAMMPSTextDumpImporter::saveToStream(ObjectSaveStream& stream, bool excludeRecomputableData)
{
	ParticleImporter::saveToStream(stream, excludeRecomputableData);

	stream.beginChunk(0x01);
	_customColumnMapping.saveToStream(stream);
	stream.endChunk();
}

/******************************************************************************
 * Loads the class' contents from the given stream.
 *****************************************************************************/
void LAMMPSTextDumpImporter::loadFromStream(ObjectLoadStream& stream)
{
	ParticleImporter::loadFromStream(stream);

	stream.expectChunk(0x01);
	_customColumnMapping.loadFromStream(stream);
	stream.closeChunk();
}

/******************************************************************************
 * Creates a copy of this object.
 *****************************************************************************/
OORef<RefTarget> LAMMPSTextDumpImporter::clone(bool deepCopy, CloneHelper& cloneHelper) const
{
	// Let the base class create an instance of this class.
	OORef<LAMMPSTextDumpImporter> clone = static_object_cast<LAMMPSTextDumpImporter>(ParticleImporter::clone(deepCopy, cloneHelper));
	clone->_customColumnMapping = this->_customColumnMapping;
	return clone;
}


}	// End of namespace
}	// End of namespace
