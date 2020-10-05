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
#include "PDBImporter.h"

#include <3rdparty/gemmi/pdb.hpp>

namespace gemmi { namespace pdb_impl {
template<>
inline size_t copy_line_from_stream<Ovito::CompressedTextReader&>(char* line, int size, Ovito::CompressedTextReader& stream) 
{
	// Return no line if end of file has been reached.
	if(stream.eof())
		return 0;

	// Read a single line form the input stream.
	const char* src_line = stream.readLine();

	// Stop reading the file when ENDMDL marker is reached. We don't want Gemmi to read all frames of a trajectory file.
	if(is_record_type(src_line, "ENDMDL")) {
		return 0;
	}

	// Copy line contents to output buffer.
	size_t len = qstrlen(src_line);
	qstrncpy(line, src_line, size);

	if(is_record_type(src_line, "ATOM") || is_record_type(src_line, "HETATM")) {
		// Some PDB files have an ATOM or HETATM line that is shorter than what Gemmi's parser expects.
		// Pad such lines by appending  additional whitespace.
		if(len < 66 && len >= 55 && size > 66) {
			while(len < 66)
				line[len++] = ' ';
			line[len] = '\0';
		}

		// Gemmi expects atom names to start at column index 12. Some files have one extra space at this positions and the 
		// name actually begins at position 13. Make the parser happy by moving the text by one positon to the left.
		// For example, turn " Au " into "Au  ", but preserve " CA " or " HE ".
		if(len >= 16 && size >= 16 && line[12] == ' ' && line[13] >= 'A' && line[13] <= 'Z' && line[14] >= 'a' && line[14] <= 'z' && line[15] == ' ') {
			line[12] = line[13];
			line[13] = line[14];
			line[14] = ' ';
			line[15] = ' ';
		}
		// Some files have 2 extra spaces at this positions and the name actually begins at position 14. Make the parser happy by moving the text by two characters to the left.
		// For example, turn "  O " into "O   ":
		else if(len >= 16 && size >= 16 && line[12] == ' ' && line[13] == ' ' && line[14] >= 'A' && line[14] <= 'Z') {
			line[12] = line[14];
			line[13] = line[15];
			line[14] = ' ';
			line[15] = ' ';
		}
		// Some files have a digit prepended to the element name. Remove it so that Gemmi can recognize the element correctly.
		// For example, turn "1HH1" into " HH1":
		else if(len >= 16 && size >= 16 && line[12] >= '1' && line[12] <= '9' && line[13] >= 'A' && line[13] <= 'Z') {
			line[12] = ' ';
		}
	}

	// Return line length (up to maximum) to caller.
	return std::min(len, (size_t)(size - 1));
}
}}

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(PDBImporter);

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool PDBImporter::OOMetaClass::checkFileFormat(const FileHandle& file) const
{
	// Open input file.
	CompressedTextReader stream(file);

	// Read up to 40 lines from the beginning of the file.
	for(int i = 0; i < 40 && !stream.eof(); i++) {
		stream.readLine(256);
		if(qstrlen(stream.line()) > 83 && !stream.lineStartsWithToken("TITLE"))
			return false;
		if(qstrlen(stream.line()) >= 7 && stream.line()[6] != ' ' && std::find(stream.line(), stream.line()+6, ' ') != stream.line()+6)
			return false;
		if(stream.lineStartsWithToken("HEADER") || stream.lineStartsWithToken("ATOM") || stream.lineStartsWith("HETATM"))
			return true;
	}

	return false;
}

/******************************************************************************
* Scans the data file and builds a list of source frames.
******************************************************************************/
void PDBImporter::FrameFinder::discoverFramesInFile(QVector<FileSourceImporter::Frame>& frames)
{
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Scanning PDB file %1").arg(stream.filename()));
	setProgressMaximum(stream.underlyingSize());

	Frame frame(fileHandle());
	frame.byteOffset = stream.byteOffset();
	frame.lineNumber = stream.lineNumber();
	while(!stream.eof()) {

		if(isCanceled())
			return;

		stream.readLine();

		if(!setProgressValueIntermittent(stream.underlyingByteOffset()))
			return;

		if(stream.lineStartsWithToken("ENDMDL")) {
			frames.push_back(frame);
			frame.byteOffset = stream.byteOffset();
			frame.lineNumber = stream.lineNumber();
		}
	}

	if(frames.empty()) {
		// It's not a trajectory file. Report just a single frame.
		frames.push_back(Frame(fileHandle()));
	}
}

/******************************************************************************
* Parses the given input file.
******************************************************************************/
FileSourceImporter::FrameDataPtr PDBImporter::FrameLoader::loadFile()
{
	// Open file for reading.
	CompressedTextReader stream(fileHandle());
	setProgressText(tr("Reading PDB file %1").arg(fileHandle().toString()));

	// Jump to byte offset.
	if(frame().byteOffset != 0)
		stream.seek(frame().byteOffset, frame().lineNumber);

	// Create the destination container for loaded data.
	std::shared_ptr<ParticleFrameData> frameData = std::make_shared<ParticleFrameData>();

#if 0
	// Parse metadata records.
	int numAtoms = 0;
	bool hasSimulationCell = false;
	while(!stream.eof()) {

		if(isCanceled())
			return {};

		stream.readLine();
		int lineLength = qstrlen(stream.line());
		if(lineLength < 3 || (lineLength > 83 && !stream.lineStartsWithToken("TITLE")))
			throw Exception(tr("Invalid line length detected in Protein Data Bank (PDB) file at line %1").arg(stream.lineNumber()));

		// Parse simulation cell.
		if(stream.lineStartsWithToken("CRYST1")) {
			FloatType a,b,c,alpha,beta,gamma;
			if(sscanf(stream.line(), "CRYST1 " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " "
					FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING " " FLOATTYPE_SCANF_STRING, &a, &b, &c, &alpha, &beta, &gamma) != 6)
				throw Exception(tr("Invalid simulation cell in Protein Data Bank (PDB) file at line %1").arg(stream.lineNumber()));
			AffineTransformation cell = AffineTransformation::Identity();
			if(alpha == 90 && beta == 90 && gamma == 90) {
				cell(0,0) = a;
				cell(1,1) = b;
				cell(2,2) = c;
			}
			else if(alpha == 90 && beta == 90) {
				gamma *= FLOATTYPE_PI / 180;
				cell(0,0) = a;
				cell(0,1) = b * cos(gamma);
				cell(1,1) = b * sin(gamma);
				cell(2,2) = c;
			}
			else {
				alpha *= FLOATTYPE_PI / 180;
				beta *= FLOATTYPE_PI / 180;
				gamma *= FLOATTYPE_PI / 180;
				FloatType v = a*b*c*sqrt(1.0 - cos(alpha)*cos(alpha) - cos(beta)*cos(beta) - cos(gamma)*cos(gamma) + 2.0 * cos(alpha) * cos(beta) * cos(gamma));
				cell(0,0) = a;
				cell(0,1) = b * cos(gamma);
				cell(1,1) = b * sin(gamma);
				cell(0,2) = c * cos(beta);
				cell(1,2) = c * (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma);
				cell(2,2) = v / (a*b*sin(gamma));
			}
			frameData->simulationCell().setMatrix(cell);
			hasSimulationCell = true;
		}
		else if(stream.lineStartsWithToken("ATOM") || stream.lineStartsWith("HETATM")) {
			// Count atoms.
			numAtoms++;
		}
		else if(stream.lineStartsWithToken("END") || stream.lineStartsWithToken("ENDMDL")) {
			// Stop
			break;
		}
	}

	setProgressMaximum(numAtoms);

	// Jump back to beginning of file.
	stream.seek(frame().byteOffset, frame().lineNumber);

	// Create the particle properties.
	PropertyAccess<Point3> posProperty = frameData->particles().createStandardProperty<ParticlesObject>(numAtoms, ParticlesObject::PositionProperty, false);
	PropertyAccess<int> typeProperty = frameData->particles().createStandardProperty<ParticlesObject>(numAtoms, ParticlesObject::TypeProperty, false);
	PropertyContainerImportData::TypeList* typeList = frameData->particles().createPropertyTypesList(typeProperty, ParticleType::OOClass());

	// Parse atoms.
	size_t atomIndex = 0;
	Point3* p = posProperty.begin();
	int* a = typeProperty.begin();
	PropertyAccess<qlonglong> particleIdentifierProperty;
	PropertyAccess<qlonglong> moleculeIdentifierProperty;
	PropertyAccess<int> moleculeTypeProperty;
	PropertyContainerImportData::TypeList* moleculeTypeList = nullptr;
	while(!stream.eof() && atomIndex < numAtoms) {
		if(!setProgressValueIntermittent(atomIndex))
			return {};

		stream.readLine();
		int lineLength = qstrlen(stream.line());
		if(lineLength < 3 || (lineLength > 83 && !stream.lineStartsWithToken("TITLE")))
			throw Exception(tr("Invalid line length detected in Protein Data Bank (PDB) file at line %1").arg(stream.lineNumber()));

		// Parse atom definition.
		if(stream.lineStartsWithToken("ATOM") || stream.lineStartsWith("HETATM")) {
			char atomType[4];
			int atomTypeLength = 0;
			for(const char* c = stream.line() + 76; c <= stream.line() + std::min(77, lineLength); ++c)
				if(*c > ' ') atomType[atomTypeLength++] = *c;
			if(atomTypeLength == 0) {
				for(const char* c = stream.line() + 12; c <= stream.line() + std::min(15, lineLength); ++c)
					if(*c > ' ') atomType[atomTypeLength++] = *c;
			}
			*a = typeList->addTypeName(atomType, atomType + atomTypeLength);
#ifdef FLOATTYPE_FLOAT
			if(lineLength <= 30 || sscanf(stream.line() + 30, "%8g%8g%8g", &p->x(), &p->y(), &p->z()) != 3)
#else
			if(lineLength <= 30 || sscanf(stream.line() + 30, "%8lg%8lg%8lg", &p->x(), &p->y(), &p->z()) != 3)
#endif
				throw Exception(tr("Invalid atom coordinates (line %1): %2").arg(stream.lineNumber()).arg(stream.lineString()));

			// Parse atom ID (serial number).
			qlonglong atomSerialNumber;
			if(sscanf(stream.line() + 6, "%5llu", &atomSerialNumber) == 1) {
				if(!particleIdentifierProperty) {
					particleIdentifierProperty = frameData->particles().createStandardProperty<ParticlesObject>(numAtoms, ParticlesObject::IdentifierProperty, true);
				}
				particleIdentifierProperty[atomIndex] = atomSerialNumber;
			}
			else if(particleIdentifierProperty && qstrncmp(stream.line() + 6, "*****", 5) == 0) {
				// This is special handling for large PDB files with more than 99,999 atoms.
				// Some codes replace the 5 digits in the 'atom serial number' column with the string '*****' in this case.
				// We we encounter this case, we simply assign consecutive IDs to the atoms.
				particleIdentifierProperty[atomIndex] = atomIndex + 1;
			}

			// Parse molecule ID (residue sequence number).
			qlonglong residueSequenceNumber;
			if(sscanf(stream.line() + 22, "%4llu", &residueSequenceNumber) == 1) {
				if(!moleculeIdentifierProperty) {
					moleculeIdentifierProperty = frameData->particles().createStandardProperty<ParticlesObject>(numAtoms, ParticlesObject::MoleculeProperty, true);
				}
				moleculeIdentifierProperty[atomIndex] = residueSequenceNumber;
			}

			// Parse molecule type.
			char moleculeType[3];
			int moleculeTypeLength = 0;
			for(const char* c = stream.line() + 17; c <= stream.line() + std::min(19, lineLength); ++c)
				if(*c != ' ') moleculeType[moleculeTypeLength++] = *c;
			if(moleculeTypeLength != 0) {
				if(!moleculeTypeProperty) {
					moleculeTypeProperty = frameData->particles().createStandardProperty<ParticlesObject>(numAtoms, ParticlesObject::MoleculeTypeProperty, true);
					moleculeTypeList = frameData->particles().createPropertyTypesList(moleculeTypeProperty, ElementType::OOClass());
				}
				moleculeTypeProperty[atomIndex] = moleculeTypeList->addTypeName(moleculeType, moleculeType + moleculeTypeLength);
			}

			atomIndex++;
			++p;
			++a;
		}
	}

	// Parse bonds.
	PropertyAccess<ParticleIndexPair> bondTopologyProperty;
	while(!stream.eof()) {

		if(isCanceled())
			return {};

		stream.readLine();
		int lineLength = qstrlen(stream.line());
		if(lineLength < 3 || (lineLength > 83 && !stream.lineStartsWithToken("TITLE")))
			throw Exception(tr("Invalid line length detected in Protein Data Bank (PDB) file at line %1").arg(stream.lineNumber()));

		// Parse bonds.
		if(stream.lineStartsWith("CONECT")) {

			qlonglong atomSerialNumber;
			if(lineLength <= 8 || !particleIdentifierProperty) 
				throw Exception(tr("Invalid CONECT record (line %1): %2").arg(stream.lineNumber()).arg(stream.lineString()));
			const char* atomidstring = stream.line() + 6; // Skip the CONECT keyword.
			bool startfound = false;
			bool firstnumber = true;
			int startnumber;
			int ndigits;
			int start;
			qlonglong atomIndex1;
			qlonglong atomIndex2;
			int i = 0;
			while(atomidstring[i] >= ' ' && i < 67) { // Loop over all characters in the string.
				if(atomidstring[i] == ' ' || atomidstring[i+1] < ' ') {
					if(atomidstring[i] != ' ') i++;
					if(startfound) { // We have found the end of a number
						start = startnumber;
						ndigits = i - startnumber;
						if(ndigits <= 5 || numAtoms > 99999) { // Normal case: space-separated numbers
							startfound = false;
						}
						else { // Most likely we have a classic pdb file where the atomID fields are stuck together 
							ndigits -= 5 * ((ndigits-1)/5);
							startnumber += ndigits;
							// Repeat cycle so the next number of the concatenated numbers is read.
							if(atomidstring[i] == ' ') i -= 1;
							else i -= 2;
						}
						if(!parseInt64(atomidstring + start, atomidstring + start + ndigits, atomSerialNumber))
							throw Exception(tr("Invalid CONECT record (line %1): %2").arg(ndigits).arg(stream.lineNumber()).arg(stream.lineString()));
						if(firstnumber) {
							firstnumber = false;
							atomIndex1 = std::distance(particleIdentifierProperty.cbegin(), std::find(particleIdentifierProperty.cbegin(), particleIdentifierProperty.cend(), atomSerialNumber));
							if((size_t)atomIndex1 >= particleIdentifierProperty.size())
								throw Exception(tr("Nonexistent atom ID %1 encountered in line %2 of PDB file.").arg(atomSerialNumber).arg(stream.lineNumber()));
						} 
						else {
					        atomIndex2 = std::distance(particleIdentifierProperty.cbegin(), std::find(particleIdentifierProperty.cbegin(), particleIdentifierProperty.cend(), atomSerialNumber));
							if((size_t)atomIndex2 >= particleIdentifierProperty.size())
								throw Exception(tr("Nonexistent atom ID %1 encountered in line %2 of PDB file.").arg(atomSerialNumber).arg(stream.lineNumber()));
							if(!bondTopologyProperty)
								bondTopologyProperty = frameData->bonds().createStandardProperty<BondsObject>(1, BondsObject::TopologyProperty, false);
							else
					        	bondTopologyProperty.storage()->resize(bondTopologyProperty.size() + 1, true);
					        bondTopologyProperty[bondTopologyProperty.size() - 1] = ParticleIndexPair{{atomIndex1, atomIndex2}};
						}
					}
				}
				else if(!startfound) { // We have found the start of a number
					startnumber = i;
					startfound = true;
				}
				else if(atomidstring[i] < ' ') {
					break;
				}
				i++;
			}
		}
		else if(stream.lineStartsWithToken("END") || stream.lineStartsWithToken("ENDMDL")) {
			break;
		}
	}
	// Detect if there are more simulation frames following in the file.
	for(int i = 0; i < 18; i++) {
		if(stream.eof()) break;
		stream.readLine();
		if(stream.lineStartsWithToken("MODEL") || stream.lineStartsWithToken("REMARK") || stream.lineStartsWithToken("TITLE")) {
			frameData->signalAdditionalFrames();
			break;
		}
	}

	// If file does not contain any simulation cell info,
	// compute bounding box of atoms and use it as an adhoc simulation cell.
	if(!hasSimulationCell && numAtoms > 0) {
		Box3 boundingBox;
		boundingBox.addPoints(posProperty);
		frameData->simulationCell().setPbcFlags(false, false, false);
		frameData->simulationCell().setMatrix(AffineTransformation(
				Vector3(boundingBox.sizeX(), 0, 0),
				Vector3(0, boundingBox.sizeY(), 0),
				Vector3(0, 0, boundingBox.sizeZ()),
				boundingBox.minc - Point3::Origin()));
	}


	if(bondTopologyProperty)
		frameData->generateBondPeriodicImageProperty();

#else

	try {
		// Parse the PDB file's contents.
		gemmi::Structure structure = gemmi::pdb_impl::read_pdb_from_line_input(stream, qPrintable(frame().sourceFile.path()), gemmi::PdbReadOptions());
		if(isCanceled()) return {};

		structure.merge_chain_parts();
		if(isCanceled()) return {};

		if(structure.models.empty())
			throw Exception(tr("PDB parsing error: No structural models."));
		const gemmi::Model& model = structure.models.back();

		// Count total number of atoms.
		size_t natoms = 0;
		for(const gemmi::Chain& chain : model.chains) {
			for(const gemmi::Residue& residue : chain.residues) {
				natoms += residue.atoms.size();
			}
		}

		// Allocate property arrays for atoms.
		PropertyAccess<Point3> posProperty = frameData->particles().createStandardProperty<ParticlesObject>(natoms, ParticlesObject::PositionProperty, false);
		PropertyAccess<int> typeProperty = frameData->particles().createStandardProperty<ParticlesObject>(natoms, ParticlesObject::TypeProperty, false);
		PropertyContainerImportData::TypeList* typeList = frameData->particles().createPropertyTypesList(typeProperty, ParticleType::OOClass());
		Point3* posIter = posProperty.begin();
		int* typeIter = typeProperty.begin();

		// Transfer atomic data from Gemmi to OVITO data structures.
		bool hasOccupancy = false;
		for(const gemmi::Chain& chain : model.chains) {
			for(const gemmi::Residue& residue : chain.residues) {
				if(isCanceled()) return {};
				for(const gemmi::Atom& atom : residue.atoms) {
					// Atomic position.
					*posIter++ = Point3(atom.pos.x, atom.pos.y, atom.pos.z);

					// Atomic type.
					*typeIter++ = atom.element.ordinal();
					if(!typeList->hasTypeId(atom.element.ordinal()))
						typeList->addNamedTypeId(atom.element.ordinal(), QString::fromStdString(atom.element.name()), false);

					// Check for presence of occupancy values.
					if(atom.occ != 0 && atom.occ != 1) hasOccupancy = true;
				}
			}
		}
		if(isCanceled()) return {};

		// Parse the optional site occupancy information.
		if(hasOccupancy) {
			PropertyAccess<FloatType> occupancyProperty = frameData->particles().addProperty(std::make_shared<PropertyStorage>(natoms, PropertyStorage::Float, 1, 0, QStringLiteral("Occupancy"), false));
			FloatType* occupancyIter = occupancyProperty.begin();
			for(const gemmi::Chain& chain : model.chains) {
				for(const gemmi::Residue& residue : chain.residues) {
					for(const gemmi::Atom& atom : residue.atoms) {
						*occupancyIter++ = atom.occ;
					}
				}
			}
			OVITO_ASSERT(occupancyIter == occupancyProperty.end());
		}

		// Since we created particle types on the go while reading the particles, the assigned particle type IDs
		// depend on the storage order of particles in the file We rather want a well-defined particle type ordering, that's
		// why we sort them now.
		typeList->sortTypesById();

		// Parse unit cell.
		if(structure.cell.is_crystal()) {
			// Process periodic unit cell definition.
			AffineTransformation cell = AffineTransformation::Identity();
			if(structure.cell.alpha == 90 && structure.cell.beta == 90 && structure.cell.gamma == 90) {
				cell(0,0) = structure.cell.a;
				cell(1,1) = structure.cell.b;
				cell(2,2) = structure.cell.c;
			}
			else if(structure.cell.alpha == 90 && structure.cell.beta == 90) {
				FloatType gamma = qDegreesToRadians(structure.cell.gamma);
				cell(0,0) = structure.cell.a;
				cell(0,1) = structure.cell.b * std::cos(gamma);
				cell(1,1) = structure.cell.b * std::sin(gamma);
				cell(2,2) = structure.cell.c;
			}
			else {
				FloatType alpha = qDegreesToRadians(structure.cell.alpha);
				FloatType beta = qDegreesToRadians(structure.cell.beta);
				FloatType gamma = qDegreesToRadians(structure.cell.gamma);
				FloatType v = structure.cell.a * structure.cell.b * structure.cell.c * sqrt(1.0 - std::cos(alpha)*std::cos(alpha) - std::cos(beta)*std::cos(beta) - std::cos(gamma)*std::cos(gamma) + 2.0 * std::cos(alpha) * std::cos(beta) * std::cos(gamma));
				cell(0,0) = structure.cell.a;
				cell(0,1) = structure.cell.b * std::cos(gamma);
				cell(1,1) = structure.cell.b * std::sin(gamma);
				cell(0,2) = structure.cell.c * std::cos(beta);
				cell(1,2) = structure.cell.c * (std::cos(alpha) - std::cos(beta)*std::cos(gamma)) / std::sin(gamma);
				cell(2,2) = v / (structure.cell.a * structure.cell.b * std::sin(gamma));
			}
			frameData->simulationCell().setMatrix(cell);
		}
		else if(posProperty.size() != 0) {
			// Use bounding box of atomic coordinates as non-periodic simulation cell.
			Box3 boundingBox;
			boundingBox.addPoints(posProperty);
			frameData->simulationCell().setPbcFlags(false, false, false);
			frameData->simulationCell().setMatrix(AffineTransformation(
					Vector3(boundingBox.sizeX(), 0, 0),
					Vector3(0, boundingBox.sizeY(), 0),
					Vector3(0, 0, boundingBox.sizeZ()),
					boundingBox.minc - Point3::Origin()));
		}
		frameData->setStatus(tr("Number of atoms: %1").arg(natoms));
	}
	catch(const std::exception& e) {
		throw Exception(tr("PDB file error: %1").arg(e.what()));
	}

	// Check if more frames are following in the trajectory file.
	if(!stream.eof()) {
		stream.readLine();
		if(!stream.eof()) {
			frameData->signalAdditionalFrames();
		}
	}

#endif

	return frameData;
}

}	// End of namespace
}	// End of namespace
