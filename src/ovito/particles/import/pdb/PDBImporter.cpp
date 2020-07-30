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
#include <ovito/particles/import/ParticleFrameData.h>
#include <ovito/core/utilities/io/CompressedTextReader.h>
#include "PDBImporter.h"

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(PDBImporter);

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool PDBImporter::OOMetaClass::checkFileFormat(const FileHandle& file) const
{
	// Open input file.
	CompressedTextReader stream(file);

	// Read the first N lines.
	for(int i = 0; i < 20 && !stream.eof(); i++) {
		stream.readLine(86);
		if(qstrlen(stream.line()) > 83 && !stream.lineStartsWithToken("TITLE"))
			return false;
		if(qstrlen(stream.line()) >= 7 && stream.line()[6] != ' ')
			return false;
		if(stream.lineStartsWithToken("HEADER") || stream.lineStartsWithToken("ATOM") || QString::compare(stream.lineString().left(6),"HETATM",Qt::CaseInsensitive) ==0)
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
		int lineLength = qstrlen(stream.line());
		if(lineLength < 3 || (lineLength > 83 && !stream.lineStartsWithToken("TITLE")))
			throw Exception(tr("Invalid line length detected in Protein Data Bank (PDB) file at line %1").arg(stream.lineNumber()));

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
		else if(stream.lineStartsWithToken("ATOM") ||  QString::compare(stream.lineString().left(6),"HETATM",Qt::CaseInsensitive) ==0) {
			// Count atoms.
			numAtoms++;
		}
		else if(stream.lineStartsWithToken("TER") || stream.lineStartsWithToken("END") || stream.lineStartsWithToken("ENDMDL")) {
			// Stop
			break;
		}
	}

	setProgressMaximum(numAtoms);

	// Jump back to beginning of file.
	stream.seek(frame().byteOffset, frame().lineNumber);

	// Create the particle properties.
	PropertyAccess<Point3> posProperty = frameData->addParticleProperty(ParticlesObject::OOClass().createStandardStorage(numAtoms, ParticlesObject::PositionProperty, false));
	PropertyAccess<int> typeProperty = frameData->addParticleProperty(ParticlesObject::OOClass().createStandardStorage(numAtoms, ParticlesObject::TypeProperty, false));
	ParticleFrameData::TypeList* typeList = frameData->createPropertyTypesList(typeProperty);

	// Parse atoms.
	size_t atomIndex = 0;
	Point3* p = posProperty.begin();
	int* a = typeProperty.begin();
	PropertyAccess<qlonglong> particleIdentifierProperty;
	PropertyAccess<qlonglong> moleculeIdentifierProperty;
	PropertyAccess<int> moleculeTypeProperty;
	ParticleFrameData::TypeList* moleculeTypeList = nullptr;
	while(!stream.eof() && atomIndex < numAtoms) {
		if(!setProgressValueIntermittent(atomIndex))
			return {};

		stream.readLine();
		int lineLength = qstrlen(stream.line());
		if(lineLength < 3 || (lineLength > 83 && !stream.lineStartsWithToken("TITLE")))
			throw Exception(tr("Invalid line length detected in Protein Data Bank (PDB) file at line %1").arg(stream.lineNumber()));

		// Parse atom definition.
		if(stream.lineStartsWithToken("ATOM") ||  QString::compare(stream.lineString().left(6),"HETATM",Qt::CaseInsensitive) ==0) {
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
					particleIdentifierProperty = frameData->addParticleProperty(ParticlesObject::OOClass().createStandardStorage(numAtoms, ParticlesObject::IdentifierProperty, true));
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
					moleculeIdentifierProperty = frameData->addParticleProperty(ParticlesObject::OOClass().createStandardStorage(numAtoms, ParticlesObject::MoleculeProperty, true));
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
					moleculeTypeProperty = frameData->addParticleProperty(ParticlesObject::OOClass().createStandardStorage(numAtoms, ParticlesObject::MoleculeTypeProperty, true));
					moleculeTypeList = frameData->createPropertyTypesList(moleculeTypeProperty);
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
		// linestartsWithToken assumes a space follows after the to be compared string. In the PDB file format other character are allowed to occur directly after the keyword and this function therefore does not work. 
		if(QString::compare(stream.lineString().left(6),"CONECT",Qt::CaseInsensitive) ==0) {

			qlonglong atomSerialNumber;
			if(lineLength <= 8 || !particleIdentifierProperty) throw Exception(tr("Invalid CONECT record (line %1): %2").arg(stream.lineNumber()).arg(stream.lineString()));
			QString atomidstring = stream.lineString().mid(6); // skip the first 6 char as they hold CONECT
			bool startfound = false;
			bool firstnumber = true;
			int startnumber;
			int ndigits;
			int start;
			size_t atomIndex1;
			size_t atomIndex2;
			bool strtointOK;
			int i = 0;
			while (i < atomidstring.size()){ // loop over all characters in the string
				if (atomidstring[i] == " " || i+1 == atomidstring.size() ){  
					if (startfound){ // we have found the end of a number
						start = startnumber;
						ndigits = i-startnumber;
						if (ndigits<=5 || numAtoms > 99999){ // case normal space separated number
							startfound = false;
						} else { // Most likely we have a classic pdb file where the atomID fields are stuck together 
							ndigits = ndigits-(5*floor((ndigits-1.0)/5.0));
							startnumber = startnumber+ndigits;
							i=i-1; // repeat cycle so the next number of the concatenated numbers is read.
						}
						atomSerialNumber = atomidstring.mid(start,ndigits).toLongLong(&strtointOK,10);
						if (!strtointOK){
							throw Exception(tr("Invalid CONECT record (line %1): %2").arg(ndigits).arg(stream.lineNumber()).arg(stream.lineString()));
						}
						if (firstnumber){
							firstnumber = false;
							atomIndex1 = boost::find(particleIdentifierProperty, atomSerialNumber) - particleIdentifierProperty.cbegin();
							if(atomIndex1 >= particleIdentifierProperty.size())
								throw Exception(tr("Nonexistent atom ID:%1 encountered in line %2 of PDB file.").arg(atomIndex1).arg(stream.lineNumber()));
						} else {
					        atomIndex2 = boost::find(particleIdentifierProperty, atomSerialNumber) - particleIdentifierProperty.cbegin();
							if(atomIndex1 >= particleIdentifierProperty.size() || atomIndex2 >= particleIdentifierProperty.size())
								throw Exception(tr("Nonexistent atom ID:%1 encountered in line %2 of PDB file.").arg(atomIndex2).arg(stream.lineNumber()));
							if(!bondTopologyProperty) {
								bondTopologyProperty = frameData->addBondProperty(BondsObject::OOClass().createStandardStorage(1, BondsObject::TopologyProperty, false));
					        }else {
					        	bondTopologyProperty.storage()->resize(bondTopologyProperty.size() + 1, true);
					        }
					        bondTopologyProperty[bondTopologyProperty.size() - 1][0] = atomIndex1;
				            bondTopologyProperty[bondTopologyProperty.size() - 1][1] = atomIndex2;
						}
					}
				}else { 
					if (startfound ==false) { // We have found the start of a number
						startnumber = i;
						startfound = true;
					}
				}
				i++;
			}	
		}
		else if(stream.lineStartsWithToken("END") || stream.lineStartsWithToken("TER") || stream.lineStartsWithToken("ENDMDL")) {
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

	// If file does not contains simulation cell info,
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

	frameData->setStatus(tr("Number of atoms: %1").arg(numAtoms));
	return frameData;
}


}	// End of namespace
}	// End of namespace
