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
#include "XTCImporter.h"

#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(XTCImporter);

class XTCFile
{
public:

	struct Frame {
		int step;
		float time;
		float prec;
		Matrix_3<float> cell;
		std::vector<Point_3<float>> xyz;
	};

	~XTCFile() {
		close();
	}

	void open(const char* filename) {
		close();
		int returnCode = read_xtc_natoms(filename, &_numAtoms);
		if(returnCode != exdrOK || _numAtoms <= 0)
			throw Exception(XTCImporter::tr("Error opening XTC file (error code %1).").arg(returnCode));
		_file = xdrfile_open(filename, "r");
		if(!_file)
			throw Exception(XTCImporter::tr("Error opening XTC file."));
		_eof = false;
	}

	void close() {
		if(_file) {
			if(xdrfile_close(_file) != exdrOK)
				qWarning() << "XTCImporter: Failure reported by xdrfile_close()";
			_file = nullptr;
		}
	}

	Frame read() {
		OVITO_ASSERT(_file);
		OVITO_ASSERT(!_eof);
		Frame frame;
		frame.xyz.resize(_numAtoms);
		int returnCode = read_xtc(_file, _numAtoms, &frame.step, &frame.time, 
			reinterpret_cast<matrix&>(frame.cell), 
			reinterpret_cast<rvec*>(frame.xyz.data()), &frame.prec);
		if(returnCode != exdrOK && returnCode != exdrENDOFFILE)
			throw Exception(XTCImporter::tr("Error reading XTC file (code %1).").arg(returnCode));
		if(returnCode == exdrENDOFFILE)
			_eof = true;
		return frame;
	}

	void seek(qint64 offset) {
		int returnCode = xdr_seek(_file, offset, SEEK_SET);
        if(returnCode != exdrOK)
			throw Exception(XTCImporter::tr("Error seeking in XTC file (code %1).").arg(returnCode));
	}

	qint64 byteOffset() const {
		OVITO_ASSERT(_file);
		return xdr_tell(_file);
	}

	bool eof() const { return _eof; }

private:

	XDRFILE* _file = nullptr;
	int _numAtoms = 0;
	bool _eof = false;
};

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool XTCImporter::OOMetaClass::checkFileFormat(const FileHandle& file) const
{
	return file.sourceUrl().fileName().endsWith(QStringLiteral(".xtc"), Qt::CaseInsensitive);
}

/******************************************************************************
* Scans the data file and builds a list of source frames.
******************************************************************************/
void XTCImporter::FrameFinder::discoverFramesInFile(QVector<FileSourceImporter::Frame>& frames)
{
	setProgressText(tr("Scanning file %1").arg(fileHandle().toString()));
	setProgressMaximum(QFileInfo(fileHandle().localFilePath()).size());

	// Open XTC file for reading.
	XTCFile file;
	file.open(fileHandle().localFilePath().toUtf8().constData());

	Frame frame(fileHandle());
	while(!file.eof() && !isCanceled()) {
		frame.byteOffset = file.byteOffset();
		if(!setProgressValue(frame.byteOffset))
			return;

		// Parse trajectory frame.
		XTCFile::Frame xtcFrame = file.read();
		if(file.eof()) break;

		// Create a new record for the timestep.
		frame.label = tr("Timestep %1").arg(xtcFrame.step);
		frames.push_back(frame);
	}
}

/******************************************************************************
* Parses the given input file.
******************************************************************************/
FileSourceImporter::FrameDataPtr XTCImporter::FrameLoader::loadFile()
{
	// Open file for reading.
	setProgressText(tr("Reading XTC file %1").arg(fileHandle().toString()));

	// Create the destination container for loaded data.
	auto frameData = std::make_shared<ParticleFrameData>();

	// Open XTC file for reading.
	XTCFile file;
	file.open(fileHandle().localFilePath().toUtf8().constData());

	// Seek to byte offset of requested trajectory frame.
	if(frame().byteOffset != 0)
		file.seek(frame().byteOffset);

	// Read trajectory frame data.
	XTCFile::Frame xtcFrame = file.read();

	// Transfer atomic coordinates to property storage. Also convert from nanometer units to angstroms.
	size_t numParticles = xtcFrame.xyz.size();
	PropertyAccess<Point3> posProperty = frameData->particles().createStandardProperty<ParticlesObject>(numParticles, ParticlesObject::PositionProperty, false);
	std::transform(xtcFrame.xyz.cbegin(), xtcFrame.xyz.cend(), posProperty.begin(), [](const Point_3<float>& p) {
		return static_cast<Point3>(p * 10.0f);
	});

	// Convert cell vectors from nanometers to angstroms.
	frameData->simulationCell().setMatrix(AffineTransformation(static_cast<Matrix3>(xtcFrame.cell * 10.0f)));

	frameData->attributes().insert(QStringLiteral("Timestep"), QVariant::fromValue(xtcFrame.step));
	frameData->attributes().insert(QStringLiteral("Time"), QVariant::fromValue((FloatType)xtcFrame.time));

	return frameData;
}

}	// End of namespace
}	// End of namespace
