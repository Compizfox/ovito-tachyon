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

#pragma once


#include <ovito/particles/Particles.h>
#include <ovito/particles/import/ParticleImporter.h>
#include <ovito/particles/import/ParticleFrameData.h>

#include <QXmlDefaultHandler>

namespace Ovito { namespace Particles { OVITO_BEGIN_INLINE_NAMESPACE(Import) OVITO_BEGIN_INLINE_NAMESPACE(Formats)

/**
 * \brief File parser for data files of the GALAMOST MD code.
 */
class OVITO_GALAMOST_EXPORT GALAMOSTImporter : public ParticleImporter
{
	/// Defines a metaclass specialization for this importer type.
	class OOMetaClass : public ParticleImporter::OOMetaClass
	{
	public:
		/// Inherit standard constructor from base meta class.
		using ParticleImporter::OOMetaClass ::OOMetaClass;

		/// Returns the file filter that specifies the files that can be imported by this service.
		virtual QString fileFilter() const override { return QStringLiteral("*.xml"); }

		/// Returns the filter description that is displayed in the drop-down box of the file dialog.
		virtual QString fileFilterDescription() const override { return tr("GALAMOST Files"); }

		/// Checks if the given file has format that can be read by this importer.
		virtual bool checkFileFormat(QFileDevice& input, const QUrl& sourceLocation) const override;
	};

	OVITO_CLASS_META(GALAMOSTImporter, OOMetaClass)
	Q_OBJECT

public:

	/// \brief Constructor.
	Q_INVOKABLE GALAMOSTImporter(DataSet *dataset) : ParticleImporter(dataset) {}

	/// Returns the title of this object.
	virtual QString objectTitle() const override { return tr("GALAMOST"); }

	/// Creates an asynchronous loader object that loads the data for the given frame from the external file.
	virtual std::shared_ptr<FileSourceImporter::FrameLoader> createFrameLoader(const Frame& frame, const QString& localFilename) override {
		return std::make_shared<FrameLoader>(frame, localFilename);
	}

private:

	/// The format-specific task object that is responsible for reading an input file in a separate thread.
	class FrameLoader : public FileSourceImporter::FrameLoader, protected QXmlDefaultHandler
	{
	public:

		/// Constructor.
		FrameLoader(const FileSourceImporter::Frame& frame, const QString& filename)
			: FileSourceImporter::FrameLoader(frame, filename) {}

	protected:

		/// Loads the frame data from the given file.
		virtual FrameDataPtr loadFile(QFile& file) override;

		/// Is called by the XML parser whenever a new XML element is read.
		virtual bool startElement(const QString& namespaceURI, const QString& localName, const QString& qName, const QXmlAttributes& atts) override;

		/// Is called by the XML parser whenever it has parsed an end element tag.
		virtual bool endElement(const QString& namespaceURI, const QString& localName, const QString& qName) override;

		/// Is called by the XML parser whenever it has parsed a chunk of character data.
		virtual bool characters(const QString& ch) override;

		/// Is called when a non-recoverable error is encountered during parsing of the XML file.
		virtual bool fatalError(const QXmlParseException& exception) override;

	private:

		/// Container for the parsed particle data.
		std::shared_ptr<ParticleFrameData> _frameData;

		/// The dimensionality of the dataset.
		int _dimensions = 3;

		/// The number of atoms.
		size_t _natoms = 0;

		/// The particle/bond property whose currently being parsed.
		PropertyPtr _currentProperty;

		/// Buffer for text data read from XML file.
		QString _characterData;

		/// The number of <configuration> elements that have been parsed so far.
		size_t _numConfigurationsRead = 0;
	};
};

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace
}	// End of namespace