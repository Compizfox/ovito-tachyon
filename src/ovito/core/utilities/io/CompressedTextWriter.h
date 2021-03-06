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

#pragma once


#include <ovito/core/Core.h>
#ifdef OVITO_ZLIB_SUPPORT
	#include <ovito/core/utilities/io/gzdevice/GzipIODevice.h>
#endif

namespace Ovito {

/**
 * \brief A helper class for writing text-based files that are compressed (gzip format).
 *
 * If the destination filename has a .gz suffix, this output stream class compresses the
 * text data on the fly if.
 *
 * \sa CompressedTextReader
 */
class OVITO_CORE_EXPORT CompressedTextWriter : public QObject
{
public:

	/// Opens the given output file device for writing. The output file's name determines if the data gets compressed.
	/// \param output The underlying Qt output device from to data should be written.
	/// \param context The DataSet which provides a context for error messages generated by this writer.
	/// \throw Exception if an I/O error has occurred.
	CompressedTextWriter(QFileDevice& output, DataSet* context = nullptr);

	/// Returns the name of the output file.
	const QString& filename() const { return _filename; }

	/// Returns the underlying I/O device.
	QFileDevice& device() { return _device; }

	/// Returns whether data written to this stream is being compressed.
	bool isCompressed() const { return _stream != &_device; }

	/// Writes an integer number to the text-based output file.
	CompressedTextWriter& operator<<(qint32 i);

	/// Writes an unsigned integer number to the text-based output file.
	CompressedTextWriter& operator<<(quint32 i);

	/// Writes a 64-bit integer number to the text-based output file.
	CompressedTextWriter& operator<<(qint64 i);

	/// Writes a 64-bit unsigned integer number to the text-based output file.
	CompressedTextWriter& operator<<(quint64 i);

#if (!defined(Q_OS_WIN) && (QT_POINTER_SIZE != 4)) || defined(Q_OS_WASM)
	/// Writes an unsigned integer number to the text-based output file.
	CompressedTextWriter& operator<<(size_t i);
#endif

	/// Writes a floating-point number to the text-based output file.
	CompressedTextWriter& operator<<(FloatType f);

	/// Writes a text string to the text-based output file.
	CompressedTextWriter& operator<<(const char* s) {
		if(_stream->write(s) == -1)
			reportWriteError();
		return *this;
	}

	/// Writes a single character to the text-based output file.
	CompressedTextWriter& operator<<(char c) {
		if(!_stream->putChar(c))
			reportWriteError();
		return *this;
	}

	/// Writes a Qt string string to the text-based output file.
	CompressedTextWriter& operator<<(const QString& s) { return *this << s.toLocal8Bit().constData(); }

	/// Returns the current output precision for floating-point numbers.
	unsigned int floatPrecision() const { return _floatPrecision; }

	/// Changes the output precision for floating-point numbers.
	void setFloatPrecision(unsigned int precision) {
		_floatPrecision = std::min(precision, (unsigned int)std::numeric_limits<FloatType>::max_digits10);
	}

private:

	/// Throws an exception to report an I/O error.
	void reportWriteError();

	/// The name of the output file (if known).
	QString _filename;

	/// The underlying output device.
	QFileDevice& _device;

#ifdef OVITO_ZLIB_SUPPORT
	/// The compression filter stream.
	GzipIODevice _compressor;
#endif

	/// The output stream.
	QIODevice* _stream;

	/// The dataset to which this writer belongs. Is used when throwing exceptions.
	DataSet* _context;

	/// The output precision for floating-point numbers.
	unsigned int _floatPrecision = 10;

	Q_OBJECT
};

}	// End of namespace
