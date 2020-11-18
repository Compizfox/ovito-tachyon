////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2013 Alexander Stukowski
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

#include <ovito/core/Core.h>
#include <ovito/core/rendering/FrameBuffer.h>
#ifdef OVITO_VIDEO_OUTPUT_SUPPORT
	#include <ovito/core/utilities/io/video/VideoEncoder.h>
#endif

namespace Ovito {

#define IMAGE_FORMAT_FILE_FORMAT_VERSION		1

/******************************************************************************
* Constructor.
******************************************************************************/
FrameBuffer::FrameBuffer(int width, int height, QObject* parent) : QObject(parent), _image(width, height, QImage::Format_ARGB32)
{
	_info.setImageWidth(width);
	_info.setImageHeight(height);
	clear();
}

/******************************************************************************
* Detects the file format based on the filename suffix.
******************************************************************************/
bool ImageInfo::guessFormatFromFilename()
{
	if(filename().endsWith(QStringLiteral(".png"), Qt::CaseInsensitive)) {
		setFormat("png");
		return true;
	}
	else if(filename().endsWith(QStringLiteral(".jpg"), Qt::CaseInsensitive) || filename().endsWith(QStringLiteral(".jpeg"), Qt::CaseInsensitive)) {
		setFormat("jpg");
		return true;
	}
#ifdef OVITO_VIDEO_OUTPUT_SUPPORT
	for(const auto& videoFormat : VideoEncoder::supportedFormats()) {
		for(const QString& extension : videoFormat.extensions) {
			if(filename().endsWith(QStringLiteral(".") + extension, Qt::CaseInsensitive)) {
				setFormat(videoFormat.name);
				return true;
			}
		}
	}
#endif

	return false;
}

/******************************************************************************
* Returns whether the selected file format is a video format.
******************************************************************************/
bool ImageInfo::isMovie() const
{
#ifdef OVITO_VIDEO_OUTPUT_SUPPORT
	for(const auto& videoFormat : VideoEncoder::supportedFormats()) {
		if(format() == videoFormat.name)
			return true;
	}
#endif

	return false;
}

/******************************************************************************
* Writes an ImageInfo to an output stream.
******************************************************************************/
SaveStream& operator<<(SaveStream& stream, const ImageInfo& i)
{
	stream.beginChunk(IMAGE_FORMAT_FILE_FORMAT_VERSION);
	stream << i._imageWidth;
	stream << i._imageHeight;
	stream << i._filename;
	stream << i._format;
	stream.endChunk();
	return stream;
}

/******************************************************************************
* Reads an ImageInfo from an input stream.
******************************************************************************/
LoadStream& operator>>(LoadStream& stream, ImageInfo& i)
{
	stream.expectChunk(IMAGE_FORMAT_FILE_FORMAT_VERSION);
	stream >> i._imageWidth;
	stream >> i._imageHeight;
	stream >> i._filename;
	stream >> i._format;
	stream.closeChunk();
	return stream;
}

/******************************************************************************
* Clears the framebuffer with a uniform color.
******************************************************************************/
void FrameBuffer::clear(const ColorA& color)
{
	_image.fill(color);
}

/******************************************************************************
* Removes unnecessary pixels along the outer edges of the image.
******************************************************************************/
void FrameBuffer::autoCrop()
{
	QImage image = this->image().convertToFormat(QImage::Format_ARGB32);
	if(image.width() <= 0 || image.height() <= 0) return;

	auto determineCropRect = [&image](QRgb backgroundColor) -> QRect {
		int x1 = 0, y1 = 0;
		int x2 = image.width() - 1, y2 = image.height() - 1;
		bool significant;
		for(;; x1++) {
			significant = false;
			for(int y = y1; y <= y2; y++) {
				if(reinterpret_cast<const QRgb*>(image.constScanLine(y))[x1] != backgroundColor) {
					significant = true;
					break;
				}
			}
			if(significant || x1 > x2)
				break;
		}
		for(; x2 >= x1; x2--) {
			significant = false;
			for(int y = y1; y <= y2; y++) {
				if(reinterpret_cast<const QRgb*>(image.constScanLine(y))[x2] != backgroundColor) {
					significant = true;
					break;
				}
			}
			if(significant || x1 > x2)
				break;
		}
		for(;; y1++) {
			significant = false;
			const QRgb* s = reinterpret_cast<const QRgb*>(image.constScanLine(y1));
			for(int x = x1; x <= x2; x++) {
				if(s[x] != backgroundColor) {
					significant = true;
					break;
				}
			}
			if(significant || y1 >= y2)
				break;
		}
		for(; y2 >= y1; y2--) {
			significant = false;
			const QRgb* s = reinterpret_cast<const QRgb*>(image.constScanLine(y2));
			for(int x = x1; x <= x2; x++) {
				if(s[x] != backgroundColor) {
					significant = true;
					break;
				}
			}
			if(significant || y1 > y2)
				break;
		}
		return QRect(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
	};

	// Use the pixel colors in the four corners of the images as candidate background colors.
	// Compute the crop rect for each candidate color and use the one that leads
	// to the smallest crop rect.
	QRect cropRect = determineCropRect(image.pixel(0,0));
	QRect r;
	r = determineCropRect(image.pixel(image.width()-1,0));
	if(r.width()*r.height() < cropRect.width()*cropRect.height()) cropRect = r;
	r = determineCropRect(image.pixel(image.width()-1,image.height()-1));
	if(r.width()*r.height() < cropRect.width()*cropRect.height()) cropRect = r;
	r = determineCropRect(image.pixel(0,image.height()-1));
	if(r.width()*r.height() < cropRect.width()*cropRect.height()) cropRect = r;

	if(cropRect != image.rect() && cropRect.width() > 0 && cropRect.height() > 0) {
		_image = _image.copy(cropRect);
		update();
	}
}

}	// End of namespace
