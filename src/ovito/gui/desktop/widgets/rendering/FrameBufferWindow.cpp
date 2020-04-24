////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2018 Alexander Stukowski
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

#include <ovito/gui/desktop/GUI.h>
#include <ovito/gui/desktop/dialogs/SaveImageFileDialog.h>
#include "FrameBufferWindow.h"

namespace Ovito {

/******************************************************************************
* Constructor.
******************************************************************************/
FrameBufferWindow::FrameBufferWindow(QWidget* parent) :
	QMainWindow(parent, (Qt::WindowFlags)(Qt::Tool | Qt::CustomizeWindowHint | Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint))
{
	// Note: This line has been disabled, because it leads to sporadic program crashes (Qt 5.12.5).
//	setAttribute(Qt::WA_MacAlwaysShowToolWindow);

	class MyScrollArea : public QScrollArea {
	public:
		MyScrollArea(QWidget* parent) : QScrollArea(parent) {}
		virtual QSize sizeHint() const override {
			int f = 2 * frameWidth();
			QSize sz(f, f);
			if(widget())
				sz += widget()->sizeHint();
			return sz;
		}
	};

	QScrollArea* scrollArea = new MyScrollArea(this);
	_frameBufferWidget = new FrameBufferWidget();
	scrollArea->setWidget(_frameBufferWidget);
	setCentralWidget(scrollArea);

	QToolBar* toolBar = addToolBar(tr("Frame Buffer"));
	toolBar->addAction(QIcon(":/gui/framebuffer/save_picture.bw.svg"), tr("Save to file"), this, SLOT(saveImage()));
	toolBar->addAction(QIcon(":/gui/framebuffer/copy_picture_to_clipboard.bw.svg"), tr("Copy to clipboard"), this, SLOT(copyImageToClipboard()));
	toolBar->addSeparator();
	toolBar->addAction(QIcon(":/gui/framebuffer/auto_crop.bw.svg"), tr("Auto-crop image"), this, SLOT(autoCrop()));

	// Disable context menu in toolbar.
	setContextMenuPolicy(Qt::NoContextMenu);
}

/******************************************************************************
* Creates a frame buffer of the requested size and adjusts the size of the window.
******************************************************************************/
const std::shared_ptr<FrameBuffer>& FrameBufferWindow::createFrameBuffer(int w, int h)
{
	// Allocate and resize frame buffer and frame buffer window if necessary.
	if(!frameBuffer()) {
		setFrameBuffer(std::make_shared<FrameBuffer>(w, h));
	}
	if(frameBuffer()->size() != QSize(w, h)) {
		frameBuffer()->setSize(QSize(w, h));
		frameBuffer()->clear();
		resize(sizeHint());
	}
	return frameBuffer();
}

/******************************************************************************
* Shows and activates the frame buffer window.
******************************************************************************/
void FrameBufferWindow::showAndActivateWindow()
{
	if(isHidden()) {
		// Center frame buffer window in main window.
		if(parentWidget()) {
			QSize s = frameGeometry().size();
			move(parentWidget()->geometry().center() - QPoint(s.width() / 2, s.height() / 2));
		}
		show();
		updateGeometry();
		update();
	}
	activateWindow();
}

/******************************************************************************
* This opens the file dialog and lets the suer save the current contents of the frame buffer
* to an image file.
******************************************************************************/
void FrameBufferWindow::saveImage()
{
	if(!frameBuffer())
		return;

	SaveImageFileDialog fileDialog(this, tr("Save image"));
	if(fileDialog.exec()) {
		QString imageFilename = fileDialog.imageInfo().filename();
		if(!frameBuffer()->image().save(imageFilename, fileDialog.imageInfo().format())) {
			Exception ex(tr("Failed to save image to file '%1'.").arg(imageFilename));
			ex.reportError();
		}
	}
}

/******************************************************************************
* This copies the current image to the clipboard.
******************************************************************************/
void FrameBufferWindow::copyImageToClipboard()
{
	if(!frameBuffer())
		return;

	QClipboard *clipboard = QApplication::clipboard();
	clipboard->setImage(frameBuffer()->image());
}

/******************************************************************************
* Removes unnecessary pixels along the outer edges of the rendered image.
******************************************************************************/
void FrameBufferWindow::autoCrop()
{
	if(frameBuffer())
		frameBuffer()->autoCrop();
}

}	// End of namespace
