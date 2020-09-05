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


#include <ovito/gui/desktop/GUI.h>

namespace Ovito {

class OverlayListModel;	// Defined in OverlayListModel.h

class OverlayAction : public QAction
{
	Q_OBJECT

public:

	/// Constructs an action for a built-in layer class.
	static OverlayAction* createForClass(OvitoClassPtr clazz);

	/// Constructs an action for a Python overlay script.
	static OverlayAction* createForScript(const QString& fileName, const QDir& directory);

	/// Returns the overlay class descriptor if this action represents a built-in overlay type. 
	OvitoClassPtr layerClass() const { return _layerClass; }

	/// The absolute path of the modifier script if this action represents a Python-based modifier function.
	const QString& scriptPath() const { return _scriptPath; }

private:

	/// The Ovito class descriptor of the viewport layer subclass. 
	OvitoClassPtr _layerClass = nullptr;

	/// The path to the overlay script on disk.
	QString _scriptPath;
};

/**
 * A Qt list model that list all available viewport layer types.
 */
class OverlayTypesModel : public QAbstractListModel
{
	Q_OBJECT

public:

	/// Constructor.
	OverlayTypesModel(QObject* parent, MainWindow* mainWindow, OverlayListModel* overlayListModel);

	/// Returns the number of rows in the model.
	virtual int rowCount(const QModelIndex& parent) const override;

	/// Returns the data associated with a list item.
	virtual QVariant data(const QModelIndex& index, int role) const override;

	/// Returns the flags for an item.
	virtual Qt::ItemFlags flags(const QModelIndex& index) const override;

	/// Returns the action that belongs to the given model index.
	OverlayAction* actionFromIndex(int index) const;

	/// Returns the action that belongs to the given model index.
	OverlayAction* actionFromIndex(const QModelIndex& index) const { return actionFromIndex(index.row()); }

private Q_SLOTS:

	/// Signal handler that inserts the selected viewport layer into the active viewport.
	void insertViewportLayer();

private:

	/// The list of viewport layer actions.
	std::vector<OverlayAction*> _actions;

	/// The context main window.
	MainWindow* _mainWindow;

	/// The model representing the viewport layers of the active viewport.
	OverlayListModel* _overlayListModel;

	/// The list of directories searched for user-defined viewport layer scripts.
	QVector<QDir> _layerScriptDirectories;
};

}	// End of namespace
