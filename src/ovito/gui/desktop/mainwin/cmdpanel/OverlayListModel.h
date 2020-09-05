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
#include <ovito/core/oo/RefTarget.h>
#include <ovito/core/oo/RefTargetListener.h>
#include "OverlayListItem.h"

namespace Ovito {

/**
 * A Qt model class used to populate the QListView widget on the viewports overlay page of the command panel.
 */
class OverlayListModel : public QAbstractListModel
{
	Q_OBJECT

public:

	/// Constructor.
	OverlayListModel(QObject* parent);

	/// Returns the number of list items.
	virtual int rowCount(const QModelIndex& parent = QModelIndex()) const override { return _items.size(); }

	/// Returns the data associated with a list entry.
	virtual QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;

	/// Changes the data associated with a list entry.
	virtual bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole) override;

	/// Returns the flags for an item.
	virtual Qt::ItemFlags flags(const QModelIndex& index) const override;

	/// Returns the associated selection model.
	QItemSelectionModel* selectionModel() const { return _selectionModel; }

	/// Returns the currently selected model item in list.
	OverlayListItem* selectedItem() const;

	/// Returns the currently selected index in the overlay list.
	int selectedIndex() const;

	/// Returns an item from the list model.
	OverlayListItem* item(int index) const {
		OVITO_ASSERT(index >= 0 && index < _items.size());
		return _items[index];
	}

	/// Populates the model with the given list items.
	void setItems(const QList<OORef<OverlayListItem>>& newItems);

	/// Returns the list of items.
	const QList<OORef<OverlayListItem>>& items() const { return _items; }

	/// The currently selected viewport whose overlays are mirrored by this list model.
	Viewport* selectedViewport() const { return _selectedViewport.target(); }

	/// Sets the currently selected viewport whose overlays should be mirrored by this list model.
	void setSelectedViewport(Viewport* viewport) {
		_selectedViewport.setTarget(viewport);
		refreshList();
	}

	/// Sets the overlay that should be selected on the next list update.
	void setNextToSelectObject(ViewportOverlay* obj) { _nextObjectToSelect = obj; }

Q_SIGNALS:

	/// This signal is emitted if a new list item has been selected, or if the currently
	/// selected item has changed.
	void selectedItemChanged();

public Q_SLOTS:

	/// Rebuilds the viewport overlay list.
	void refreshList();

	/// Updates the appearance of a single list item.
	void refreshItem(OverlayListItem* item);

	/// Rebuilds the list of items as soon as possible.
	void requestUpdate() {
		if(_needListUpdate) return;	// Update is already pending.
		_needListUpdate = true;
		// Invoke actual refresh function at some later time.
		QMetaObject::invokeMethod(this, "refreshList", Qt::QueuedConnection);
	}

private Q_SLOTS:

	/// Handles notification events generated by the selected viewport.
	void onViewportEvent(const ReferenceEvent& event);

private:

	/// List of visible items in the model.
	QList<OORef<OverlayListItem>> _items;

	/// Holds reference to the currently selected Viewport.
	RefTargetListener<Viewport> _selectedViewport;

	/// The item in the list that should be selected on the next list update.
	ViewportOverlay* _nextObjectToSelect = nullptr;

	/// The selection model of the list view widget.
	QItemSelectionModel* _selectionModel;

	/// Indicates that the list of items needs to be updated.
	bool _needListUpdate = false;

	// Status icons:
	QPixmap _statusInfoIcon;
	QPixmap _statusWarningIcon;
	QPixmap _statusErrorIcon;
	QPixmap _statusNoneIcon;

	/// Font used for section headers.
	QFont _sectionHeaderFont;

	/// The background brush used for list section headers.
	QBrush _sectionHeaderBackgroundBrush;

	/// The foreground brush used for list section headers.
	QBrush _sectionHeaderForegroundBrush;
};

}	// End of namespace
