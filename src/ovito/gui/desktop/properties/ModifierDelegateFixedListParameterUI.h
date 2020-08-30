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
#include "RefTargetListParameterUI.h"

namespace Ovito {

/******************************************************************************
* A list view that shows the fixed set of delegates of a MultiDelegatingModifier,
* which can each be enabled or disabled by the user.
******************************************************************************/
class OVITO_GUI_EXPORT ModifierDelegateFixedListParameterUI : public RefTargetListParameterUI
{
	Q_OBJECT
	OVITO_CLASS(ModifierDelegateFixedListParameterUI)

public:

	/// Constructor.
	ModifierDelegateFixedListParameterUI(QObject* parentEditor,
			const RolloutInsertionParameters& rolloutParams = RolloutInsertionParameters(), OvitoClassPtr defaultEditorClass = nullptr);

	/// This method is called when a new editable object has been activated.
	virtual void resetUI() override {
		RefTargetListParameterUI::resetUI();
		// Clear initial selection by default.
		listWidget()->selectionModel()->clear();
	}

protected:

	/// Returns a data item from the list data model.
	virtual QVariant getItemData(RefTarget* target, const QModelIndex& index, int role) override;

	/// Returns the model/view item flags for the given entry.
	virtual Qt::ItemFlags getItemFlags(RefTarget* target, const QModelIndex& index) override;

	/// Sets the role data for the item at index to value.
	virtual bool setItemData(RefTarget* target, const QModelIndex& index, const QVariant& value, int role) override;

	/// Returns the number of columns for the table view.
	virtual int tableColumnCount() override { return 1; }

	/// Returns the header data under the given role for the given RefTarget.
	virtual QVariant getHorizontalHeaderData(int index, int role) override {
		if(role == Qt::DisplayRole) {
			return QVariant::fromValue(tr("Data type"));
		}
		return RefTargetListParameterUI::getHorizontalHeaderData(index, role);
	}
};

}	// End of namespace
