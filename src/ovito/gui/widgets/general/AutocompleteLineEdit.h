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

/**
 * \file AutocompleteLineEdit.h
 * \brief Contains the definition of the Ovito::AutocompleteLineEdit class.
 */

#pragma once


#include <ovito/gui/GUI.h>

namespace Ovito { OVITO_BEGIN_INLINE_NAMESPACE(Gui) OVITO_BEGIN_INLINE_NAMESPACE(Widgets)

/**
 * \brief A text editor widget that provides auto-completion of words.
 */
class OVITO_GUI_EXPORT AutocompleteLineEdit : public QLineEdit
{
	Q_OBJECT

public:

	/// \brief Constructs the widget.
	AutocompleteLineEdit(QWidget* parent = nullptr);

	/// Sets the list of words that can be completed.
	void setWordList(const QStringList& words) { _wordListModel->setStringList(words); }

protected Q_SLOTS:

	/// Inserts a complete word into the text field.
	void onComplete(const QString& completion);

protected:

	/// Handles key-press events.
	virtual void keyPressEvent(QKeyEvent* event) override;

	/// Creates a list of tokens from the current text string.
	QStringList getTokenList() const;

protected:

	/// The completer object used by the widget.
	QCompleter* _completer;

	/// The list model storing the words that can be completed.
	QStringListModel* _wordListModel;

	/// Regular expression used to split a text into words.
	QRegularExpression _wordSplitter;
};

OVITO_END_INLINE_NAMESPACE
OVITO_END_INLINE_NAMESPACE
}	// End of namespace


