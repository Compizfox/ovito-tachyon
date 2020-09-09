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

#include <ovito/gui/desktop/GUI.h>
#include <ovito/gui/desktop/actions/ActionManager.h>

namespace Ovito {

void ActionManager::setupCommandSearch()
{
	// Set up QAction that activates quick search.
	QWidgetAction* commandQuickSearchAction = new QWidgetAction(this);
	commandQuickSearchAction->setText(tr("Quick Command Search"));
	commandQuickSearchAction->setObjectName(ACTION_COMMAND_QUICKSEARCH);
#ifndef Q_OS_MAC
	commandQuickSearchAction->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_Q));
#else
	commandQuickSearchAction->setShortcut(QKeySequence(Qt::CTRL | Qt::Key_P));
#endif
	commandQuickSearchAction->setStatusTip(tr("Quickly access program commands."));

	// Sort and filter input list of actions.
	class ActionListModel : public QSortFilterProxyModel {
	public:
		ActionListModel(QObject* parent, QAbstractItemModel* sourceModel) : QSortFilterProxyModel(parent) {
			setSourceModel(sourceModel);
			sort(0, Qt::AscendingOrder);
		}
		bool filterAcceptsRow(int sourceRow, const QModelIndex& sourceParent) const override {
			QAction* action = sourceModel()->index(sourceRow, 0, sourceParent).data(ActionRole).value<QAction*>();
			return action->isVisible();
		}
	};

	// Subclass QLineEdit.
	class SearchField : public QLineEdit {
	public:
		SearchField(ActionManager* actionManager) : _actionManager(actionManager) {
			_completer = new QCompleter(this);
			_completer->setCompletionMode(QCompleter::PopupCompletion);
			_completer->setCaseSensitivity(Qt::CaseInsensitive);
			_completer->setFilterMode(Qt::MatchContains);
			_completer->setModel(new ActionListModel(_completer, actionManager));
			_completer->setCompletionRole(SearchTextRole);
			_completer->setWidget(this);

			class ItemDelegate : public QStyledItemDelegate 
			{
			public:
				ItemDelegate() {
					_tooltipFont = QGuiApplication::font();
				}
			protected:
				QFont _tooltipFont;
				virtual void paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index) const override {
					QStyleOptionViewItem options = option;
					initStyleOption(&options, index);
					options.features |= QStyleOptionViewItem::HasDecoration;
					options.decorationSize = static_cast<const QAbstractItemView*>(option.widget)->iconSize();

					// Draw list item without text content.
					QString text = std::move(options.text);
					options.text.clear();
					options.widget->style()->drawControl(QStyle::CE_ItemViewItem, &options, painter, options.widget);

					// Draw shortcut text.
					options.rect.adjust(0, 4, 0, -4);
					options.backgroundBrush = {};
#ifndef Q_OS_WIN
					// Override text color for highlighted items.
					if(options.state & QStyle::State_Selected)
						options.palette.setColor(QPalette::Text, options.palette.color(QPalette::Active, QPalette::HighlightedText));
#endif
					options.state.setFlag(QStyle::State_Selected, false);
					options.state.setFlag(QStyle::State_MouseOver, false);
					options.icon = {};
					options.displayAlignment = Qt::AlignRight | Qt::AlignVCenter;
					QKeySequence keySequence = index.data(ShortcutRole).value<QKeySequence>();
					if(!keySequence.isEmpty()) {
						options.text = keySequence.toString(QKeySequence::NativeText) + QStringLiteral(" ");
						QFont oldFont = std::move(options.font);
						options.font = _tooltipFont;
						options.widget->style()->drawControl(QStyle::CE_ItemViewItem, &options, painter, options.widget);
						options.font = std::move(oldFont);
						QFontMetrics fm(_tooltipFont);
						options.rect.setWidth(options.rect.width() - fm.boundingRect(options.text).width());
					}

					// Draw first line of text.
					options.displayAlignment = Qt::AlignLeft | Qt::AlignTop;
					options.text = std::move(text);
					options.widget->style()->drawControl(QStyle::CE_ItemViewItem, &options, painter, options.widget);

					// Draw second line of text.
					options.text = index.data(Qt::StatusTipRole).toString();
					options.displayAlignment = Qt::AlignLeft | Qt::AlignBottom;
					options.font = _tooltipFont;
					options.widget->style()->drawControl(QStyle::CE_ItemViewItem, &options, painter, options.widget);
				}

				virtual QSize sizeHint(const QStyleOptionViewItem& option, const QModelIndex& index) const override {
					QStyleOptionViewItem options = option;
					initStyleOption(&options, index);
					QSize size = options.widget->style()->sizeFromContents(QStyle::CT_ItemViewItem, &options, QSize(), options.widget);
					QFontMetrics fm1(options.font);
					QFontMetrics fm2(_tooltipFont);
					size.setHeight(std::max(size.height(), fm1.height() + fm2.height() + 8));
					return size;
				}
			};
			static_cast<QListView*>(_completer->popup())->setUniformItemSizes(true);
			static_cast<QListView*>(_completer->popup())->setItemDelegate(new ItemDelegate());
			_completer->popup()->setIconSize(QSize(44, 32));

			connect(_completer, qOverload<const QModelIndex&>(&QCompleter::activated), actionManager, &ActionManager::onQuickSearchCommandSelected);
			connect(_completer, qOverload<const QModelIndex&>(&QCompleter::activated), this, &QLineEdit::clear);
		}
		void showPopup() {
			if(_completer->popup()->isVisible() == false && text().isEmpty())
				_actionManager->updateActionStates();
			_completer->setCompletionPrefix(text().trimmed());
			_completer->popup()->setCurrentIndex(_completer->completionModel()->index(0,0));
			QRect rect = this->rect();
			rect.setWidth(rect.width() * 2);
			if(layoutDirection() == Qt::RightToLeft)
				rect.setLeft(rect.left() - rect.width() / 2);
			_completer->complete(rect);
		}
	protected:
		ActionManager* _actionManager;
		QCompleter* _completer;
		virtual void keyPressEvent(QKeyEvent* event) override {
			if(_completer->popup()->isVisible()) {
				if(event->key() == Qt::Key_Enter || event->key() == Qt::Key_Return || event->key() == Qt::Key_Escape || event->key() == Qt::Key_Tab) {
					event->ignore();
					return;
				}
			}
			else {
				if(event->key() == Qt::Key_Escape) {
					event->ignore();
					clearFocus();
					return;
				}
			}

			QLineEdit::keyPressEvent(event);

			if(event->key() != Qt::Key_Control && event->key() != Qt::Key_Shift && event->key() != Qt::Key_Meta && event->key() != Qt::Key_Alt) {
				showPopup();
			}
		}		
		virtual void focusInEvent(QFocusEvent* event) override {
			QLineEdit::focusInEvent(event);
			if(event->reason() == Qt::MouseFocusReason || event->reason() == Qt::ShortcutFocusReason || event->reason() == Qt::OtherFocusReason) {
				showPopup();
			}
		}
		virtual void focusOutEvent(QFocusEvent* event) override {
			QLineEdit::focusOutEvent(event);
			clear();
		}
		virtual bool event(QEvent* event) override {
			// This is required to forward key up/down input events to the popup listbox.
			if(event->type() == QEvent::ShortcutOverride) {
				if(_completer->popup() && _completer->popup()->isVisible()) {
					event->accept();
					return true;
				}
			}
			return QLineEdit::event(event);
		}
	};

	// Set up the command quick search field.
	SearchField* commandQuickSearchInputField = new SearchField(this);
	commandQuickSearchInputField->setPlaceholderText(tr("Quick command search (%1)").arg(commandQuickSearchAction->shortcut().toString(QKeySequence::NativeText)));
	commandQuickSearchInputField->setMaximumSize(QSize(260, commandQuickSearchInputField->maximumHeight()));
	commandQuickSearchAction->setDefaultWidget(commandQuickSearchInputField);

	// Set input focus to search field when action shortcut is triggered.  
	connect(commandQuickSearchAction, &QAction::triggered, commandQuickSearchInputField, [commandQuickSearchInputField]() {
		commandQuickSearchInputField->setFocus(Qt::ShortcutFocusReason);
		commandQuickSearchInputField->showPopup();
	});

	addAction(commandQuickSearchAction);
}

// Is called when the user selects a command in the quick search field.
void ActionManager::onQuickSearchCommandSelected(const QModelIndex& index)
{
	QAction* action = index.data(ActionRole).value<QAction*>();
	if(action && action->isEnabled())
		action->trigger();
}

}	// End of namespace