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

class PipelineListModel;	// Defined in PipelineListModel.h

class ModifierAction : public QAction
{
	Q_OBJECT

public:

	/// Constructs an action for a built-in modifier class.
	static ModifierAction* createForClass(ModifierClassPtr clazz);

	/// Constructs an action for a modifier template.
	static ModifierAction* createForTemplate(const QString& templateName);

	/// Constructs an action for a Python modifier script.
	static ModifierAction* createForScript(const QString& fileName, const QDir& directory);

	/// Returns the modifier's category.
	const QString& category() const { return _category; }

	/// Returns the modifier class descriptor if this action represents a built-in modifier. 
	ModifierClassPtr modifierClass() const { return _modifierClass; }

	/// The name of the modifier template if this action represents a saved modifier template.
	const QString& templateName() const { return _templateName; }

	/// The absolute path of the modifier script if this action represents a Python-based modifier function.
	const QString& scriptPath() const { return _scriptPath; }

	/// Updates the actions enabled/disabled state depending on the current data pipeline.
	bool updateState(const PipelineFlowState& input);

private:

	/// The Ovito class descriptor of the modifier subclass. 
	ModifierClassPtr _modifierClass = nullptr;

	/// The modifier's category.
	QString _category;

	/// The name of the modifier template.
	QString _templateName;

	/// The path to the modifier script on disk.
	QString _scriptPath;
};

/**
 * A Qt list model that list all available modifier types that are applicable to the current data pipeline.
 */
class ModifierListModel : public QAbstractListModel
{
	Q_OBJECT

public:

	/// Constructor.
	ModifierListModel(QObject* parent, MainWindow* mainWindow, PipelineListModel* pipelineListModel);

	/// Returns the number of rows in the model.
	virtual int rowCount(const QModelIndex& parent) const override;

	/// Returns the data associated with a list item.
	virtual QVariant data(const QModelIndex& index, int role) const override;

	/// Returns the flags for an item.
	virtual Qt::ItemFlags flags(const QModelIndex& index) const override;

	/// Returns the action that belongs to the given model index.
	ModifierAction* actionFromIndex(int index) const;

	/// Returns the action that belongs to the given model index.
	ModifierAction* actionFromIndex(const QModelIndex& index) const { return actionFromIndex(index.row()); }

	/// Returns the index of the modifier category to which the given list model index belongs.
	int categoryIndexFromListIndex(int index) const;

	/// Returns the list model index where the given modifier category starts.
	int listIndexFromCategoryIndex(int categoryIndex) const;

	/// Returns the category index for the modifier templates.
	int modifierTemplatesCategory() const { return _actionsPerCategory.size() - 2; }

	/// Returns the category index for the Python modifiers.
	int modifierScriptsCategory() const { return _actionsPerCategory.size() - 1; }

	/// Returns whether sorting of available modifiers into categories is enabled.
	bool useCategories() const { return _useCategories; }

	/// Sets whether available modifiers are storted by category instead of name.
	void setUseCategories(bool on);

	/// Returns whether sorting of available modifiers into categories is enabled globally for the application.
	static bool useCategoriesGlobal();

	/// Sets whether available modifiers are storted by category gloablly for the application.
	static void setUseCategoriesGlobal(bool on);

public Q_SLOTS:

	/// Updates the enabled/disabled state of all modifier actions based on the current pipeline.
	void updateActionState();

private Q_SLOTS:

	/// Signal handler that inserts the selected modifier into the current pipeline.
	void insertModifier();

	/// Rebuilds the list of actions for the modifier templates.
	void refreshModifierTemplates();

private:

	/// The complete list of modifier actions, sorted alphabetically.
	std::vector<ModifierAction*> _allActions;

	/// The list of modifier actions, sorted by category.
	std::vector<std::vector<ModifierAction*>> _actionsPerCategory;

	/// The list of modifier categories.
	std::vector<QString> _categoryNames;

	/// The context main window.
	MainWindow* _mainWindow;

	/// The model representing the current data pipeline.
	PipelineListModel* _pipelineListModel;

	/// The list of directories searched for user-defined modifier scripts.
	QVector<QDir> _modifierScriptDirectories;

	/// The font used for category header items.
	QFont _categoryFont;

	/// Colors used for category header items.
	QBrush _categoryBackgroundBrush{Qt::lightGray, Qt::Dense4Pattern};
	QBrush _categoryForegroundBrush{Qt::blue};

	/// Controls the sorting of available modifiers into categories.
	bool _useCategories = useCategoriesGlobal();
};

}	// End of namespace
