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
#include <ovito/gui/desktop/mainwin/MainWindow.h>
#include <ovito/gui/desktop/actions/ActionManager.h>
#include <ovito/gui/desktop/dataset/GuiDataSetContainer.h>
#include <ovito/core/app/PluginManager.h>
#include <ovito/core/dataset/pipeline/Modifier.h>
#include <ovito/core/dataset/DataSetContainer.h>
#include <ovito/core/dataset/pipeline/ModifierTemplates.h>
#include "ModifierListModel.h"
#include "PipelineListModel.h"

namespace Ovito {

/******************************************************************************
* Constructs an action for a built-in modifier class.
******************************************************************************/
ModifierAction* ModifierAction::createForClass(ModifierClassPtr clazz)
{
	ModifierAction* action = new ModifierAction();
	action->_modifierClass = clazz;
	action->_category = clazz->modifierCategory();

	// Generate a unique identifier for the action:
	action->setObjectName(QStringLiteral("InsertModifier.%1.%2").arg(clazz->pluginId(), clazz->name()));

	// Set the action's UI display name.
	action->setText(clazz->displayName());

	// Give the modifier a status bar text.
	QString description = clazz->descriptionString();
	action->setStatusTip(!description.isEmpty() ? std::move(description) : tr("Insert this modifier into the data pipeline."));

	// Give the action an icon.
	static QIcon icon(":/gui/actions/modify/modifier_action_icon.svg");
	action->setIcon(icon);

	// Modifiers without a category are moved into the "Other" category.
	if(action->_category.isEmpty())
		action->_category = tr("Other");
	
	return action;
}

/******************************************************************************
* Constructs an action for a modifier template.
******************************************************************************/
ModifierAction* ModifierAction::createForTemplate(const QString& templateName)
{
	ModifierAction* action = new ModifierAction();
	action->_templateName = templateName;

	// Generate a unique identifier for the action:
	action->setObjectName(QStringLiteral("InsertModifierTemplate.%1").arg(templateName));

	// Set the action's UI display name.
	action->setText(templateName);

	// Give the modifier a status bar text.
	action->setStatusTip(tr("Insert this modifier template into the data pipeline."));

	// Give the action an icon.
	static QIcon icon(":/gui/actions/modify/modifier_action_icon.svg");
	action->setIcon(icon);
	
	return action;
}

/******************************************************************************
* Constructs an action for a Python modifier script.
******************************************************************************/
ModifierAction* ModifierAction::createForScript(const QString& fileName, const QDir& directory)
{
	ModifierAction* action = new ModifierAction();
	action->_scriptPath = directory.filePath(fileName);

	// Generate a unique identifier for the action:
	action->setObjectName(QStringLiteral("InsertModifierScript.%1").arg(action->_scriptPath));

	// Set the action's UI display name. Chop of ".py" extension of filename.
#if QT_VERSION >= QT_VERSION_CHECK(5, 10, 0)
	action->setText(fileName.chopped(3));
#else
	action->setText(fileName.left(fileName.size() - 3));
#endif

	// Give the modifier a status bar text.
	action->setStatusTip(tr("Insert this Python modifier into the data pipeline."));

	// Give the action an icon.
	static QIcon icon(":/gui/actions/modify/modifier_action_icon.svg");
	action->setIcon(icon);
	
	return action;
}

/******************************************************************************
* Updates the actions enabled/disabled state depending on the current data pipeline.
******************************************************************************/
bool ModifierAction::updateState(const PipelineFlowState& input)
{
	bool enable = input.data() && (!modifierClass() || modifierClass()->isApplicableTo(*input.data()));
	if(isEnabled() != enable) {
		setEnabled(enable);
		return true;
	}
	return false;
}

/******************************************************************************
* Constructor.
******************************************************************************/
ModifierListModel::ModifierListModel(QObject* parent, MainWindow* mainWindow, PipelineListModel* pipelineListModel) : QAbstractListModel(parent), _mainWindow(mainWindow), _pipelineListModel(pipelineListModel)
{
	// Update the state of this model's actions whenther the ActionManager requests it.
	connect(_mainWindow->actionManager(), &ActionManager::actionUpdateRequested, this, &ModifierListModel::updateActionState);

	// Enumerate all built-in modifier classes.
	std::vector<ModifierAction*> actions;
	for(ModifierClassPtr clazz : PluginManager::instance().metaclassMembers<Modifier>()) {

		// Skip modifiers that want to be hidden from the user.
		// Do not add it to the list of available modifiers.
		if(clazz->modifierCategory() == QStringLiteral("-"))
			continue;

		// Create action for the modifier class.
		ModifierAction* action = ModifierAction::createForClass(clazz);
		actions.push_back(action);

		// Register it with the global ActionManager.
		_mainWindow->actionManager()->addAction(action);
		OVITO_ASSERT(action->parent() == _mainWindow->actionManager());

		// Handle the insertion action.
		connect(action, &QAction::triggered, this, &ModifierListModel::insertModifier);
	}

	// Order actions list by category name.
	std::sort(actions.begin(), actions.end(), [](ModifierAction* a, ModifierAction* b) { return QString::localeAwareCompare(a->category(), b->category()) < 0; });

	// Sort actions into categories.
	for(ModifierAction* action : actions) {
		if(_categoryNames.empty() || _categoryNames.back() != action->category()) {
			_categoryNames.push_back(action->category());
			_actions.emplace_back();
		}
		_actions.back().push_back(action);
	}

	// Sort actions by name within each category.
	for(std::vector<ModifierAction*>& actions : _actions)
		std::sort(actions.begin(), actions.end(), [](ModifierAction* a, ModifierAction* b) { return QString::localeAwareCompare(a->text(), b->text()) < 0; });

	// Create category for modifier templates.
	_categoryNames.push_back(tr("Modifier templates"));
	_actions.emplace_back();
	for(const QString& templateName : ModifierTemplates::get()->templateList()) {
		// Create action for the modifier template.
		ModifierAction* action = ModifierAction::createForTemplate(templateName);
		_actions.back().push_back(action);

		// Register it with the global ActionManager.
		_mainWindow->actionManager()->addAction(action);
		OVITO_ASSERT(action->parent() == _mainWindow->actionManager());

		// Handle the action.
		connect(action, &QAction::triggered, this, &ModifierListModel::insertModifier);
	}

	// Listen for changes to the underlying modifier template list.
	connect(ModifierTemplates::get(), &QAbstractItemModel::rowsInserted, this, &ModifierListModel::refreshModifierTemplates);
	connect(ModifierTemplates::get(), &QAbstractItemModel::rowsRemoved, this, &ModifierListModel::refreshModifierTemplates);
	connect(ModifierTemplates::get(), &QAbstractItemModel::modelReset, this, &ModifierListModel::refreshModifierTemplates);
	connect(ModifierTemplates::get(), &QAbstractItemModel::dataChanged, this, &ModifierListModel::refreshModifierTemplates);

	// Add the built-in extension script directory.
	QDir prefixDir(QCoreApplication::applicationDirPath());
	_modifierScriptDirectories.push_back(prefixDir.absolutePath() + QChar('/') + QStringLiteral(OVITO_SCRIPT_EXTENSIONS_RELATIVE_PATH) + QStringLiteral("/modifiers"));
	// Add the user extension script directory.
	_modifierScriptDirectories.push_back(QDir::homePath() + QStringLiteral("/.config/Ovito/scripts/modifiers"));
	for(QDir& dir : _modifierScriptDirectories)
		dir.makeAbsolute();

	// Create category for script modifiers.
#ifndef OVITO_BUILD_BASIC
	_categoryNames.push_back(tr("Python modifiers"));
#else
	_categoryNames.push_back(tr("Python modifiers (Pro)"));
#endif
	_actions.emplace_back();

	// Load user-defined Python script modifiers.
	for(const QDir& scriptsDirectory : _modifierScriptDirectories) {
		QStringList scriptFiles = scriptsDirectory.entryList(QStringList() << QStringLiteral("*.py"), QDir::Files, QDir::Name);
		for(const QString& fileName : scriptFiles) {

			// Create action for the modifier script.
			ModifierAction* action = ModifierAction::createForScript(fileName, scriptsDirectory);
			_actions.back().push_back(action);

			// Register it with the global ActionManager.
			_mainWindow->actionManager()->addAction(action);
			OVITO_ASSERT(action->parent() == _mainWindow->actionManager());

			// Handle the action.
			connect(action, &QAction::triggered, this, &ModifierListModel::insertModifier);
		}
	}

	// Define font, colors, etc.
	_categoryFont = QGuiApplication::font();
	_categoryFont.setBold(true);
#ifndef Q_OS_WIN
	if(_categoryFont.pixelSize() < 0)
		_categoryFont.setPointSize(_categoryFont.pointSize() * 4 / 5);
	else
		_categoryFont.setPixelSize(_categoryFont.pixelSize() * 4 / 5);
#endif
}

/******************************************************************************
* Returns the action that belongs to the given model index.
******************************************************************************/
ModifierAction* ModifierListModel::actionFromIndex(int index) const
{
	if(index == 0) return nullptr;
	index--;

	for(const auto& categoryActions : _actions) {
		if(!categoryActions.empty()) {
			if(index == 0) return nullptr;
			index--;
			if(index < categoryActions.size())
				return categoryActions[index];
			index -= categoryActions.size();
		}
	}

	return nullptr;
}

/******************************************************************************
* Returns the index of the modifier category to which the given list model index belongs.
******************************************************************************/
int ModifierListModel::categoryIndexFromListIndex(int index) const
{
	if(index == 0) return -1;
	index--;

	int categoryIndex = 0;
	for(const auto& categoryActions : _actions) {
		if(index == 0)
			return categoryIndex;
		if(!categoryActions.empty())
			index -= categoryActions.size() + 1;
		categoryIndex++;
	}

	return -1;
}

/******************************************************************************
* Returns the list model index where the given modifier category starts.
******************************************************************************/
int ModifierListModel::listIndexFromCategoryIndex(int categoryIndex) const
{
	int index = 1;

	for(const auto& categoryActions : _actions) {
		if(categoryIndex == 0)
			return index;
		if(!categoryActions.empty())
			index += categoryActions.size() + 1;
		categoryIndex--;
	}

	OVITO_ASSERT(false);
	return -1;
}

/******************************************************************************
* Returns the number of rows in the model.
******************************************************************************/
int ModifierListModel::rowCount(const QModelIndex& parent) const 
{
	int sum = 1; // First entry is the "Add modification..." item.

	for(const auto& categoryActions : _actions)
		if(!categoryActions.empty())
			sum += categoryActions.size() + 1; 	// Take into account category header.

	return sum;
}

/******************************************************************************
* Returns the data associated with a list item.
******************************************************************************/
QVariant ModifierListModel::data(const QModelIndex& index, int role) const
{
	if(role == Qt::DisplayRole) {
		if(ModifierAction* action = actionFromIndex(index)) {
			return action->text();
		}
		else {
			int categoryIndex = categoryIndexFromListIndex(index.row());
			if(categoryIndex < 0)
				return tr("Add modification...");
			else
				return _categoryNames[categoryIndex];
		}
	}
	else if(role == Qt::FontRole) {
		if(categoryIndexFromListIndex(index.row()) >= 0)
			return _categoryFont;
	}
	else if(role == Qt::ForegroundRole) {
		if(categoryIndexFromListIndex(index.row()) >= 0)
			return _categoryForegroundBrush;
	}
	else if(role == Qt::BackgroundRole) {
		if(categoryIndexFromListIndex(index.row()) >= 0)
			return _categoryBackgroundBrush;
	}
	else if(role == Qt::TextAlignmentRole) {
		if(categoryIndexFromListIndex(index.row()) >= 0)
			return Qt::AlignCenter;
	}
	return {};
}

/******************************************************************************
* Returns the flags for an item.
******************************************************************************/
Qt::ItemFlags ModifierListModel::flags(const QModelIndex& index) const
{
	if(categoryIndexFromListIndex(index.row()) >= 0)
		return Qt::ItemIsEnabled;
	else if(ModifierAction* action = actionFromIndex(index))
		return action->isEnabled() ? (Qt::ItemIsEnabled | Qt::ItemIsSelectable) : Qt::NoItemFlags;

	return QAbstractListModel::flags(index);
}

/******************************************************************************
* Signal handler that inserts the selected modifier into the current pipeline.
******************************************************************************/
void ModifierListModel::insertModifier()
{
	// Get the action that emitted the signal.
	ModifierAction* action = qobject_cast<ModifierAction*>(sender());
	OVITO_ASSERT(action);

	// Get the current dataset.
	DataSet* dataset = _mainWindow->datasetContainer().currentSet();

	// Instantiate the new modifier(s) and insert them into the pipeline.
	UndoableTransaction::handleExceptions(dataset->undoStack(), tr("Insert modifier"), [&]() {

		if(action->modifierClass()) {
			// Create an instance of the modifier.
			OORef<Modifier> modifier = static_object_cast<Modifier>(action->modifierClass()->createInstance(dataset));
			// Initialize parameters to user-defined defaults.
			modifier->loadUserDefaults();
			// Insert modifier into the data pipeline.
			_pipelineListModel->applyModifiers({modifier});
		}
		else if(!action->templateName().isEmpty()) {
			// Load modifier template from the store.
			QVector<OORef<Modifier>> modifierSet = ModifierTemplates::get()->instantiateTemplate(action->templateName(), dataset);
			// Insert modifier(s) into the data pipeline.
			_pipelineListModel->applyModifiers(modifierSet);
		}
		else if(!action->scriptPath().isEmpty()) {

			// Load source code from template file.
			QFile file(action->scriptPath());
			if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
				throw Exception(tr("Failed to open Python file '%1' for reading: %2").arg(action->scriptPath()).arg(file.errorString()));
			QString scriptCode = tr("# This is a copy of the template file '%1'.\n# Feel free to modify the code below as needed.\n\n").arg(QDir::toNativeSeparators(action->scriptPath())) + QString::fromUtf8(file.readAll());

			// Get the PythonScriptModifier modifier class.
			if(OvitoClassPtr clazz = PluginManager::instance().findClass({}, QStringLiteral("PythonScriptModifier"))) {
				if(!clazz->isAbstract() && clazz->isDerivedFrom(Modifier::OOClass())) {
					const Modifier::OOMetaClass* modifierClass = static_cast<const Modifier::OOMetaClass*>(clazz);

					// Instantiate the PythonScriptModifier class.
					OORef<Modifier> modifier = static_object_cast<Modifier>(modifierClass->createInstance(dataset));
					OVITO_CHECK_OBJECT_POINTER(modifier);
					modifier->setTitle(action->text());
					modifier->loadUserDefaults();

					// Assign the script code.
					const PropertyFieldDescriptor* scriptPropertyField = modifierClass->findPropertyField("script");
					OVITO_ASSERT(scriptPropertyField);
					modifier->setPropertyFieldValue(*scriptPropertyField, QVariant::fromValue(scriptCode));

					// Insert modifier(s) into the data pipeline.
					_pipelineListModel->applyModifiers({modifier});
				}
			}
		}

	});
}

/******************************************************************************
* Rebuilds the list of actions for the modifier templates.
******************************************************************************/
void ModifierListModel::refreshModifierTemplates()
{
	std::vector<ModifierAction*>& actions = _actions[modifierTemplatesCategory()];

	// Discard old list of actions.
	int startIndex = listIndexFromCategoryIndex(modifierTemplatesCategory());
	if(!actions.empty()) {
		beginRemoveRows(QModelIndex(), startIndex, startIndex + actions.size());
		for(ModifierAction* action : actions)
			_mainWindow->actionManager()->deleteAction(action);
		actions.clear();
		endRemoveRows();
	}

	// Create new actions for the modifier templates.
	int count = ModifierTemplates::get()->templateList().size();
	if(count != 0) {
		beginInsertRows(QModelIndex(), startIndex, startIndex + count);
		for(const QString& templateName : ModifierTemplates::get()->templateList()) {
			// Create action for the modifier template.
			ModifierAction* action = ModifierAction::createForTemplate(templateName);
			actions.push_back(action);

			// Register it with the global ActionManager.
			_mainWindow->actionManager()->addAction(action);
			OVITO_ASSERT(action->parent() == _mainWindow->actionManager());

			// Handle the action.
			connect(action, &QAction::triggered, this, &ModifierListModel::insertModifier);
		}
		endInsertRows();
	}
}

/******************************************************************************
* Updates the enabled/disabled state of all modifier actions based on the current pipeline.
******************************************************************************/
void ModifierListModel::updateActionState()
{
	// Retrieve the input pipeline state, which a newly inserted modifier would be applied to.
	// This is used to determine which modifiers are applicable.
	PipelineFlowState inputState;

	// Get the selected item in the pipeline editor.
	PipelineListItem* currentItem = _pipelineListModel->selectedItem();
	while(currentItem && currentItem->parent()) {
		currentItem = currentItem->parent();
	}

	// Evaluate pipeline at the selected stage.
	if(currentItem) {
		if(DataSet* dataset = _pipelineListModel->datasetContainer().currentSet()) {
			if(ModifierApplication* modApp = dynamic_object_cast<ModifierApplication>(currentItem->object())) {
				inputState = modApp->evaluateSynchronous(dataset->animationSettings()->time());
			}
			else if(PipelineObject* pipelineObject = dynamic_object_cast<PipelineObject>(currentItem->object())) {
				inputState = pipelineObject->evaluateSynchronous(dataset->animationSettings()->time());
			}
			else if(PipelineSceneNode* pipeline = _pipelineListModel->selectedNode()) {
				inputState = pipeline->evaluatePipelineSynchronous(false);
			}
		}
	}

	// Update the actions.
	int row = 1;
	for(const auto& categoryActions : _actions) {
		if(!categoryActions.empty()) 
			row++;
		for(ModifierAction* action : categoryActions) {
			if(action->updateState(inputState))
				Q_EMIT dataChanged(index(row), index(row));
			row++;
		}
	}
}

}	// End of namespace
