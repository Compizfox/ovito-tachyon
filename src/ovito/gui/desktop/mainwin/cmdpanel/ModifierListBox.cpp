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
#include <ovito/core/app/PluginManager.h>
#include <ovito/core/dataset/pipeline/Modifier.h>
#include <ovito/core/dataset/pipeline/PipelineObject.h>
#include <ovito/core/dataset/pipeline/ModifierTemplates.h>
#include <ovito/core/dataset/scene/PipelineSceneNode.h>
#include <ovito/core/dataset/DataSetContainer.h>
#include <ovito/core/dataset/animation/AnimationSettings.h>
#include "ModifierListBox.h"
#include "PipelineListModel.h"

namespace Ovito {

/******************************************************************************
* Initializes the widget.
******************************************************************************/
ModifierListBox::ModifierListBox(QWidget* parent, PipelineListModel* pipelineList) : QComboBox(parent),
		_pipelineList(pipelineList)
{
	setSizeAdjustPolicy(QComboBox::AdjustToContents);

	// A category of modifiers.
	struct ModifierCategory {
		QString name;
		QVector<ModifierClassPtr> modifierClasses;
	};

	QVector<ModifierCategory> modifierCategories;
	ModifierCategory otherCategory;
	otherCategory.name = tr("Others");

	// Retrieve all installed modifier classes.
	for(ModifierClassPtr clazz : PluginManager::instance().metaclassMembers<Modifier>()) {
		// Sort modifiers into categories.
		QString categoryName = clazz->modifierCategory();
		if(categoryName == QStringLiteral("-")) {
			// This modifier requests to be hidden from the user.
			// Do not add it to the list of available modifiers.
		}
		else if(!categoryName.isEmpty()) {
			// Check if category has already been created.
			bool found = false;
			for(auto& category : modifierCategories) {
				if(category.name == categoryName) {
					category.modifierClasses.push_back(clazz);
					found = true;
					break;
				}
			}
			// Create a new category.
			if(!found) {
				ModifierCategory category;
				category.name = categoryName;
				category.modifierClasses.push_back(clazz);
				modifierCategories.push_back(category);
			}
		}
		else {
			// Insert modifiers that don't have category information to the "Other" category.
			otherCategory.modifierClasses.push_back(clazz);
		}
	}

	// Sort category list.
	std::sort(modifierCategories.begin(), modifierCategories.end(), [](const ModifierCategory& a, const ModifierCategory& b) {
		return QString::compare(a.name, b.name, Qt::CaseInsensitive) < 0;
	} );

	// Assign modifiers that haven't been assigned to a category yet to the "Other" category.
	if(!otherCategory.modifierClasses.isEmpty())
		modifierCategories.push_back(otherCategory);

	// Sort modifiers in each category.
	for(auto& category : modifierCategories) {
		std::sort(category.modifierClasses.begin(), category.modifierClasses.end(), [](ModifierClassPtr a, ModifierClassPtr b) {
			return QString::compare(a->displayName(), b->displayName(), Qt::CaseInsensitive) < 0;
		} );
	}

	// Define font colors etc.
	QFont categoryFont = font();
	categoryFont.setBold(true);
	if(categoryFont.pixelSize() < 0)
		categoryFont.setPointSize(categoryFont.pointSize() * 4 / 5);
	else
		categoryFont.setPixelSize(categoryFont.pixelSize() * 4 / 5);
	QBrush categoryBackgroundBrush(Qt::lightGray, Qt::Dense4Pattern);
	QBrush categoryForegroundBrush(Qt::blue);

	// Populate item model.
	_model = new QStandardItemModel(this);

	// Lists starts with the special "Add modification...", which is used just
	// as a label for the combo box.
	QStandardItem* titleItem = new QStandardItem(tr("Add modification..."));
	titleItem->setFlags(Qt::ItemIsEnabled);
	_model->appendRow(titleItem);

	QStandardItem* mruListItem = new QStandardItem(tr("Most recently used modifiers"));
	mruListItem->setFont(categoryFont);
	mruListItem->setBackground(categoryBackgroundBrush);
	mruListItem->setForeground(categoryForegroundBrush);
	mruListItem->setFlags(Qt::ItemIsEnabled);
	mruListItem->setTextAlignment(Qt::AlignCenter);
	_model->appendRow(mruListItem);

	// Create items for all modifiers and the category titles.
	for(const ModifierCategory& category : modifierCategories) {
		if(category.modifierClasses.empty()) continue;

		QStandardItem* categoryItem = new QStandardItem(category.name);
		categoryItem->setFont(categoryFont);
		categoryItem->setBackground(categoryBackgroundBrush);
		categoryItem->setForeground(categoryForegroundBrush);
		categoryItem->setFlags(Qt::ItemIsEnabled);
		categoryItem->setTextAlignment(Qt::AlignCenter);
		_model->appendRow(categoryItem);

		for(ModifierClassPtr descriptor : category.modifierClasses) {
			QStandardItem* modifierItem = new QStandardItem("   " + descriptor->displayName());
			modifierItem->setData(QVariant::fromValue(descriptor), Qt::UserRole);
			_model->appendRow(modifierItem);
			_modifierItems.push_back(modifierItem);
		}
	}

	// Create category for modifier templates.
	QStandardItem* categoryItem = new QStandardItem(tr("Modifier templates"));
	categoryItem->setFont(categoryFont);
	categoryItem->setBackground(categoryBackgroundBrush);
	categoryItem->setForeground(categoryForegroundBrush);
	categoryItem->setFlags(Qt::ItemIsEnabled);
	categoryItem->setTextAlignment(Qt::AlignCenter);
	_model->appendRow(categoryItem);

	// Create category for user-defined modifier scripts.
	categoryItem = new QStandardItem(tr("Python modifiers"));
	categoryItem->setFont(categoryFont);
	categoryItem->setBackground(categoryBackgroundBrush);
	categoryItem->setForeground(categoryForegroundBrush);
	categoryItem->setFlags(Qt::ItemIsEnabled);
	categoryItem->setTextAlignment(Qt::AlignCenter);
	_model->appendRow(categoryItem);

	// Append the "Show all modifiers" item at the end of the list.
	QStandardItem* showAllItem = new QStandardItem(tr("Show all modifiers..."));
	QFont boldFont = font();
	boldFont.setBold(true);
	showAllItem->setFont(boldFont);
	showAllItem->setTextAlignment(Qt::AlignCenter);
	_model->appendRow(showAllItem);

	// Filler item to workaround bug in Qt which doesn't fully show all items in the drop-down menu.
	QStandardItem* fillerItem = new QStandardItem();
	fillerItem->setFlags(Qt::ItemIsEnabled);
	_model->appendRow(fillerItem);

	// Expand list when the "Show all modifiers" entry is selected and update MRU list.
    connect(this, (void (QComboBox::*)(int))&QComboBox::activated, this, [this](int index) {
		if(!showAllModifiers() && index >= count() - 2 && !itemData(index).isValid()) {
			_showAllModifiers = true;
			showPopup();
		}
		else {
			if(itemData(index).isValid()) {
				updateMRUList(itemText(index));
			}
		}
	}, Qt::QueuedConnection);

	// Set up filter model.
	class MyFilterModel : public QSortFilterProxyModel {
	public:
		MyFilterModel(QObject* parent) : QSortFilterProxyModel(parent) {}
	protected:
		virtual bool filterAcceptsRow(int source_row, const QModelIndex& source_parent) const override {
			// Delegate to parent class.
			return static_cast<ModifierListBox*>(parent())->filterAcceptsRow(source_row, source_parent);
		}
		virtual bool lessThan(const QModelIndex& source_left, const QModelIndex& source_right) const override {
			// Delegate to parent class.
			return static_cast<ModifierListBox*>(parent())->filterSortLessThan(source_left, source_right);
		}
	};
	_filterModel = new MyFilterModel(this);
	_filterModel->setDynamicSortFilter(false);
	_filterModel->sort(0);
	_filterModel->setSourceModel(_model);
	setModel(_filterModel);

    connect(this, (void (QComboBox::*)(int))&QComboBox::activated, this, &ModifierListBox::listItemSelected);

	// Add the built-in extension script directory.
	QDir prefixDir(QCoreApplication::applicationDirPath());
	_modifierScriptDirectories.push_back(prefixDir.absolutePath() + QChar('/') + QStringLiteral(OVITO_SCRIPT_EXTENSIONS_RELATIVE_PATH) + QStringLiteral("/modifiers"));

	// Add the user extension script directory.
	_modifierScriptDirectories.push_back(QDir::homePath() + QStringLiteral("/.config/Ovito/scripts/modifiers"));

	for(QDir& dir : _modifierScriptDirectories)
		dir.makeAbsolute();
}

/******************************************************************************
* Filters the full list of modifiers to show only most recently used ones.
******************************************************************************/
bool ModifierListBox::filterAcceptsRow(int source_row, const QModelIndex& source_parent)
{
	if(showAllModifiers()) {
		// Don't show the "Most recently used" entry if all modifier are shown.
		if(source_row == 1)
			return false;
		// Don't show the "Show all modifiers" entry if all are already shown.
		if(source_row >= _model->rowCount(source_parent) - 2)
			return false;
		// Don't show the modifier templates category if there are no templates.
		if(_numModifierTemplates == 0 && source_row == _model->rowCount(source_parent) - 4 - _numModifierScripts)
			return false;
		// Don't show the user-defined modifier scripts category if there are no scripts.
		if(_numModifierScripts == 0 && source_row == _model->rowCount(source_parent) - 3)
			return false;
		return true;
	}
	// Always show the "Add modification..." entry.
	if(source_row == 0)
		return true;

	// Always show the "Most recently used" entry.
	if(source_row == 1)
		return true;

	// Show the "Show all modifiers" entry
	if(source_row >= _model->rowCount(source_parent) - 2)
		return true;

	// Don't show modifier categories.
	if(!_model->index(source_row, 0, source_parent).data(Qt::UserRole).isValid())
		return false;

	// Only show modifiers from MRU list.
	QString modifierName = _model->index(source_row, 0, source_parent).data().toString();
	return _mostRecentlyUsedModifiers.contains(modifierName);
}

/******************************************************************************
* Determines the sort order of the modifier list.
******************************************************************************/
bool ModifierListBox::filterSortLessThan(const QModelIndex& source_left, const QModelIndex& source_right)
{
	if(showAllModifiers() || source_left.row() <= 1 || source_right.row() <= 1 || source_left.row() >= _model->rowCount() - 2 || source_right.row() >= _model->rowCount() - 2) {
		return source_left.row() < source_right.row();
	}
	else {
		return source_left.data().toString().localeAwareCompare(source_right.data().toString()) < 0;
	}
}

/******************************************************************************
* Updates the MRU list after the user has selected a modifier.
******************************************************************************/
void ModifierListBox::updateMRUList(const QString& selectedModifierName)
{
	QSettings settings;
	settings.beginGroup("core/modifier/mru/");
	if(!settings.value("enable_mru", false).toBool())
		return;

	int index = _mostRecentlyUsedModifiers.indexOf(selectedModifierName);
	if(index >= 0) {
		_mostRecentlyUsedModifiers.removeAt(index);
	}
	else if(_mostRecentlyUsedModifiers.size() >= _maxMRUSize) {
		_mostRecentlyUsedModifiers.pop_back();
	}
	_mostRecentlyUsedModifiers.push_front(selectedModifierName);

	// Store MRU list in application settings.
	settings.setValue("list", QVariant::fromValue(_mostRecentlyUsedModifiers));

	// Update list of modifiers shown in the combo box.
	_filterModel->invalidate();
}

/******************************************************************************
* Updates the list box of modifier classes that can be applied to the current selected
* item in the modification list.
******************************************************************************/
void ModifierListBox::updateApplicableModifiersList()
{
	// Always select the "Add modification..." entry by default.
	setCurrentIndex(0);

	// Should we show an MRU list?
	QSettings settings;
	settings.beginGroup("core/modifier/mru/");
	if(settings.value("enable_mru", false).toBool())
		_mostRecentlyUsedModifiers = settings.value("list").toStringList();
	else
		_mostRecentlyUsedModifiers.clear();
	settings.endGroup();

	// Retrieve the input state which a newly inserted modifier would be applied to.
	// This is used to filter the list of available modifiers.
	PipelineListItem* currentItem = _pipelineList->selectedItem();
	while(currentItem && currentItem->parent()) {
		currentItem = currentItem->parent();
	}
	DataSet* dataset = _pipelineList->datasetContainer().currentSet();
	if(!dataset) return;

	PipelineFlowState inputState;
	if(ModifierApplication* modApp = dynamic_object_cast<ModifierApplication>(currentItem->object())) {
		inputState = modApp->evaluateSynchronous(dataset->animationSettings()->time());
	}
	else if(PipelineObject* pipelineObject = dynamic_object_cast<PipelineObject>(currentItem->object())) {
		inputState = pipelineObject->evaluateSynchronous(dataset->animationSettings()->time());
	}
	else if(PipelineSceneNode* node = _pipelineList->selectedNode()) {
		inputState = node->evaluatePipelineSynchronous(false);
	}

	// Update state of combo box items.
	for(QStandardItem* item : _modifierItems) {
		ModifierClassPtr modifierClass = item->data(Qt::UserRole).value<ModifierClassPtr>();
		OVITO_ASSERT(modifierClass);
		item->setEnabled(inputState.data() && modifierClass->isApplicableTo(*inputState.data()));
	}

	// Load custom modifier templates.
	ModifierTemplates modifierTemplates;
	int numCustom = 0;
	int currentRowIndex = _model->rowCount() - 3 - _numModifierTemplates - _numModifierScripts;
	for(const QString& name : modifierTemplates.templateList()) {
		QStandardItem* modifierItem;
		if(numCustom < _numModifierTemplates) {
			modifierItem = _model->item(currentRowIndex);
		}
		else {
			modifierItem = new QStandardItem();
			_model->insertRow(currentRowIndex, modifierItem);
		}
		modifierItem->setText(QStringLiteral("   ") + name);
		modifierItem->setData(QVariant::fromValue(name), Qt::UserRole);
		numCustom++;
		currentRowIndex++;
	}
	// Remove excess list items.
	if(numCustom < _numModifierTemplates)
		_model->removeRows(currentRowIndex, _numModifierTemplates - numCustom);
	_numModifierTemplates = numCustom;
	currentRowIndex++;

	// Load user-defined Python script modifiers.
	int numScripts = 0;
	for(const QDir& scriptsDirectory : _modifierScriptDirectories) {
		QStringList scriptFiles = scriptsDirectory.entryList(QStringList() << QStringLiteral("*.py"), QDir::Files, QDir::Name);
		for(const QString& fileName : scriptFiles) {
			QStandardItem* modifierItem;
			if(numScripts < _numModifierScripts) {
				modifierItem = _model->item(currentRowIndex);
			}
			else {
				modifierItem = new QStandardItem();
				_model->insertRow(currentRowIndex, modifierItem);
			}
			modifierItem->setText(QStringLiteral("   ") + fileName.chopped(3));
			modifierItem->setData(QVariant::fromValue(scriptsDirectory.filePath(fileName)), Qt::UserRole);
			numScripts++;
			currentRowIndex++;	
		}
	}
	// Remove excess list items.
	if(numScripts < _numModifierScripts)
		_model->removeRows(currentRowIndex, _numModifierScripts - numScripts);
	_numModifierScripts = numScripts;
}

/******************************************************************************
* Is called when the user has selected an item in the modifier class list.
******************************************************************************/
void ModifierListBox::listItemSelected(int index)
{
	QModelIndex sourceIndex = _filterModel->mapToSource(_filterModel->index(index, 0));

	DataSet* dataset = _pipelineList->datasetContainer().currentSet();

	QVariant itemData = sourceIndex.data(Qt::UserRole);
	if(ModifierClassPtr modifierClass = itemData.value<ModifierClassPtr>()) {
		UndoableTransaction::handleExceptions(dataset->undoStack(), tr("Apply modifier"), [&]() {
			// Create an instance of the modifier class.
			OORef<Modifier> modifier = static_object_cast<Modifier>(modifierClass->createInstance(dataset));
			OVITO_CHECK_OBJECT_POINTER(modifier);
			// Load user-defined default parameter values.
			modifier->loadUserDefaults();
			// Apply it.
			Q_EMIT applyModifiers({modifier});
		});
	}
	else if(sourceIndex.row() >= _model->rowCount() - 3 - _numModifierScripts - _numModifierTemplates && sourceIndex.row() < _model->rowCount() - 3 - _numModifierScripts) {
		QString templateName = itemData.toString();
		OVITO_ASSERT(!templateName.isEmpty());

		UndoableTransaction::handleExceptions(dataset->undoStack(), tr("Insert modifier template"), [&]() {
			// Load modifier template from the store.
			ModifierTemplates modifierTemplates;
			QVector<OORef<Modifier>> modifierSet = modifierTemplates.instantiateTemplate(templateName, dataset);
			Q_EMIT applyModifiers(modifierSet);
		});
	}
	else if(sourceIndex.row() >= _model->rowCount() - 2 - _numModifierScripts && sourceIndex.row() < _model->rowCount() - 2) {
		QString scriptPath = itemData.toString();
		OVITO_ASSERT(!scriptPath.isEmpty());

		UndoableTransaction::handleExceptions(dataset->undoStack(), tr("Insert Python modifier"), [&]() {

			// Load source code from template file.
			QFile file(scriptPath);
			if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
				throw Exception(tr("Failed to open Python file '%1' for reading: %2").arg(scriptPath).arg(file.errorString()));
			QString scriptCode = tr("# This is a copy of the template file '%1'.\n# Feel free to modify the code below as needed.\n\n").arg(QDir::toNativeSeparators(scriptPath)) + QString::fromUtf8(file.readAll());

			// Get the PythonScriptModifier modifier class.
			if(OvitoClassPtr clazz = PluginManager::instance().findClass({}, QStringLiteral("PythonScriptModifier"))) {
				if(!clazz->isAbstract() && clazz->isDerivedFrom(Modifier::OOClass())) {
					const Modifier::OOMetaClass* modifierClass = static_cast<const Modifier::OOMetaClass*>(clazz);

					// Instantiate the PythonScriptModifier class.
					OORef<Modifier> modifier = static_object_cast<Modifier>(modifierClass->createInstance(dataset));
					OVITO_CHECK_OBJECT_POINTER(modifier);
					modifier->setTitle(sourceIndex.data(Qt::DisplayRole).toString().trimmed());
					modifier->loadUserDefaults();

					// Assign the script code.
					const PropertyFieldDescriptor* scriptPropertyField = modifierClass->findPropertyField("script");
					OVITO_ASSERT(scriptPropertyField);
					modifier->setPropertyFieldValue(*scriptPropertyField, QVariant::fromValue(scriptCode));

					// Insert modifier into data pipeline.
					Q_EMIT applyModifiers({modifier});
				}
			}
		});
	}

	setCurrentIndex(0);
}

}	// End of namespace
