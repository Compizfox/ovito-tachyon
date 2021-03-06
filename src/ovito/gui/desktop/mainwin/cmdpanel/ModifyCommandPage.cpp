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
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include <ovito/core/dataset/pipeline/Modifier.h>
#include <ovito/core/dataset/pipeline/ModifierTemplates.h>
#include <ovito/core/dataset/scene/PipelineSceneNode.h>
#include <ovito/core/dataset/scene/SelectionSet.h>
#include <ovito/core/viewport/ViewportConfiguration.h>
#include <ovito/core/dataset/UndoStack.h>
#include <ovito/core/dataset/DataSetContainer.h>
#include <ovito/core/app/Application.h>
#include <ovito/core/oo/CloneHelper.h>
#include <ovito/gui/desktop/actions/ActionManager.h>
#include <ovito/gui/desktop/app/GuiApplication.h>
#include <ovito/gui/desktop/mainwin/MainWindow.h>
#include <ovito/gui/desktop/dialogs/ModifierTemplatesPage.h>
#include <ovito/gui/desktop/widgets/selection/SceneNodeSelectionBox.h>
#include "ModifyCommandPage.h"
#include "PipelineListModel.h"
#include "ModifierListModel.h"

#include <QtNetwork>

namespace Ovito {

/******************************************************************************
* Initializes the modify page.
******************************************************************************/
ModifyCommandPage::ModifyCommandPage(MainWindow* mainWindow, QWidget* parent) : QWidget(parent),
		_datasetContainer(mainWindow->datasetContainer()), _actionManager(mainWindow->actionManager())
{
	QGridLayout* layout = new QGridLayout(this);
	layout->setContentsMargins(2,2,2,2);
	layout->setSpacing(4);
	layout->setColumnStretch(0,1);

	SceneNodeSelectionBox* nodeSelBox = new SceneNodeSelectionBox(_datasetContainer, this);
	layout->addWidget(nodeSelBox, 0, 0, 1, 1);
	QMenu* pipelineMenu = new QMenu(this);
	pipelineMenu->addAction(_actionManager->getAction(ACTION_EDIT_RENAME_PIPELINE));
	pipelineMenu->addAction(_actionManager->getAction(ACTION_EDIT_CLONE_PIPELINE));
	pipelineMenu->addSeparator();
	pipelineMenu->addAction(_actionManager->getAction(ACTION_EDIT_DELETE));
	pipelineMenu->addSeparator();

	// Set up the 'Precompute all frames' menu action.
	QAction* precomputeFramesAction = pipelineMenu->addAction(QIcon(":/gui/actions/file/cache_pipeline_output.svg"), tr("Precompute all frames"));
	precomputeFramesAction->setCheckable(true);
	connect(pipelineMenu, &QMenu::aboutToShow, this, [this,precomputeFramesAction]() {
		if(_datasetContainer.currentSet() && _datasetContainer.currentSet()->selection()->nodes().empty() == false) {
			if(PipelineSceneNode* pipeline = dynamic_object_cast<PipelineSceneNode>(_datasetContainer.currentSet()->selection()->nodes().front())) {
				precomputeFramesAction->setChecked(pipeline->pipelineTrajectoryCachingEnabled());
				precomputeFramesAction->setEnabled(true);
				return;
			}
		}
		precomputeFramesAction->setChecked(false);
		precomputeFramesAction->setEnabled(false);
	});
	connect(precomputeFramesAction, &QAction::triggered, this, [this,precomputeFramesAction]() {
		if(_datasetContainer.currentSet() && _datasetContainer.currentSet()->selection()->nodes().empty() == false) {
			if(PipelineSceneNode* pipeline = dynamic_object_cast<PipelineSceneNode>(_datasetContainer.currentSet()->selection()->nodes().front())) {
				pipeline->setPipelineTrajectoryCachingEnabled(precomputeFramesAction->isChecked());
			}
		}
	});

	QToolButton* pipelineMenuButton = new QToolButton(this);
	pipelineMenuButton->setStyleSheet(
		"QToolButton { padding: 0px; margin: 0px; border: none; background-color: transparent; } "
		"QToolButton::menu-indicator { image: none; } ");
	pipelineMenuButton->setPopupMode(QToolButton::InstantPopup);
	pipelineMenuButton->setIcon(QIcon(":/gui/actions/edit/pipeline_menu.svg"));
	pipelineMenuButton->setMenu(pipelineMenu);
	pipelineMenuButton->setEnabled(nodeSelBox->isEnabled());
	pipelineMenuButton->setToolTip(tr("Pipeline menu"));
	layout->addWidget(pipelineMenuButton, 0, 1, 1, 1);
	connect(nodeSelBox, &SceneNodeSelectionBox::enabledChanged, pipelineMenuButton, &QToolButton::setEnabled);

	_pipelineListModel = new PipelineListModel(_datasetContainer, this);
	class ModifierListBox : public QComboBox {
	public:
		using QComboBox::QComboBox;
		virtual void showPopup() override {
			static_cast<ModifierListModel*>(model())->updateActionState();
			QComboBox::showPopup();
		}
	};
	_modifierSelector = new ModifierListBox(this);
    layout->addWidget(_modifierSelector, 1, 0, 1, 1);
	_modifierSelector->setSizeAdjustPolicy(QComboBox::AdjustToContents);
	_modifierSelector->setModel(new ModifierListModel(this, mainWindow, _pipelineListModel));
	_modifierSelector->setMaxVisibleItems(0xFFFF);
    connect(_modifierSelector, qOverload<int>(&QComboBox::activated), this, [this](int index) {
		QComboBox* selector = static_cast<QComboBox*>(sender());
		if(QAction* action = static_cast<ModifierListModel*>(selector->model())->actionFromIndex(index))
			action->trigger();
		selector->setCurrentIndex(0);
	});

	class PipelineListView : public QListView {
	public:
		PipelineListView(QWidget* parent) : QListView(parent) {}
		virtual QSize sizeHint() const override { return QSize(256, 260); }
	protected:
		virtual bool edit(const QModelIndex& index, QAbstractItemView::EditTrigger trigger, QEvent* event) override {
			// Avoid triggering edit mode when user clicks the check box next to a list item.
			if(trigger == QAbstractItemView::SelectedClicked && event->type() == QEvent::MouseButtonRelease) {
				QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);
				if(mouseEvent->pos().x() < visualRect(index).left() + 50)
					trigger = QAbstractItemView::NoEditTriggers;
			}
			return QListView::edit(index, trigger, event);
		}
	};

	QSplitter* splitter = new QSplitter(Qt::Vertical);
	splitter->setChildrenCollapsible(false);

	QWidget* upperContainer = new QWidget();
	splitter->addWidget(upperContainer);
	QHBoxLayout* subLayout = new QHBoxLayout(upperContainer);
	subLayout->setContentsMargins(0,0,0,0);
	subLayout->setSpacing(2);

	_pipelineWidget = new PipelineListView(upperContainer);
	_pipelineWidget->setDragDropMode(QAbstractItemView::InternalMove);
	_pipelineWidget->setDragEnabled(true);
	_pipelineWidget->setAcceptDrops(true);
	_pipelineWidget->setDragDropOverwriteMode(false);
	_pipelineWidget->setDropIndicatorShown(true);
	_pipelineWidget->setEditTriggers(QAbstractItemView::SelectedClicked);
	_pipelineWidget->setModel(_pipelineListModel);
	_pipelineWidget->setSelectionModel(_pipelineListModel->selectionModel());
	subLayout->addWidget(_pipelineWidget);

	// Listen to selection changes in the pipeline editor list widget.
	connect(_pipelineListModel, &PipelineListModel::selectedItemChanged, this, &ModifyCommandPage::onSelectedItemChanged);

	// Double-click on a modifier or visual element toggles the enabled state of the element.
	connect(_pipelineWidget, &PipelineListView::doubleClicked, this, &ModifyCommandPage::onModifierStackDoubleClicked);

	QToolBar* editToolbar = new QToolBar(this);
	editToolbar->setOrientation(Qt::Vertical);
#ifndef Q_OS_MACX
	editToolbar->setStyleSheet("QToolBar { padding: 0px; margin: 0px; border: 0px none black; }");
#endif
	subLayout->addWidget(editToolbar);

	QAction* deleteModifierAction = _actionManager->createCommandAction(ACTION_MODIFIER_DELETE, tr("Delete Modifier"), ":/gui/actions/modify/delete_modifier.bw.svg", tr("Delete the selected modifier from the pipeline."));
	connect(deleteModifierAction, &QAction::triggered, this, &ModifyCommandPage::onDeleteModifier);
	editToolbar->addAction(deleteModifierAction);

	editToolbar->addSeparator();

	QAction* moveModifierUpAction = _actionManager->createCommandAction(ACTION_MODIFIER_MOVE_UP, tr("Move Modifier Up"), ":/gui/actions/modify/modifier_move_up.bw.svg", tr("Move the selected modifier up in the pipeline."));
	connect(moveModifierUpAction, &QAction::triggered, this, &ModifyCommandPage::onModifierMoveUp);
	editToolbar->addAction(moveModifierUpAction);
	QAction* moveModifierDownAction = mainWindow->actionManager()->createCommandAction(ACTION_MODIFIER_MOVE_DOWN, tr("Move Modifier Down"), ":/gui/actions/modify/modifier_move_down.bw.svg", tr("Move the selected modifier down in the pipeline."));
	connect(moveModifierDownAction, &QAction::triggered, this, &ModifyCommandPage::onModifierMoveDown);
	editToolbar->addAction(moveModifierDownAction);

	QAction* toggleModifierStateAction = _actionManager->createCommandAction(ACTION_MODIFIER_TOGGLE_STATE, tr("Enable/Disable Modifier"), {}, tr("Turn the selected modifier in the pipeline on or off."));
	toggleModifierStateAction->setCheckable(true);
	QIcon toggleStateActionIcon(QString(":/gui/actions/modify/modifier_enabled_large.png"));
	toggleStateActionIcon.addFile(QString(":/gui/actions/modify/modifier_disabled_large.png"), QSize(), QIcon::Normal, QIcon::On);
	toggleModifierStateAction->setIcon(toggleStateActionIcon);
	connect(toggleModifierStateAction, &QAction::triggered, this, &ModifyCommandPage::onModifierToggleState);

	editToolbar->addSeparator();
	QAction* makeElementIndependentAction = _actionManager->createCommandAction(ACTION_PIPELINE_MAKE_INDEPENDENT, tr("Replace With Independent Copy"), ":/gui/actions/modify/make_element_independent.bw.svg", tr("Duplicate an entry that is shared by multiple pipelines."));
	connect(makeElementIndependentAction, &QAction::triggered, this, &ModifyCommandPage::onMakeElementIndependent);
	editToolbar->addAction(makeElementIndependentAction);

	QAction* manageModifierTemplatesAction = _actionManager->createCommandAction(ACTION_MODIFIER_MANAGE_TEMPLATES, tr("Manage Modifier Templates..."), ":/gui/actions/modify/modifier_save_preset.bw.svg", tr("Open the dialog that lets you manage the saved modifier templates."));
	connect(manageModifierTemplatesAction, &QAction::triggered, [mainWindow]() {
		ApplicationSettingsDialog dlg(mainWindow, &ModifierTemplatesPage::OOClass());
		dlg.exec();
	});
	editToolbar->addAction(manageModifierTemplatesAction);

	layout->addWidget(splitter, 2, 0, 1, 2);
	layout->setRowStretch(2, 1);

	// Create the properties panel.
	_propertiesPanel = new PropertiesPanel(nullptr, mainWindow);
	_propertiesPanel->setFrameStyle(QFrame::NoFrame | QFrame::Plain);
	splitter->addWidget(_propertiesPanel);
	splitter->setStretchFactor(1,1);

	connect(&_datasetContainer, &DataSetContainer::selectionChangeComplete, this, &ModifyCommandPage::onSelectionChangeComplete);
	updateActions(nullptr);

	// Create About panel.
	createAboutPanel();
}

/******************************************************************************
* This is called after all changes to the selection set have been completed.
******************************************************************************/
void ModifyCommandPage::onSelectionChangeComplete(SelectionSet* newSelection)
{
	// Make sure we get informed about any future changes of the selection set.
	_pipelineListModel->refreshList();
}

/******************************************************************************
* Is called when a new modification list item has been selected, or if the currently
* selected item has changed.
******************************************************************************/
void ModifyCommandPage::onSelectedItemChanged()
{
	PipelineListItem* currentItem = pipelineListModel()->selectedItem();
	RefTarget* editObject = nullptr;

	if(currentItem != nullptr) {
		_aboutRollout->hide();
		editObject = currentItem->object();
		if(currentItem->isSubObject())
			pipelineListModel()->setNextSubObjectTitleToSelect(currentItem->title());
	}

	if(editObject != _propertiesPanel->editObject()) {
		_propertiesPanel->setEditObject(editObject);

		// Request a viewport update whenever a new item in the pipeline editor is selected, 
		// because the currently selected modifier may be rendering gizmos in the viewports. 
		if(_datasetContainer.currentSet())
			_datasetContainer.currentSet()->viewportConfig()->updateViewports();

		updateActions(currentItem);
	}

	// Whenever no object is selected, show the About Panel containing information about the program.
	if(currentItem == nullptr)
		_aboutRollout->show();
}

/******************************************************************************
* Updates the state of the actions that can be invoked on the currently selected item.
******************************************************************************/
void ModifyCommandPage::updateActions(PipelineListItem* currentItem)
{
	QAction* deleteModifierAction = _actionManager->getAction(ACTION_MODIFIER_DELETE);
	QAction* moveModifierUpAction = _actionManager->getAction(ACTION_MODIFIER_MOVE_UP);
	QAction* moveModifierDownAction = _actionManager->getAction(ACTION_MODIFIER_MOVE_DOWN);
	QAction* toggleModifierStateAction = _actionManager->getAction(ACTION_MODIFIER_TOGGLE_STATE);
	QAction* makeElementIndependentAction = _actionManager->getAction(ACTION_PIPELINE_MAKE_INDEPENDENT);

	_modifierSelector->setEnabled(currentItem != nullptr);
	RefTarget* currentObject = currentItem ? currentItem->object() : nullptr;

	if(ModifierApplication* modApp = dynamic_object_cast<ModifierApplication>(currentObject)) {
		deleteModifierAction->setEnabled(true);
		moveModifierDownAction->setEnabled(
			dynamic_object_cast<ModifierApplication>(modApp->input()) != nullptr &&
			modApp->input()->isPipelineBranch(true) == false);
		int index = std::find(_pipelineListModel->items().begin(), _pipelineListModel->items().end(), currentItem) - _pipelineListModel->items().begin();
		moveModifierUpAction->setEnabled(
			index > 0 &&
			dynamic_object_cast<ModifierApplication>(_pipelineListModel->item(index - 1)->object()) != nullptr &&
			modApp->isPipelineBranch(true) == false);
		toggleModifierStateAction->setEnabled(true);
		toggleModifierStateAction->setChecked(modApp->modifier() && modApp->modifier()->isEnabled() == false);
	}
	else {
		deleteModifierAction->setEnabled(false);
		moveModifierUpAction->setEnabled(false);
		moveModifierDownAction->setEnabled(false);
		toggleModifierStateAction->setChecked(false);
		toggleModifierStateAction->setEnabled(false);
	}

	makeElementIndependentAction->setEnabled(PipelineListModel::isSharedObject(currentObject));
}

/******************************************************************************
* Handles the ACTION_MODIFIER_DELETE command.
******************************************************************************/
void ModifyCommandPage::onDeleteModifier()
{
	// Get the currently selected modifier.
	PipelineListItem* selectedItem = pipelineListModel()->selectedItem();
	if(!selectedItem) return;

	OORef<ModifierApplication> modApp = dynamic_object_cast<ModifierApplication>(selectedItem->object());
	if(!modApp) return;

	UndoableTransaction::handleExceptions(_datasetContainer.currentSet()->undoStack(), tr("Delete modifier"), [modApp, this]() {
		auto dependentsList = modApp->dependents();
		for(RefMaker* dependent : dependentsList) {
			if(ModifierApplication* precedingModApp = dynamic_object_cast<ModifierApplication>(dependent)) {
				if(precedingModApp->input() == modApp) {
					precedingModApp->setInput(modApp->input());
					pipelineListModel()->setNextToSelectObject(modApp->input());
				}
			}
			else if(PipelineSceneNode* pipeline = dynamic_object_cast<PipelineSceneNode>(dependent)) {
				if(pipeline->dataProvider() == modApp) {
					pipeline->setDataProvider(modApp->input());
					pipelineListModel()->setNextToSelectObject(pipeline->dataProvider());
				}
			}
		}
		OORef<Modifier> modifier = modApp->modifier();
		modApp->setInput(nullptr);
		modApp->setModifier(nullptr);

		// Delete modifier if there are no more applications left.
		if(modifier->modifierApplications().empty())
			modifier->deleteReferenceObject();
	});
}

/******************************************************************************
* This called when the user double clicks on an item in the modifier stack.
******************************************************************************/
void ModifyCommandPage::onModifierStackDoubleClicked(const QModelIndex& index)
{
	PipelineListItem* item = pipelineListModel()->item(index.row());
	OVITO_CHECK_OBJECT_POINTER(item);

	if(ModifierApplication* modApp = dynamic_object_cast<ModifierApplication>(item->object())) {
		// Toggle enabled state of modifier.
		UndoableTransaction::handleExceptions(_datasetContainer.currentSet()->undoStack(), tr("Toggle modifier state"), [modApp]() {
			modApp->modifier()->setEnabled(!modApp->modifier()->isEnabled());
		});
	}

	if(DataVis* vis = dynamic_object_cast<DataVis>(item->object())) {
		// Toggle enabled state of vis element.
		UndoableTransaction::handleExceptions(_datasetContainer.currentSet()->undoStack(), tr("Toggle visual element"), [vis]() {
			vis->setEnabled(!vis->isEnabled());
		});
	}
}

/******************************************************************************
* Handles the ACTION_MODIFIER_MOVE_UP command, which moves the selected
* modifier up one position in the stack.
******************************************************************************/
void ModifyCommandPage::onModifierMoveUp()
{
	PipelineListItem* selectedItem = pipelineListModel()->selectedItem();
	if(!selectedItem) return;
	OORef<ModifierApplication> modApp = dynamic_object_cast<ModifierApplication>(selectedItem->object());
	if(!modApp) return;

	UndoableTransaction::handleExceptions(_datasetContainer.currentSet()->undoStack(), tr("Move modifier up"), [modApp]() {
		OVITO_ASSERT(modApp->isPipelineBranch(true) == false);
		for(RefMaker* dependent : modApp->dependents()) {
			if(OORef<ModifierApplication> predecessor = dynamic_object_cast<ModifierApplication>(dependent)) {
				if(predecessor->pipelines(true).empty()) continue;
				const auto dependentsList2 = predecessor->dependents();
				for(RefMaker* dependent2 : dependentsList2) {
					if(ModifierApplication* predecessor2 = dynamic_object_cast<ModifierApplication>(dependent2)) {
						predecessor2->setInput(modApp);
					}
					else if(PipelineSceneNode* predecessor2 = dynamic_object_cast<PipelineSceneNode>(dependent2)) {
						predecessor2->setDataProvider(modApp);
					}
				}
				predecessor->setInput(modApp->input());
				modApp->setInput(predecessor);
				break;
			}
		}
	});
}

/******************************************************************************
* Handles the ACTION_MODIFIER_MOVE_DOWN command, which moves the selected
* modifier down one position in the stack.
******************************************************************************/
void ModifyCommandPage::onModifierMoveDown()
{
	PipelineListItem* selectedItem = pipelineListModel()->selectedItem();
	if(!selectedItem) return;
	OORef<ModifierApplication> modApp = dynamic_object_cast<ModifierApplication>(selectedItem->object());
	if(!modApp) return;

	UndoableTransaction::handleExceptions(_datasetContainer.currentSet()->undoStack(), tr("Move modifier down"), [modApp]() {
		if(OORef<ModifierApplication> successor = dynamic_object_cast<ModifierApplication>(modApp->input())) {
			OVITO_ASSERT(successor->isPipelineBranch(true) == false);
			const auto dependentsList = modApp->dependents();
			for(RefMaker* dependent : dependentsList) {
				if(ModifierApplication* predecessor = dynamic_object_cast<ModifierApplication>(dependent)) {
					predecessor->setInput(successor);
				}
				else if(PipelineSceneNode* predecessor = dynamic_object_cast<PipelineSceneNode>(dependent)) {
					predecessor->setDataProvider(successor);
				}
			}
			modApp->setInput(successor->input());
			successor->setInput(modApp);
		}
	});
}

/******************************************************************************
* Handles the ACTION_MODIFIER_TOGGLE_STATE command, which toggles the
* enabled/disable state of the selected modifier.
******************************************************************************/
void ModifyCommandPage::onModifierToggleState(bool newState)
{
	// Get the selected modifier from the modifier stack box.
	QModelIndexList selection = _pipelineWidget->selectionModel()->selectedRows();
	if(selection.empty())
		return;

	onModifierStackDoubleClicked(selection.front());
}

/******************************************************************************
* Handles the ACTION_PIPELINE_MAKE_INDEPENDENT command.
******************************************************************************/
void ModifyCommandPage::onMakeElementIndependent()
{
	// Get the currently selected item.
	PipelineListItem* selectedItem = pipelineListModel()->selectedItem();
	if(!selectedItem) return;

	if(DataVis* visElement = dynamic_object_cast<DataVis>(selectedItem->object())) {
		UndoableTransaction::handleExceptions(_datasetContainer.currentSet()->undoStack(), tr("Make visual element independent"), [this,visElement]() {
			PipelineSceneNode* node = pipelineListModel()->selectedNode();
			DataVis* replacementVisElement = node->makeVisElementIndependent(visElement);
			pipelineListModel()->setNextToSelectObject(replacementVisElement);
		});
	}
	else if(PipelineObject* selectedPipelineObj = dynamic_object_cast<PipelineObject>(selectedItem->object())) {
		UndoableTransaction::handleExceptions(_datasetContainer.currentSet()->undoStack(), tr("Make pipeline element independent"), [this,selectedPipelineObj]() {
			CloneHelper cloneHelper;
			OORef<PipelineObject> currentObj = pipelineListModel()->selectedNode()->dataProvider();
			ModifierApplication* predecessorModApp = nullptr;
			// Go up the pipeline, starting at the node, until we reach the selected pipeline object.
			// Duplicate all shared ModifierApplications to remove pipeline branches.
			// When ariving at the selected modifier application, duplicate the modifier too
			// in case it is shared by multiple pipelines.
			while(currentObj) {
				PipelineObject* nextObj = nullptr;
				if(ModifierApplication* modApp = dynamic_object_cast<ModifierApplication>(currentObj)) {
					if(modApp->pipelines(true).size() > 1) {
						OORef<ModifierApplication> clonedModApp = cloneHelper.cloneObject(modApp, false);
						if(predecessorModApp)
							predecessorModApp->setInput(clonedModApp);
						else
							pipelineListModel()->selectedNode()->setDataProvider(clonedModApp);
						predecessorModApp = clonedModApp;
						pipelineListModel()->setNextToSelectObject(clonedModApp);
					}
					else {
						predecessorModApp = modApp;
					}
					if(currentObj == selectedPipelineObj) {
						if(predecessorModApp->modifier()) {
							QSet<PipelineSceneNode*> pipelines;
							for(ModifierApplication* modApp : predecessorModApp->modifier()->modifierApplications())
								pipelines.unite(modApp->pipelines(true));
							if(pipelines.size() > 1)
								predecessorModApp->setModifier(cloneHelper.cloneObject(predecessorModApp->modifier(), true));
						}
						break;
					}
					currentObj = predecessorModApp->input();
				}
				else if(currentObj == selectedPipelineObj) {
					if(currentObj->pipelines(true).size() > 1) {
						OORef<PipelineObject> clonedObject = cloneHelper.cloneObject(currentObj, false);
						if(predecessorModApp)
							predecessorModApp->setInput(clonedObject);
						else
							pipelineListModel()->selectedNode()->setDataProvider(clonedObject);
					}
					break;
				}
				else {
					OVITO_ASSERT(false);
					break;
				}
			}
		});
	}
}

/******************************************************************************
* Creates the rollout panel that shows information about the application
* whenever no object is selected.
******************************************************************************/
void ModifyCommandPage::createAboutPanel()
{
	QWidget* rollout = new QWidget();
	QVBoxLayout* layout = new QVBoxLayout(rollout);
	layout->setContentsMargins(8,8,8,8);

	QTextBrowser* aboutLabel = new QTextBrowser(rollout);
	aboutLabel->setObjectName("AboutLabel");
	aboutLabel->setOpenExternalLinks(true);
	aboutLabel->setMinimumHeight(600);
	aboutLabel->setFrameStyle(QFrame::NoFrame | QFrame::Plain);
	aboutLabel->viewport()->setAutoFillBackground(false);
	layout->addWidget(aboutLabel);

	QByteArray newsPage;
#if !defined(OVITO_BUILD_APPSTORE_VERSION)
	QSettings settings;
	if(settings.value("updates/check_for_updates", true).toBool()) {
		// Retrieve cached news page from settings store.
		newsPage = settings.value("news/cached_webpage").toByteArray();
	}
	if(newsPage.isEmpty()) {
		QResource res("/gui/mainwin/command_panel/about_panel.html");
		newsPage = QByteArray((const char *)res.data(), (int)res.size());
	}
#else
	QResource res("/gui/mainwin/command_panel/about_panel_no_updates.html");
	newsPage = QByteArray((const char *)res.data(), (int)res.size());
#endif

	_aboutRollout = _propertiesPanel->addRollout(rollout, Application::applicationName());
	showProgramNotice(QString::fromUtf8(newsPage.constData()));

#if !defined(OVITO_BUILD_APPSTORE_VERSION)
	if(settings.value("updates/check_for_updates", true).toBool()) {
		constexpr int INSTALLATION_ID_LENGTH = 18;

		// Retrieve/generate unique installation id.
		QByteArray id;
		if(settings.value("updates/transmit_id", true).toBool()) {
			if(settings.contains("installation/id")) {
				id = QByteArray::fromHex(settings.value("installation/id").toString().toLatin1());
				if(id == QByteArray(INSTALLATION_ID_LENGTH, '\0') || id.size() != INSTALLATION_ID_LENGTH)
					id.clear();
			}
			if(id.isEmpty()) {
				// Generate a new unique ID.
				id.fill('0', INSTALLATION_ID_LENGTH);
				std::random_device rdev;
				std::uniform_int_distribution<> rdist(0, 0xFF);
				for(auto& c : id)
					c = (char)rdist(rdev);
				settings.setValue("installation/id", QVariant::fromValue(QString::fromLatin1(id.toHex())));
			}
		}
		else {
			id.fill(0, INSTALLATION_ID_LENGTH);
		}

		QString operatingSystemString;
#if defined(Q_OS_MACOS)
		operatingSystemString = QStringLiteral("macosx");
#elif defined(Q_OS_WIN) || defined(Q_OS_CYGWIN)
		operatingSystemString = QStringLiteral("win");
#elif defined(Q_OS_LINUX) || defined(Q_OS_UNIX)
		operatingSystemString = QStringLiteral("linux");
#else
		operatingSystemString = QStringLiteral("other");
#endif

		QString programEdition;
#if defined(OVITO_BUILD_BASIC)
		programEdition = QStringLiteral("basic/");
#elif defined(OVITO_BUILD_PRO)
		programEdition = QStringLiteral("pro/");
#endif

		// Fetch newest web page from web server.
		QString urlString = QString("https://www.ovito.org/appnews/v%1.%2.%3/%4?ovito=%5&OS=%6%7")
				.arg(Application::applicationVersionMajor())
				.arg(Application::applicationVersionMinor())
				.arg(Application::applicationVersionRevision())
				.arg(programEdition)
				.arg(QString(id.toHex()))
				.arg(operatingSystemString)
				.arg(QT_POINTER_SIZE*8);
		QNetworkAccessManager* networkAccessManager = Application::instance()->networkAccessManager();
		QNetworkReply* networkReply = networkAccessManager->get(QNetworkRequest(QUrl(urlString)));
		connect(networkReply, &QNetworkReply::finished, this, &ModifyCommandPage::onWebRequestFinished);
	}
#endif
}

/******************************************************************************
* Is called by the system when fetching the news web page from the server is
* completed.
******************************************************************************/
void ModifyCommandPage::onWebRequestFinished()
{
	QNetworkReply* reply = qobject_cast<QNetworkReply*>(QObject::sender());
	if(reply->error() == QNetworkReply::NoError) {
		QByteArray page = reply->readAll();
		reply->close();
		if(page.startsWith("<html><!--OVITO-->")) {

			showProgramNotice(QString::fromUtf8(page.constData()));

			QSettings settings;
			settings.setValue("news/cached_webpage", page);
		}
	}
	reply->deleteLater();
}

/******************************************************************************
* Displays the given HTML page content in the About pane.
******************************************************************************/
void ModifyCommandPage::showProgramNotice(const QString& htmlPage)
{
	QString finalText = htmlPage;

#ifdef OVITO_EXPIRATION_DATE
	QDate expirationDate = QDate::fromString(QStringLiteral(OVITO_EXPIRATION_DATE), Qt::ISODate);
	QDate currentDate = QDate::currentDate();
	QString expirationNotice;
	if(currentDate < expirationDate.addMonths(-1)) {
		expirationNotice = tr("<h4>Preview build</h4><p>You are using a preview version of %1, which will expire on %2. "
				"The final release of %1 will be made available on our website <a href=\"https://www.ovito.org/\">www.ovito.org</a>.</p>")
			.arg(Application::applicationName())
			.arg(expirationDate.toString(Qt::SystemLocaleShortDate));
	}
	else if(currentDate <= expirationDate) {
		expirationNotice = tr("<h4>Preview build: Expiration date approaching</h4><p style=\"background-color: rgb(230,180,180);\">You are using a preview version of %1, which will expire on %2. "
				"The final program release is now available. Please visit our website <a href=\"https://www.ovito.org/\">www.ovito.org</a>. "
				"<br>This preview release will stop working in %3 days!</p>")
			.arg(Application::applicationName())
			.arg(expirationDate.toString(Qt::SystemLocaleShortDate))
			.arg(currentDate.daysTo(expirationDate));
	}
	else {
		expirationNotice = tr("<h4>Preview build</h4><p style=\"background-color: rgb(230,180,180);\">This preview version of %1 has expired on %2 and will no longer work. "
				"The final program release is now available, please visit our website <a href=\"https://www.ovito.org/\">www.ovito.org</a>. ")
			.arg(Application::applicationName())
			.arg(expirationDate.toString(Qt::SystemLocaleShortDate));
	}
	finalText.replace(QStringLiteral("<p>&nbsp;</p>"), expirationNotice);
#endif

	QTextBrowser* aboutLabel = _aboutRollout->findChild<QTextBrowser*>("AboutLabel");
	OVITO_CHECK_POINTER(aboutLabel);
	aboutLabel->setHtml(finalText);
}


}	// End of namespace
