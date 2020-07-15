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
#include <ovito/gui/desktop/properties/BooleanParameterUI.h>
#include <ovito/gui/desktop/properties/IntegerParameterUI.h>
#include <ovito/gui/desktop/dialogs/AnimationSettingsDialog.h>
#include <ovito/core/dataset/animation/AnimationSettings.h>
#include <ovito/core/dataset/io/FileSource.h>
#include "FileSourcePlaybackRateEditor.h"

namespace Ovito {

IMPLEMENT_OVITO_CLASS(FileSourcePlaybackRateEditor);

/******************************************************************************
* Sets up the UI of the editor.
******************************************************************************/
void FileSourcePlaybackRateEditor::createUI(const RolloutInsertionParameters& rolloutParams)
{
	// Create another rollout for animation control.
	QWidget* rollout = createRollout(tr("Animation"), rolloutParams, "data_sources.html#scene_objects.file_source.configure_playback");

	// Create the rollout contents.
	QVBoxLayout* layout = new QVBoxLayout(rollout);
	layout->setContentsMargins(4,4,4,4);
	layout->setSpacing(8);

	_trajectoryModeBtn = new QRadioButton(tr("Animated trajectory:"));
	layout->addWidget(_trajectoryModeBtn);

	QGridLayout* sublayout = new QGridLayout();
	sublayout->setContentsMargins(40,0,0,0);
	sublayout->setSpacing(4);
	sublayout->setColumnStretch(2, 1);
	layout->addLayout(sublayout);

	sublayout->addWidget(new QLabel(tr("Show")), 0, 0, Qt::AlignRight);
	IntegerParameterUI* playbackSpeedNumeratorUI = new IntegerParameterUI(this, PROPERTY_FIELD(FileSource::playbackSpeedNumerator));
	playbackSpeedNumeratorUI->setEnabled(false);
	sublayout->addLayout(playbackSpeedNumeratorUI->createFieldLayout(), 0, 1);
	sublayout->addWidget(new QLabel(tr("trajectory frame(s)")), 0, 2);
	sublayout->addWidget(new QLabel(tr("per")), 1, 0, Qt::AlignRight);
	IntegerParameterUI* playbackSpeedDenominatorUI = new IntegerParameterUI(this, PROPERTY_FIELD(FileSource::playbackSpeedDenominator));
	playbackSpeedDenominatorUI->setEnabled(false);
	sublayout->addLayout(playbackSpeedDenominatorUI->createFieldLayout(), 1, 1);
	sublayout->addWidget(new QLabel(tr("animation frame(s)")), 1, 2);
	sublayout->addWidget(new QLabel(tr("starting at animation frame")), 2, 0, Qt::AlignRight);
	IntegerParameterUI* playbackStartUI = new IntegerParameterUI(this, PROPERTY_FIELD(FileSource::playbackStartTime));
	playbackStartUI->setEnabled(false);
	sublayout->addLayout(playbackStartUI->createFieldLayout(), 2, 1);

	sublayout->setRowMinimumHeight(3, 10);
	sublayout->addWidget(new QLabel(tr("Trajectory length:")), 4, 0, Qt::AlignRight);
	_numTrajectoryFramesDisplay = new QLabel();
	sublayout->addWidget(_numTrajectoryFramesDisplay, 4, 1);
	sublayout->addWidget(new QLabel(tr("Animation length:")), 5, 0, Qt::AlignRight);
	_numAnimationFramesDisplay = new QLabel();
	sublayout->addWidget(_numAnimationFramesDisplay, 5, 1);

	QPushButton* animSettingsBtn = new QPushButton(tr("Animation settings..."));
	sublayout->addWidget(animSettingsBtn, 5, 2);
	connect(animSettingsBtn, &QPushButton::clicked, this, [this]() {
		AnimationSettingsDialog(dataset()->animationSettings(), container()->window()).exec();
		updateInformation();
	});

	connect(_trajectoryModeBtn, &QRadioButton::toggled, playbackSpeedNumeratorUI, &IntegerParameterUI::setEnabled);
	connect(_trajectoryModeBtn, &QRadioButton::toggled, playbackSpeedDenominatorUI, &IntegerParameterUI::setEnabled);
	connect(_trajectoryModeBtn, &QRadioButton::toggled, playbackStartUI, &IntegerParameterUI::setEnabled);
	connect(_trajectoryModeBtn, &QRadioButton::clicked, this, [&](bool checked) {
		if(checked)
			undoableTransaction(tr("Change trajectory playback"), [&]() {
				if(FileSource* fileSource = static_object_cast<FileSource>(editObject()))
					fileSource->setRestrictToFrame(-1);
			});
	});

	_staticFrameModeBtn = new QRadioButton(tr("Static frame:"));
	layout->addSpacing(12);
	layout->addWidget(_staticFrameModeBtn);

	sublayout = new QGridLayout();
	sublayout->setContentsMargins(40,0,0,0);
	sublayout->setSpacing(4);
	sublayout->setColumnStretch(1, 1);
	layout->addLayout(sublayout);

	_staticFrameNumberUI = new IntegerParameterUI(this, PROPERTY_FIELD(FileSource::restrictToFrame));
	_staticFrameNumberUI->setEnabled(false);
	sublayout->addLayout(_staticFrameNumberUI->createFieldLayout(), 0, 0);
	_framesListBox = new QComboBox();
	_framesListBox->setEditable(false);
	_framesListBox->setEnabled(false);
	sublayout->addWidget(_framesListBox, 0, 1);
	// To improve performance of drop-down list display:
	_framesListBox->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
	static_cast<QListView*>(_framesListBox->view())->setUniformItemSizes(true);
	static_cast<QListView*>(_framesListBox->view())->setLayoutMode(QListView::Batched);
	_framesListModel = new QStringListModel(this);
	_framesListBox->setModel(_framesListModel);
	connect(_framesListBox, (void (QComboBox::*)(int))&QComboBox::activated, this, [&](int index) {
		undoableTransaction(tr("Change trajectory playback"), [&]() {
			if(FileSource* fileSource = static_object_cast<FileSource>(editObject()))
				fileSource->setRestrictToFrame(index);
		});
	});

	connect(_staticFrameModeBtn, &QRadioButton::toggled, _framesListBox, &QComboBox::setEnabled);
	connect(_staticFrameModeBtn, &QRadioButton::toggled, _staticFrameNumberUI, &IntegerParameterUI::setEnabled);
	connect(_staticFrameModeBtn, &QRadioButton::clicked, this, [&](bool checked) {
		if(checked)
			undoableTransaction(tr("Change trajectory playback"), [&]() {
				if(FileSource* fileSource = static_object_cast<FileSource>(editObject()))
					fileSource->setRestrictToFrame(fileSource->dataCollectionFrame());
			});
	});

	// Whenever a new FileSource gets loaded into the editor:
	connect(this, &PropertiesEditor::contentsReplaced, this, [this, con1 = QMetaObject::Connection(), con2 = QMetaObject::Connection()](RefTarget* editObject) mutable {
		disconnect(con1);
		disconnect(con2);

		// Update displayed information.
		updateFramesList();
		updateInformation();

		// Update the frames list displayed in the UI whenever it changes.
		con1 = editObject ? connect(static_object_cast<FileSource>(editObject), &FileSource::framesListChanged, this, &FileSourcePlaybackRateEditor::updateFramesList) : QMetaObject::Connection();

		// Update the information display when animation interval changes.
		con2 = editObject ? connect(editObject->dataset()->animationSettings(), &AnimationSettings::intervalChanged, this, &FileSourcePlaybackRateEditor::updateInformation) : QMetaObject::Connection();
	});

	connect(this, &PropertiesEditor::contentsChanged, this, &FileSourcePlaybackRateEditor::updateInformation);
}

/******************************************************************************
* Updates the displayed information in the dialog.
******************************************************************************/
void FileSourcePlaybackRateEditor::updateInformation()
{
	FileSource* fileSource = static_object_cast<FileSource>(editObject());
	if(!fileSource) return; 

	_numTrajectoryFramesDisplay->setText(tr("%n frame(s)", nullptr, fileSource->frames().size()));

	AnimationSettings* anim = fileSource->dataset()->animationSettings();
	_numAnimationFramesDisplay->setText(tr("%n frame(s)", nullptr, anim->lastFrame() - anim->firstFrame() + 1) + 
		(!anim->autoAdjustInterval() ? tr(" (fixed)") : QString()));

	if(fileSource->restrictToFrame() < 0) {
		_staticFrameModeBtn->setChecked(false);
		_trajectoryModeBtn->setChecked(true);
		_framesListBox->setCurrentIndex(fileSource->dataCollectionFrame());
	}
	else {
		_trajectoryModeBtn->setChecked(false);
		_staticFrameModeBtn->setChecked(true);
		_framesListBox->setCurrentIndex(fileSource->restrictToFrame());
	}
}

/******************************************************************************
* Updates the list of trajectory frames displayed in the UI.
******************************************************************************/
void FileSourcePlaybackRateEditor::updateFramesList()
{
	FileSource* fileSource = static_object_cast<FileSource>(editObject());
	if(!fileSource) return;

	QStringList stringList;
	stringList.reserve(fileSource->frames().size());
	for(const FileSourceImporter::Frame& frame : fileSource->frames())
		stringList.push_back(frame.label);
	_framesListModel->setStringList(std::move(stringList));
	if(fileSource->restrictToFrame() >= 0)
		_framesListBox->setCurrentIndex(fileSource->restrictToFrame());
	else
		_framesListBox->setCurrentIndex(fileSource->dataCollectionFrame());
	_staticFrameNumberUI->spinner()->setMinValue(0);
	_staticFrameNumberUI->spinner()->setMaxValue(fileSource->frames().size() - 1);
}

}	// End of namespace
