////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2018 Alexander Stukowski
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

#include <ovito/stdobj/gui/StdObjGui.h>
#include <ovito/stdobj/gui/io/DataSeriesPlotExporter.h>
#include <ovito/stdobj/io/DataSeriesExporter.h>
#include <ovito/gui/mainwin/MainWindow.h>
#include <ovito/gui/dialogs/FileExporterSettingsDialog.h>
#include <ovito/gui/dialogs/HistoryFileDialog.h>
#include <ovito/gui/utilities/concurrent/ProgressDialog.h>
#include <ovito/core/utilities/concurrent/AsyncOperation.h>
#include "SeriesInspectionApplet.h"

namespace Ovito { namespace StdObj {

IMPLEMENT_OVITO_CLASS(SeriesInspectionApplet);

/******************************************************************************
* Lets the applet create the UI widget that is to be placed into the data
* inspector panel.
******************************************************************************/
QWidget* SeriesInspectionApplet::createWidget(MainWindow* mainWindow)
{
	createBaseWidgets();
	_mainWindow = mainWindow;

	QSplitter* splitter = new QSplitter();
	splitter->addWidget(containerSelectionWidget());

	QWidget* rightContainer = new QWidget();
	splitter->addWidget(rightContainer);
	splitter->setStretchFactor(0, 1);
	splitter->setStretchFactor(1, 4);

	QHBoxLayout* rightLayout = new QHBoxLayout(rightContainer);
	rightLayout->setContentsMargins(0,0,0,0);
	rightLayout->setSpacing(0);

	QToolBar* toolbar = new QToolBar();
	toolbar->setOrientation(Qt::Vertical);
	toolbar->setToolButtonStyle(Qt::ToolButtonIconOnly);
	toolbar->setIconSize(QSize(22,22));
	toolbar->setStyleSheet("QToolBar { padding: 0px; margin: 0px; border: 0px none black; spacing: 0px; }");

	QActionGroup* plotTypeActionGroup = new QActionGroup(this);
	_switchToPlotAction = plotTypeActionGroup->addAction(QIcon(":/gui/mainwin/inspector/show_chart.svg"), tr("Chart view"));
	_switchToTableAction = plotTypeActionGroup->addAction(QIcon(":/gui/mainwin/inspector/table_chart.svg"), tr("Data table view"));
	toolbar->addAction(_switchToPlotAction);
	toolbar->addAction(_switchToTableAction);
	_switchToPlotAction->setCheckable(true);
	_switchToTableAction->setCheckable(true);
	_switchToPlotAction->setChecked(true);
	toolbar->addSeparator();

	_exportSeriesToFileAction = new QAction(QIcon(":/gui/actions/file/file_save_as.bw.svg"), tr("Export data plot"), this);
	connect(_exportSeriesToFileAction, &QAction::triggered, this, &SeriesInspectionApplet::exportDataToFile);
	toolbar->addAction(_exportSeriesToFileAction);

	_stackedWidget = new QStackedWidget();
	rightLayout->addWidget(_stackedWidget, 1);
	rightLayout->addWidget(toolbar, 0);

	connect(_switchToPlotAction, &QAction::triggered, this, [this]() {
		_stackedWidget->setCurrentIndex(0);
		_exportSeriesToFileAction->setToolTip(tr("Export data plot"));
	});
	connect(_switchToTableAction, &QAction::triggered, this, [this]() {
		_stackedWidget->setCurrentIndex(1);
		_exportSeriesToFileAction->setToolTip(tr("Export data to text file"));
	});

	_plotWidget = new DataSeriesPlotWidget();
	_stackedWidget->addWidget(_plotWidget);
	_stackedWidget->addWidget(tableView());

	return splitter;
}

/******************************************************************************
* Is called when the user selects a different container object from the list.
******************************************************************************/
void SeriesInspectionApplet::currentContainerChanged()
{
	PropertyInspectionApplet::currentContainerChanged();

	// Update the displayed plot.
	plotWidget()->setSeries(static_object_cast<DataSeriesObject>(selectedContainerObject()));

	// Update actions.
	_exportSeriesToFileAction->setEnabled(plotWidget()->series() != nullptr);
}

/******************************************************************************
* Selects a specific data object in this applet.
******************************************************************************/
bool SeriesInspectionApplet::selectDataObject(PipelineObject* dataSource, const QString& objectIdentifierHint, const QVariant& modeHint)
{
	// Let the base class switch to the right data series object. 
	bool result = PropertyInspectionApplet::selectDataObject(dataSource, objectIdentifierHint, modeHint);
	
	if(result) {
		// The mode hint is used to switch between plot and table view.
		int mode = modeHint.toInt();
		if(mode == 0) {
			_switchToPlotAction->trigger(); // Plot view
		}
		else {
			_switchToTableAction->trigger(); // Table view
		}
	}

	return result;
}

/******************************************************************************
* Exports the current data series to a text file.
******************************************************************************/
void SeriesInspectionApplet::exportDataToFile()
{
	const DataSeriesObject* series = plotWidget()->series();
	if(!series)
		return;

	// Let the user select a destination file.
	HistoryFileDialog dialog("export", _mainWindow, tr("Export Data Series"));
	QString filterString;
	if(_stackedWidget->currentIndex() == 0)
		filterString = QStringLiteral("%1 (%2)").arg(DataSeriesPlotExporter::OOClass().fileFilterDescription(), DataSeriesPlotExporter::OOClass().fileFilter());
	else
		filterString = QStringLiteral("%1 (%2)").arg(DataSeriesExporter::OOClass().fileFilterDescription(), DataSeriesExporter::OOClass().fileFilter());
	dialog.setNameFilter(filterString);
	dialog.setOption(QFileDialog::DontUseNativeDialog);
	dialog.setAcceptMode(QFileDialog::AcceptSave);
	dialog.setFileMode(QFileDialog::AnyFile);
	dialog.setConfirmOverwrite(true);

	// Go to the last directory used.
	QSettings settings;
	settings.beginGroup("file/export");
	QString lastExportDirectory = settings.value("last_export_dir").toString();
	if(!lastExportDirectory.isEmpty())
		dialog.setDirectory(lastExportDirectory);

	if(!dialog.exec() || dialog.selectedFiles().empty())
		return;
	QString exportFile = dialog.selectedFiles().front();

	// Remember directory for the next time...
	settings.setValue("last_export_dir", dialog.directory().absolutePath());

	// Export to selected file.
	try {
		// Create exporter service.
		OORef<FileExporter> exporter;
		if(_stackedWidget->currentIndex() == 0)
			exporter = new DataSeriesPlotExporter(series->dataset());
		else
			exporter = new DataSeriesExporter(series->dataset());

		// Load user-defined default settings.
		exporter->loadUserDefaults();

		// Pass output filename to exporter.
		exporter->setOutputFilename(exportFile);

		// Set scene node to be exported.
		exporter->setNodeToExport(currentSceneNode());

		// Set data series to be exported.
		exporter->setDataObjectToExport(DataObjectReference(&DataSeriesObject::OOClass(), series->identifier(), series->title()));

		// Let the user adjust the export settings.
		FileExporterSettingsDialog settingsDialog(_mainWindow, exporter);
		if(settingsDialog.exec() != QDialog::Accepted)
			return;

		// Show progress dialog.
		ProgressDialog progressDialog(_mainWindow, tr("File export"));

		// Let the exporter do its job.
		exporter->doExport(progressDialog.taskManager().createMainThreadOperation<>(true));
	}
	catch(const Exception& ex) {
		ex.reportError();
	}
}

}	// End of namespace
}	// End of namespace
