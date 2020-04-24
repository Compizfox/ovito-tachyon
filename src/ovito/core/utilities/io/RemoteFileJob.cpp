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

#include <ovito/core/Core.h>
#include <ovito/core/utilities/concurrent/Future.h>
#include <ovito/core/utilities/concurrent/TaskWatcher.h>
#include <ovito/core/utilities/io/FileManager.h>
#include <ovito/core/app/Application.h>
#ifdef OVITO_SSH_CLIENT
	#include <ovito/core/utilities/io/ssh/SshConnection.h>
	#include <ovito/core/utilities/io/ssh/ScpChannel.h>
	#include <ovito/core/utilities/io/ssh/LsChannel.h>
	#include <ovito/core/utilities/io/ssh/CatChannel.h>
#endif
#include "RemoteFileJob.h"

#include <QNetworkReply>

namespace Ovito {

using namespace Ovito::Ssh;

/// List SFTP jobs that are waiting to be executed.
QQueue<RemoteFileJob*> RemoteFileJob::_queuedJobs;

/// Tracks of how many jobs are currently active.
int RemoteFileJob::_numActiveJobs = 0;

/// The maximum number of simultaneous jobs at a time.
constexpr int MaximumNumberOfSimultaneousJobs = 2;

/******************************************************************************
* Constructor.
******************************************************************************/
RemoteFileJob::RemoteFileJob(QUrl url, PromiseBase& promise) :
		_url(std::move(url)), _promise(promise)
{
	// Run all event handlers of this class in the main thread.
	moveToThread(QCoreApplication::instance()->thread());

	// Start download process in the main thread.
	QMetaObject::invokeMethod(this, "start", Qt::QueuedConnection);
}

/******************************************************************************
* Opens the network connection.
******************************************************************************/
void RemoteFileJob::start()
{
	if(!_isActive) {
		// Keep a counter of active jobs.
		// If there are too many jobs active simultaneously, queue them to be executed later.
		if(_numActiveJobs >= MaximumNumberOfSimultaneousJobs) {
			_queuedJobs.enqueue(this);
			return;
		}
		else {
			_numActiveJobs++;
			_isActive = true;
		}
	}

	// This background task started to run.
	_promise.setStarted();

	// Check if process has already been canceled.
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}

	// Get notified if user has canceled the task.
	_promise.finally([this](const TaskPtr& task) {
		if(task->isCanceled())
			QMetaObject::invokeMethod(this, "connectionCanceled");
	});

	// Show task progress in the GUI.
	_promise.task()->taskManager()->registerPromise(_promise);

#ifdef OVITO_SSH_CLIENT
	if(_url.scheme() == QStringLiteral("sftp")) {
		// Handle sftp URLs.

		SshConnectionParameters connectionParams;
		connectionParams.host = _url.host();
		connectionParams.userName = _url.userName();
		connectionParams.password = _url.password();
		connectionParams.port = _url.port(0);
		_promise.setProgressText(tr("Connecting to remote host %1").arg(connectionParams.host));

		// Open connection
		_connection = Application::instance()->fileManager()->acquireSshConnection(connectionParams);
		OVITO_CHECK_POINTER(_connection);

		// Listen for signals of the connection.
		connect(_connection, &SshConnection::error, this, &RemoteFileJob::connectionError);
		connect(_connection, &SshConnection::canceled, this, &RemoteFileJob::connectionCanceled);
		connect(_connection, &SshConnection::allAuthsFailed, this, &RemoteFileJob::authenticationFailed);
		if(_connection->isConnected()) {
			QTimer::singleShot(0, this, &RemoteFileJob::connectionEstablished);
			return;
		}
		connect(_connection, &SshConnection::connected, this, &RemoteFileJob::connectionEstablished);

		// Start to connect.
		_connection->connectToHost();
	}
	else {
#endif
		// Handle http(s) URLs.

		_promise.setProgressText(tr("Downloading file %1 from %2").arg(_url.fileName()).arg(_url.host()));
		QNetworkAccessManager* networkAccessManager = Application::instance()->networkAccessManager();
		_networkReply = networkAccessManager->get(QNetworkRequest(_url));

		connect(_networkReply, &QNetworkReply::downloadProgress, this, &RemoteFileJob::networkReplyDownloadProgress);
		connect(_networkReply, &QNetworkReply::finished, this, &RemoteFileJob::networkReplyFinished);

#ifdef OVITO_SSH_CLIENT
	}
#endif
}

/******************************************************************************
* Closes the network connection.
******************************************************************************/
void RemoteFileJob::shutdown(bool success)
{
#ifdef OVITO_SSH_CLIENT
	if(_connection) {
		disconnect(_connection, nullptr, this, nullptr);
		Application::instance()->fileManager()->releaseSshConnection(_connection);
		_connection = nullptr;
	}
#endif
	if(_networkReply) {
		disconnect(_networkReply, nullptr, this, nullptr);
		_networkReply->abort();
		_networkReply->deleteLater();
		_networkReply = nullptr;
	}

	_promise.setFinished();

	// Update the counter of active jobs.
	if(_isActive) {
		_numActiveJobs--;
		_isActive = false;
	}

	// Schedule this object for deletion.
	deleteLater();

	// If there jobs waiting in the queue, execute next job.
	if(!_queuedJobs.isEmpty() && _numActiveJobs < MaximumNumberOfSimultaneousJobs) {
		RemoteFileJob* waitingJob = _queuedJobs.dequeue();
		if(!waitingJob->_promise.isCanceled()) {
			waitingJob->start();
		}
		else {
			// Skip canceled jobs.
			waitingJob->_promise.setStarted();
			waitingJob->shutdown(false);
		}
	}
}

#ifdef OVITO_SSH_CLIENT
/******************************************************************************
* Handles SSH connection errors.
******************************************************************************/
void RemoteFileJob::connectionError()
{
	_promise.setException(std::make_exception_ptr(
		Exception(tr("Cannot access URL\n\n%1\n\nSSH connection error: %2").arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded)).
			arg(_connection->errorMessage()))));

	shutdown(false);
}

/******************************************************************************
* Handles SSH authentication errors.
******************************************************************************/
void RemoteFileJob::authenticationFailed()
{
	_promise.setException(std::make_exception_ptr(
		Exception(tr("Cannot access URL\n\n%1\n\nSSH authentication failed").arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded)))));

	shutdown(false);
}
#endif

/******************************************************************************
* Handles cancelation by the user.
******************************************************************************/
void RemoteFileJob::connectionCanceled()
{
	// If user has canceled the connection, cancel the file retrieval operation as well.
	_promise.cancel();
	shutdown(false);
}

/******************************************************************************
* Handles QNetworkReply finished signals.
******************************************************************************/
void RemoteFileJob::networkReplyFinished()
{
	if(_networkReply->error() == QNetworkReply::NoError) {

		shutdown(true);
	}
	else {
		_promise.setException(std::make_exception_ptr(
			Exception(tr("Cannot access URL\n\n%1\n\n%2").arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded)).
				arg(_networkReply->errorString()))));

		shutdown(false);
	}
}

#ifdef OVITO_SSH_CLIENT
/******************************************************************************
* Handles closed SSH channel.
******************************************************************************/
void DownloadRemoteFileJob::channelClosed()
{
	if(!_promise.isFinished()) {
		_promise.setException(std::make_exception_ptr(
			Exception(tr("Cannot access URL\n\n%1\n\nSSH channel closed: %2").arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded)).
				arg(_scpChannel->errorMessage()))));
	}

	shutdown(false);
}

/******************************************************************************
* Is called when the SSH connection has been established.
******************************************************************************/
void DownloadRemoteFileJob::connectionEstablished()
{
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}

	// Open the SCP channel.
	_promise.setProgressText(tr("Opening SCP channel to remote host %1").arg(_connection->hostname()));
	_scpChannel = new ScpChannel(_connection, _url.path());
	connect(_scpChannel, &ScpChannel::receivingFile, this, &DownloadRemoteFileJob::receivingFile);
	connect(_scpChannel, &ScpChannel::receivedData, this, &DownloadRemoteFileJob::receivedData);
	connect(_scpChannel, &ScpChannel::receivedFileComplete, this, &DownloadRemoteFileJob::receivedFileComplete);
	connect(_scpChannel, &ScpChannel::error, this, &DownloadRemoteFileJob::channelError);
	connect(_scpChannel, &ScpChannel::closed, this, &DownloadRemoteFileJob::channelClosed);
	_scpChannel->openChannel();
}

/******************************************************************************
* Is called when the SCP channel failed.
******************************************************************************/
void DownloadRemoteFileJob::channelError()
{
	_promise.setException(std::make_exception_ptr(
		Exception(tr("Cannot access remote URL\n\n%1\n\n%2")
			.arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded))
			.arg(_scpChannel->errorMessage()))));

	shutdown(false);
}
#endif

/******************************************************************************
* Closes the network connection.
******************************************************************************/
void DownloadRemoteFileJob::shutdown(bool success)
{
#ifdef OVITO_SSH_CLIENT
	// Close file channel.
	if(_scpChannel) {
		disconnect(_scpChannel, nullptr, this, nullptr);
		_scpChannel->closeChannel();
		_scpChannel->deleteLater();
		_scpChannel = nullptr;
	}
#endif

	// Write all received data to the local file.
	if(success)
		storeReceivedData();

	// Close local file and clean up.
	if(_localFile) {
		// Make sure the received data was successfully written to the temporary file.
		if(_fileMapping) {
			if(!_localFile->unmap(_fileMapping)) {
				_promise.setException(std::make_exception_ptr(Exception(
					tr("Failed to write to local file %1: %2").arg(_localFile->fileName()).arg(_localFile->errorString()))));
				success = false;
			}
			_fileMapping = nullptr;
		}
		if(!_localFile->flush() || _localFile->error() != QFileDevice::NoError) {
			_promise.setException(std::make_exception_ptr(Exception(
				tr("Failed to write to local file %1: %2").arg(_localFile->fileName()).arg(_localFile->errorString()))));
			success = false;
		}
		_localFile->close();
	}
	if(_localFile && success)
		_promise.setResults(FileHandle(url(), _localFile->fileName()));
	else
		_localFile.reset();

	// Close network connection.
	RemoteFileJob::shutdown(success);

	// Hand downloaded file over to FileManager cache.
	Application::instance()->fileManager()->fileFetched(url(), _localFile.take());
}

#ifdef OVITO_SSH_CLIENT
/******************************************************************************
* Is called when the remote host starts sending the file.
******************************************************************************/
void DownloadRemoteFileJob::receivingFile(qint64 fileSize)
{
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}
	_promise.setProgressMaximum(fileSize);
	_promise.setProgressText(tr("Fetching remote file %1").arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded)));

	// Create the destination file.
	try {
		_localFile.reset(new QTemporaryFile());
		if(!_localFile->open() || !_localFile->resize(fileSize))
			throw Exception(tr("Failed to create temporary file: %1").arg(_localFile->errorString()));

		// Map the file to memory and let the SCP channel write the received data to the memory buffer.
		if(fileSize) {
			_fileMapping = _localFile->map(0, fileSize);
			if(!_fileMapping)
				throw Exception(tr("Failed to map temporary file to memory: %1").arg(_localFile->errorString()));
		}
		_scpChannel->setDestinationBuffer(reinterpret_cast<char*>(_fileMapping));
	}
    catch(Exception&) {
		_promise.captureException();
		shutdown(false);
	}
}

/******************************************************************************
* Is called after the file has been downloaded.
******************************************************************************/
void DownloadRemoteFileJob::receivedFileComplete()
{
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}
    shutdown(true);
}

/******************************************************************************
* Is called when the remote host sent some file data.
******************************************************************************/
void DownloadRemoteFileJob::receivedData(qint64 totalReceivedBytes)
{
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}
	_promise.setProgressValue(totalReceivedBytes);
}
#endif

/******************************************************************************
* Handles QNetworkReply progress signals.
******************************************************************************/
void DownloadRemoteFileJob::networkReplyDownloadProgress(qint64 bytesReceived, qint64 bytesTotal)
{
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}
	if(bytesTotal > 0) {
		_promise.setProgressMaximum(bytesTotal);
		_promise.setProgressValue(bytesReceived);
	}
	storeReceivedData();
}

/******************************************************************************
* Writes the data received from the server so far to the local file. 
******************************************************************************/
void DownloadRemoteFileJob::storeReceivedData()
{
	if(!_networkReply)
		return;

	try {
		// Create the destination file and open it for writing.
		if(!_localFile) {
			_localFile.reset(new QTemporaryFile());
			if(!_localFile->open())
				throw Exception(tr("Failed to create temporary file: %1").arg(_localFile->errorString()));
		}

		// Read data from the network stream.
		QByteArray buffer = _networkReply->read(_networkReply->bytesAvailable());

		// Write data into local file.
		if(_localFile->write(buffer) == -1)
			throw Exception(tr("Failed to write downloaded data to temporary file: %1").arg(_localFile->errorString()));
	}
    catch(Exception&) {
		_promise.captureException();
		shutdown(false);
	}
}

#ifdef OVITO_SSH_CLIENT
/******************************************************************************
* Is called when the SSH connection has been established.
******************************************************************************/
void ListRemoteDirectoryJob::connectionEstablished()
{
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}

	// Open the SCP channel.
	_promise.setProgressText(tr("Opening channel to remote host %1").arg(_connection->hostname()));
	_lsChannel = new LsChannel(_connection, _url.path());
	connect(_lsChannel, &LsChannel::error, this, &ListRemoteDirectoryJob::channelError);
	connect(_lsChannel, &LsChannel::receivingDirectory, this, &ListRemoteDirectoryJob::receivingDirectory);
	connect(_lsChannel, &LsChannel::receivedDirectoryComplete, this, &ListRemoteDirectoryJob::receivedDirectoryComplete);
	connect(_lsChannel, &LsChannel::closed, this, &ListRemoteDirectoryJob::channelClosed);
	_lsChannel->openChannel();
}

/******************************************************************************
* Is called before transmission of the directory listing begins.
******************************************************************************/
void ListRemoteDirectoryJob::receivingDirectory()
{
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}

	// Set progress text.
	_promise.setProgressText(tr("Listing remote directory %1").arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded)));
}

/******************************************************************************
* Is called when the SCP channel failed.
******************************************************************************/
void ListRemoteDirectoryJob::channelError()
{
	_promise.setException(std::make_exception_ptr(
		Exception(tr("Cannot access remote URL\n\n%1\n\n%2")
			.arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded))
			.arg(_lsChannel->errorMessage()))));

	shutdown(false);
}

/******************************************************************************
* Is called after the directory listing has been fully transmitted.
******************************************************************************/
void ListRemoteDirectoryJob::receivedDirectoryComplete(const QStringList& listing)
{
	if(_promise.isCanceled()) {
		shutdown(false);
		return;
	}

	_promise.setResults(listing);
    shutdown(true);
}
#endif

/******************************************************************************
* Closes the SSH connection.
******************************************************************************/
void ListRemoteDirectoryJob::shutdown(bool success)
{
#ifdef OVITO_SSH_CLIENT
	if(_lsChannel) {
		disconnect(_lsChannel, nullptr, this, nullptr);
		_lsChannel->closeChannel();
		_lsChannel->deleteLater();
		_lsChannel = nullptr;
	}
#endif

	RemoteFileJob::shutdown(success);
}

#ifdef OVITO_SSH_CLIENT
/******************************************************************************
* Handles closed SSH channel.
******************************************************************************/
void ListRemoteDirectoryJob::channelClosed()
{
	if(!_promise.isFinished()) {
		_promise.setException(std::make_exception_ptr(
			Exception(tr("Cannot access URL\n\n%1\n\nSSH channel closed: %2").arg(_url.toString(QUrl::RemovePassword | QUrl::PreferLocalFile | QUrl::PrettyDecoded)).
				arg(_lsChannel->errorMessage()))));
	}

	shutdown(false);
}
#endif

}	// End of namespace
