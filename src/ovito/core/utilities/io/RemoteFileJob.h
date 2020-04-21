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


#include <ovito/core/Core.h>
#include <ovito/core/utilities/concurrent/Future.h>

#include <QQueue>

namespace Ovito {

#ifdef OVITO_SSH_CLIENT
namespace Ssh {
	// These classes are defined elsewhere:
	class SshConnection;
	class ScpChannel;
	class LsChannel;
	class CatChannel;
}
#endif

/**
 * \brief Base class for background jobs that access remote files and directories via SSH.
 */
class RemoteFileJob : public QObject
{
	Q_OBJECT

public:

	/// Constructor.
	RemoteFileJob(QUrl url, PromiseBase& promise);

	/// Returns the URL being accessed.
	const QUrl& url() const { return _url; }

protected:

	/// Opens the network connection.
	Q_INVOKABLE void start();

	/// Closes the network connection.
	virtual void shutdown(bool success);

protected Q_SLOTS:

#ifdef OVITO_SSH_CLIENT
	/// Handles network connection errors.
	void connectionError();

	/// Handles network authentication errors.
	void authenticationFailed();

	/// Is called when the network connection has been established.
    virtual void connectionEstablished() = 0;
#endif

	/// Handles network connection cancelation by user.
	void connectionCanceled();

	/// Handles QNetworkReply finished signals.
	void networkReplyFinished();

	/// Handles QNetworkReply progress signals.
	virtual void networkReplyDownloadProgress(qint64 bytesReceived, qint64 bytesTotal) {}

protected:

	/// The URL of the file or directory.
	const QUrl _url;

#ifdef OVITO_SSH_CLIENT
	/// The SSH connection.
	Ovito::Ssh::SshConnection* _connection = nullptr;
#endif

	/// The Qt network request reply.
	QNetworkReply* _networkReply = nullptr;

    /// The associated asynchronous task of the job.
    PromiseBase& _promise;

    /// Indicates whether this job is currently active.
    bool _isActive = false;

    /// Queue of jobs that are waiting to be executed.
    static QQueue<RemoteFileJob*> _queuedJobs;

    /// Keeps track of how many jobs are currently active.
    static int _numActiveJobs;
};

/**
 * \brief A background jobs that downloads a file stored on a remote host to the local computer.
 */
class DownloadRemoteFileJob : public RemoteFileJob
{
	Q_OBJECT

public:

	/// Constructor.
	DownloadRemoteFileJob(QUrl url, TaskManager& taskManager) :
		RemoteFileJob(std::move(url), _promise), 
		_promise(Promise<FileHandle>::createAsynchronousOperation(taskManager, false)) {}

	/// Returns a future yielding the file downloaded by this job.
	SharedFuture<FileHandle> sharedFuture() {
		return _promise.sharedFuture();
	}

protected:

	/// Closes the network connection.
	virtual void shutdown(bool success) override;

#ifdef OVITO_SSH_CLIENT
	/// Is called when the network connection has been established.
    virtual void connectionEstablished() override;
#endif

	/// Handles QNetworkReply progress signals.
	virtual void networkReplyDownloadProgress(qint64 bytesReceived, qint64 bytesTotal) override;

	/// Writes the data received from the server so far to the local file. 
	void storeReceivedData();

protected Q_SLOTS:

#ifdef OVITO_SSH_CLIENT
    /// Is called when the remote host starts sending the file.
	void receivingFile(qint64 fileSize);

    /// Is called when the remote host sent some file data.
	void receivedData(qint64 totalReceivedBytes);

	/// Is called after the file has been downloaded.
	void receivedFileComplete();

    /// Is called when an SCP error occurs in the channel.
    void channelError();

	/// Handles SSH channel close.
	void channelClosed();
#endif

private:

#ifdef OVITO_SSH_CLIENT
	/// The SCP channel.
    Ovito::Ssh::ScpChannel* _scpChannel = nullptr;
#endif

    /// The local copy of the file.
    QScopedPointer<QTemporaryFile> _localFile;

	/// The memory-mapped destination file.
	uchar* _fileMapping = nullptr;

	/// The promise through which the result of this download job is returned.
	Promise<FileHandle> _promise;
};

/**
 * \brief A background jobs that lists the files in a directory on a remote host.
 */
class ListRemoteDirectoryJob : public RemoteFileJob
{
	Q_OBJECT

public:

	/// Constructor.
	ListRemoteDirectoryJob(QUrl url, TaskManager& taskManager) :
		RemoteFileJob(std::move(url), _promise), 
		_promise(Promise<QStringList>::createAsynchronousOperation(taskManager, false)) {}

	/// Returns a future yielding the file list downloaded by this job.
	Future<QStringList> future() {
		return _promise.future();
	}

protected:

	/// Closes the network connection.
	virtual void shutdown(bool success) override;

#ifdef OVITO_SSH_CLIENT
	/// Is called when the network connection has been established.
    virtual void connectionEstablished() override;
#endif

protected Q_SLOTS:

#ifdef OVITO_SSH_CLIENT
    /// Is called before transmission of the directory listing begins.
    void receivingDirectory();

    /// Is called after the directory listing has been fully transmitted.
    void receivedDirectoryComplete(const QStringList& listing);

    /// Is called when an error occurs in the SSH channel.
    void channelError();

	/// Handles SSH channel close.
	void channelClosed();
#endif

private:

#ifdef OVITO_SSH_CLIENT
	/// The listing channel.
    Ovito::Ssh::LsChannel* _lsChannel = nullptr;
#endif

	/// The promise through which the result of this job is returned.
	Promise<QStringList> _promise;
};

}	// End of namespace
