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


#include <ovito/particles/Particles.h>
#include <ovito/core/dataset/pipeline/Modifier.h>
#include <ovito/core/dataset/animation/AnimationSettings.h>

#include <QTcpSocket>

namespace Ovito { namespace Particles {

/**
 * \brief A modifier that updates the particle positions using real-time MD trajectory data it receives from a 
 *        remote network server process using the Interactive Molecular Dynamics (IMD) protocol.
 */
class OVITO_PARTICLES_EXPORT InteractiveMolecularDynamicsModifier : public Modifier
{
	/// Give this modifier class its own metaclass.
	class OOMetaClass : public ModifierClass
	{
	public:

		/// Inherit constructor from base class.
		using ModifierClass::ModifierClass;

		/// Asks the metaclass whether the modifier can be applied to the given input data.
		virtual bool isApplicableTo(const DataCollection& input) const override;
	};

	Q_OBJECT
	OVITO_CLASS_META(InteractiveMolecularDynamicsModifier, OOMetaClass)
	Q_CLASSINFO("DisplayName", "Interactive molecular dynamics");
	Q_CLASSINFO("Description", "Visualize live atomic trajectories from a running MD simulation as they are being calculated.");
	Q_CLASSINFO("ModifierCategory", "Visualization");

public:

	/// Constructor.
	Q_INVOKABLE InteractiveMolecularDynamicsModifier(DataSet* dataset);

	/// Destructor.
	virtual ~InteractiveMolecularDynamicsModifier();

	/// Returns the network socket for the server connection.
	const QTcpSocket& socket() const { return _socket; }

	/// Starts establishing a connection to the remote IMD server.
	void connectToServer();

	/// Disconnects the connection to the remote IMD server.
	void disconnectFromServer();

	/// Asks this object to delete itself.
	virtual void deleteReferenceObject() override;

	/// Modifies the input data synchronously.
	virtual void evaluateSynchronous(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state) override;

protected:

	/// Is called when the value of a non-animatable property field of this RefMaker has changed.
	virtual void propertyChanged(const PropertyFieldDescriptor& field) override;

private Q_SLOTS:

	/// Is called when the connection state of the TCP socket changes.
	void connectionStateChanged();

	/// Is called when a TCP connection error occurred.
	void connectionError(QAbstractSocket::SocketError socketError);

	/// Is called when the TCP socket has received a chunk of data.
	void dataReceived();

private:

	/// Is called to indicate an error state. Closes the connection to the server.
	void protocolError(const QString& errorString, PipelineStatus::StatusType statusType = PipelineStatus::Error);

	/// The header data structure used by the IMD network protocol. 
	struct IMDHeader {
		qint32 type;
		qint32 length;
	};

	/// IMD command message type enumerations.
	enum IMDType {
		IMD_DISCONNECT,   ///< close IMD connection, leaving sim running
		IMD_ENERGIES,     ///< energy data block
		IMD_FCOORDS,      ///< atom coordinates
		IMD_GO,           ///< start the simulation
		IMD_HANDSHAKE,    ///< endianism and version check message
		IMD_KILL,         ///< kill the simulation job, shutdown IMD
		IMD_MDCOMM,       ///< MDComm style force data
		IMD_PAUSE,        ///< pause the running simulation
		IMD_TRATE,        ///< set IMD update transmission rate
		IMD_IOERROR       ///< indicate an I/O error
	};

	/// IMD simulation energy report structure.
	struct IMDEnergies {
		qint32 tstep;     ///< integer timestep index
		float T;          ///< Temperature in degrees Kelvin
		float Etot;       ///< Total energy, in Kcal/mol
		float Epot;       ///< Potential energy, in Kcal/mol
		float Evdw;       ///< Van der Waals energy, in Kcal/mol
		float Eelec;      ///< Electrostatic energy, in Kcal/mol
		float Ebond;      ///< Bond energy, Kcal/mol
		float Eangle;     ///< Angle energy, Kcal/mol
		float Edihe;      ///< Dihedral energy, Kcal/mol
		float Eimpr;      ///< Improper energy, Kcal/mol
	};      


	/// The network hostname of the IMD server.
	DECLARE_MODIFIABLE_PROPERTY_FIELD_FLAGS(QString, hostName, setHostName, PROPERTY_FIELD_MEMORIZE);

	/// The network port of the IMD server.
	DECLARE_MODIFIABLE_PROPERTY_FIELD_FLAGS(int, port, setPort, PROPERTY_FIELD_MEMORIZE);

	/// Controls how often MD frames are transmitted from the server to the client.
	DECLARE_MODIFIABLE_PROPERTY_FIELD_FLAGS(int, transmissionInterval, setTransmissionInterval, PROPERTY_FIELD_MEMORIZE);

	/// The network socket for the server connection.
	QTcpSocket _socket;

	/// The current status of the modifier.
	PipelineStatus _modifierStatus;

	/// The detected endianess of the server side. 
	QSysInfo::Endian _serverEndianess;

	// The last IMD protocol header received.
	IMDHeader _header;

	/// Indicates that the connection to the IMD server has been fully established.
	bool _isConnected = false;

	/// The number of bytes to receive as part of the current message.
	qint64 _messageBytesToReceive;

	/// The last particle coordinates received from the server.
	PropertyPtr _coordinates;

	/// The number of frames received from the server so far.
	qint64 _numFramesReceived;

	/// Indicates that a pipeline update is pending, because the previous update is still in progress.
	bool _pipelineUpdatePending = false;

	/// This utility object is used to suppress preliminary viewport updates when playing an IMD trajectory.
	boost::optional<PreliminaryViewportUpdatesSuspender> _preliminaryUpdateSuspender;
};

}	// End of namespace
}	// End of namespace
