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

#include <ovito/particles/Particles.h>
#include <ovito/particles/objects/ParticlesObject.h>
#include <ovito/stdobj/properties/PropertyAccess.h>
#include <ovito/stdobj/simcell/SimulationCellObject.h>
#include <ovito/core/dataset/DataSet.h>
#include <ovito/core/utilities/units/UnitsManager.h>
#include <ovito/core/utilities/concurrent/ParallelFor.h>
#include "InteractiveMolecularDynamicsModifier.h"

#include <QtEndian>

namespace Ovito { namespace Particles {

IMPLEMENT_OVITO_CLASS(InteractiveMolecularDynamicsModifier);
DEFINE_PROPERTY_FIELD(InteractiveMolecularDynamicsModifier, hostName);
DEFINE_PROPERTY_FIELD(InteractiveMolecularDynamicsModifier, port);
DEFINE_PROPERTY_FIELD(InteractiveMolecularDynamicsModifier, transmissionInterval);
SET_PROPERTY_FIELD_LABEL(InteractiveMolecularDynamicsModifier, hostName, "Hostname");
SET_PROPERTY_FIELD_LABEL(InteractiveMolecularDynamicsModifier, port, "Port");
SET_PROPERTY_FIELD_LABEL(InteractiveMolecularDynamicsModifier, transmissionInterval, "Transmit every Nth timestep");
SET_PROPERTY_FIELD_UNITS_AND_RANGE(InteractiveMolecularDynamicsModifier, port, IntegerParameterUnit, 0x0000, 0xFFFF);
SET_PROPERTY_FIELD_UNITS_AND_MINIMUM(InteractiveMolecularDynamicsModifier, transmissionInterval, IntegerParameterUnit, 0);

/******************************************************************************
* Constructs the modifier instance.
******************************************************************************/
InteractiveMolecularDynamicsModifier::InteractiveMolecularDynamicsModifier(DataSet* dataset) : Modifier(dataset),
	_hostName(QStringLiteral("localhost")),
	_port(8888),
	_transmissionInterval(1),
	_modifierStatus(PipelineStatus::Warning, tr("IMD connection not established yet."))
{
	// Update modifier status whenever the server connection changes.
	connect(&_socket, &QAbstractSocket::stateChanged, this, &InteractiveMolecularDynamicsModifier::connectionStateChanged);

	// Update modifier status whenever the a connection error occurs.
	connect(&_socket, QOverload<QAbstractSocket::SocketError>::of(&QAbstractSocket::error), this, &InteractiveMolecularDynamicsModifier::connectionError);

	// Process incoming data. 
	connect(&_socket, &QAbstractSocket::readyRead, this, &InteractiveMolecularDynamicsModifier::dataReceived);
}

/******************************************************************************
* Destructor.
******************************************************************************/
InteractiveMolecularDynamicsModifier::~InteractiveMolecularDynamicsModifier()
{
	_socket.disconnect(this);
	_socket.abort();
}

/******************************************************************************
* Asks the modifier whether it can be applied to the given input data.
******************************************************************************/
bool InteractiveMolecularDynamicsModifier::OOMetaClass::isApplicableTo(const DataCollection& input) const
{
	return input.containsObject<ParticlesObject>();
}

/******************************************************************************
* Starts establishing a connection to the remote IMD server.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::connectToServer()
{
	if(_socket.state() == QAbstractSocket::UnconnectedState) {
		_messageBytesToReceive = 0;
		_numFramesReceived = 0;
		_pipelineUpdatePending = false;
		_socket.connectToHost(hostName(), port());
	}
}

/******************************************************************************
* Disconnects the connection to the remote IMD server.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::disconnectFromServer()
{
	_isConnected = false;
	if(_socket.state() == QAbstractSocket::ConnectedState) {
		_socket.disconnectFromHost();
	}
	else {
		_socket.abort();
	}
}

/******************************************************************************
* Is called when the connection state of the TCP socket changes.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::connectionStateChanged()
{
	if(isAboutToBeDeleted()) return;

	switch(_socket.state()) {
	case QAbstractSocket::UnconnectedState:
		if(_modifierStatus.type() != PipelineStatus::Error)
			_modifierStatus = PipelineStatus(PipelineStatus::Warning, tr("Not connected to IMD server"));
		break;
	case QAbstractSocket::HostLookupState:
		_modifierStatus = PipelineStatus(PipelineStatus::Warning, tr("Looking up host name..."));
		break;
	case QAbstractSocket::ConnectingState:
	case QAbstractSocket::BoundState:
		_modifierStatus = PipelineStatus(PipelineStatus::Warning, tr("Connecting to IMD server..."));
		break;
	case QAbstractSocket::ConnectedState:
	case QAbstractSocket::ListeningState:
		_modifierStatus = PipelineStatus(PipelineStatus::Success, tr("Connected to IMD server"));
		break;
	case QAbstractSocket::ClosingState:
		if(_modifierStatus.type() != PipelineStatus::Error)
			_modifierStatus = PipelineStatus(PipelineStatus::Success, tr("Closing connection to IMD server..."));
		_isConnected = false;
		_preliminaryUpdateSuspender.reset();
		break;
	}

	notifyTargetChanged();
}

/******************************************************************************
* Is called when a TCP connection error occurred.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::connectionError(QAbstractSocket::SocketError socketError)
{
	if(socketError != QAbstractSocket::RemoteHostClosedError)
		protocolError(_socket.errorString());
	else
		_modifierStatus = PipelineStatus(PipelineStatus::Warning, tr("IMD connection closed by remote host"));
}

/******************************************************************************
* Is called to indicate an error state. Closes the connection to the server.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::protocolError(const QString& errorString, PipelineStatus::StatusType statusType)
{
	_isConnected = false;
	_socket.abort();
	_modifierStatus = PipelineStatus(statusType, tr("IMD connection error: %1").arg(errorString));
	notifyTargetChanged();
}

/******************************************************************************
* Is called when the TCP socket has received a chunk of data.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::dataReceived()
{
	while(_socket.bytesAvailable() >= sizeof(IMDHeader)) {
		if(_messageBytesToReceive == 0) {

			// Read IMDHeader structure from input data stream.
			if(_socket.read(reinterpret_cast<char*>(&_header), sizeof(_header)) != sizeof(IMDHeader)) {
				protocolError(tr("IMD protocol failure. Did not received sufficient data."));
				break;
			}
			_header.type = qFromBigEndian(_header.type);

			// Expect handshaking code as first message.
			if(_isConnected == false && _header.type != IMD_HANDSHAKE) {
				protocolError(tr("IMD server handshaking failed. Did not receive correct handshaking code from server."));
				break;
			}

			// Handle IMD server messages.
			if(_header.type == IMD_HANDSHAKE) {
				const qint32 IMDVERSION = 2;
				if(qFromBigEndian(_header.length) == IMDVERSION) _serverEndianess = QSysInfo::BigEndian;
				else if(qFromLittleEndian(_header.length) == IMDVERSION) _serverEndianess = QSysInfo::LittleEndian;
				else {
					protocolError(tr("IMD server handshaking failed. Could not determine endianess of server side."));
					break;
				}
				_isConnected = true;

				// Send IMD_GO message back to server.
				IMDHeader header;
				header.type = qToBigEndian<qint32>(IMD_GO);
				header.length = qToBigEndian<qint32>(sizeof(IMDHeader));
				_socket.write(reinterpret_cast<const char*>(&header), sizeof(header));

				// Send IMD_TRATE message to server.
				header.type = qToBigEndian<qint32>(IMD_TRATE);
				header.length = qToBigEndian<qint32>(transmissionInterval());
				_socket.write(reinterpret_cast<const char*>(&header), sizeof(header));
			}
			else if(_header.type == IMD_FCOORDS) {
				// Receive coordinates.
				qint64 numCoords = qFromBigEndian(_header.length);
				if(numCoords <= 0) {
					protocolError(tr("Invalid number of coordinates: %1").arg(numCoords));
					break;
				}
				// Calculate number of bytes to receive.
				_messageBytesToReceive = numCoords * sizeof(float) * 3;
			}
			else if(_header.type == IMD_ENERGIES) {
				// Calculate number of bytes to receive.
				_messageBytesToReceive = sizeof(IMDEnergies);
			}
			else if(_header.type == IMD_IOERROR) {
				protocolError(tr("Received I/O error signal from server"));
				break;
			}
			else {
				protocolError(tr("Unsupported IMD message type: %1").arg(_header.type));
				break;
			}
		}
		else {
			if(_socket.bytesAvailable() < _messageBytesToReceive) 
				break;	// Wait until more data arrives.

			if(_header.type == IMD_FCOORDS) {
				// Read coordinates data from socket stream.
				qint64 numCoords = qFromBigEndian(_header.length);
				std::vector<Point_3<float>> coords(numCoords);
				OVITO_ASSERT(_messageBytesToReceive == numCoords * sizeof(Point_3<float>));
				if(_socket.read(reinterpret_cast<char*>(coords.data()), _messageBytesToReceive) != _messageBytesToReceive) {
					protocolError(tr("Could not read sufficient data from socket"));
					break;
				}
				_messageBytesToReceive = 0;

				// Convert data array into particle coordinates property.
				_coordinates = ParticlesObject::OOClass().createStandardStorage(numCoords, ParticlesObject::PositionProperty, false);
				std::transform(coords.cbegin(), coords.cend(), PropertyAccess<Point3>(_coordinates).begin(), [](const Point_3<float>& p) { return static_cast<Point3>(p); });

				// Notify pipeline system that this modifier has new results. 
				_numFramesReceived++;
				_modifierStatus = PipelineStatus(PipelineStatus::Success, tr("Connected to IMD server.\n%n frames received.", nullptr, _numFramesReceived));
				if(!_preliminaryUpdateSuspender) {
					// Suppress viewport updates that normally occur when preliminary pipeline updates are available.
					// That's because the user propably likes to see only the final pipeline output when playing 
					// an IMD trajectory pseudo-animation. 
					_preliminaryUpdateSuspender.emplace(dataset()->animationSettings());
				}
				if(!_pipelineUpdatePending) {
					_pipelineUpdatePending = true;
					// Wait until pipeline update that is caurrently in progress has completed before
					// triggering a new pipeline update.
					dataset()->whenSceneReady().finally(executor(), [&](const TaskPtr& task) {
						if(!task->isCanceled()) {
							_pipelineUpdatePending = false;
							notifyTargetChanged();
						}
					});
				}
			}
			else if(_header.type == IMD_ENERGIES) {
				// Ignore message.
			}
			else {
				protocolError(tr("Unsupported IMD message type: %1").arg(_header.type));
				break;
			}
		}
	}
}

/******************************************************************************
* Asks this object to delete itself.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::deleteReferenceObject()
{
	// Automatically disconnect from server when modifier is being deleted.
	disconnectFromServer();

	Modifier::deleteReferenceObject();
}

/******************************************************************************
* Modifies the input data synchronously.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::evaluateSynchronous(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	state.setStatus(_modifierStatus);
	if(!state || !_coordinates) return;

	// Make a modifiable copy of the particles object.
	ParticlesObject* outputParticles = state.expectMutableObject<ParticlesObject>();
	outputParticles->verifyIntegrity();

	// Try to replace particle positions with coordinates received from server.
	PropertyObject* posProperty = outputParticles->expectMutableProperty(ParticlesObject::PositionProperty);
	if(_coordinates->size() != posProperty->size())
		throwException(tr("Number of local particles (%1) does not match number of coordinates received from IMD server (%2).").arg(posProperty->size()).arg(_coordinates->size()));

	posProperty->setStorage(_coordinates);

	// Check if there are any bonds and a simulation cell with periodic boundary conditions.
	if(outputParticles->bonds()) {
		if(const SimulationCellObject* cellObj = state.getObject<SimulationCellObject>()) {
			const SimulationCell cell = cellObj->data();
			if(cell.hasPbc()) {
				if(ConstPropertyAccess<ParticleIndexPair> topologyProperty = outputParticles->bonds()->getProperty(BondsObject::TopologyProperty)) {
					outputParticles->makeBondsMutable();
					ConstPropertyAccess<Point3> positions(posProperty);
					PropertyAccess<Vector3I> periodicImageProperty = outputParticles->bonds()->createProperty(BondsObject::PeriodicImageProperty, true);
					// Recompute PBC vectors of bonds as particle may have moved over arbitrary distances.
					parallelForChunks(topologyProperty.size(), [&](size_t startIndex, size_t count) {
						for(size_t bondIndex = startIndex, endIndex = startIndex+count; bondIndex < endIndex; bondIndex++) {
							size_t index1 = topologyProperty[bondIndex][0];
							size_t index2 = topologyProperty[bondIndex][1];
							if(index1 < positions.size() && index2 < positions.size()) {
								Vector3 delta = cell.absoluteToReduced(positions[index1] - positions[index2]);
								for(size_t dim = 0; dim < 3; dim++) {
									if(cell.hasPbc(dim))
										periodicImageProperty[bondIndex][dim] = (int)std::round(delta[dim]);
								}
							}
						}
					});
				}
			}
		}
	}
}

/******************************************************************************
* Is called when the value of a non-animatable property field of this RefMaker has changed.
******************************************************************************/
void InteractiveMolecularDynamicsModifier::propertyChanged(const PropertyFieldDescriptor& field)
{
	if(field == PROPERTY_FIELD(transmissionInterval)) {
		if(_socket.state() == QAbstractSocket::ConnectedState && _isConnected) {
			// Send IMD_TRATE message to server.
			IMDHeader header;
			header.type = qToBigEndian<qint32>(IMD_TRATE);
			header.length = qToBigEndian<qint32>(transmissionInterval());
			_socket.write(reinterpret_cast<const char*>(&header), sizeof(header));
		}
	}
	else if(field == PROPERTY_FIELD(Modifier::isEnabled)) {
		// Automatically disconnect from server when user disables the modifier.
		if(!isEnabled())
			disconnectFromServer();
	}

	Modifier::propertyChanged(field);
}

}	// End of namespace
}	// End of namespace
