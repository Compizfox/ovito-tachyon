///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2013) Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  OVITO is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __OVITO_NETCDF_IMPORTER_H
#define __OVITO_NETCDF_IMPORTER_H

#include <core/Core.h>
#include <core/gui/properties/PropertiesEditor.h>
#include <plugins/particles/importer/InputColumnMappingDialog.h>
#include <plugins/particles/importer/ParticleImporter.h>

#ifdef MAKING_MODULE_NETCDF
#  define OVITO_NETCDF_EXPORT Q_DECL_EXPORT
#else
#  define OVITO_NETCDF_EXPORT Q_DECL_IMPORT
#endif

namespace NetCDF {

using namespace Ovito;
using namespace Particles;

/**
 * \brief File parser for text-based NetCDF dump simulation files.
 */
class OVITO_NETCDF_EXPORT NetCDFImporter : public ParticleImporter
{
public:

	/// \brief Constructs a new instance of this class.
	Q_INVOKABLE NetCDFImporter(DataSet *dataset) : ParticleImporter(dataset), _useCustomColumnMapping(true) {
		INIT_PROPERTY_FIELD(NetCDFImporter::_useCustomColumnMapping);
		setMultiTimestepFile(true);
	}

	/// \brief Returns the file filter that specifies the files that can be imported by this service.
	/// \return A wild-card pattern that specifies the file types that can be handled by this import class.
	virtual QString fileFilter() override { return "*"; }

	/// \brief Returns the filter description that is displayed in the drop-down box of the file dialog.
	/// \return A string that describes the file format.
	virtual QString fileFilterDescription() override { return tr("NetCDF Files"); }

	/// \brief Checks if the given file has format that can be read by this importer.
	virtual bool checkFileFormat(QIODevice& input, const QUrl& sourceLocation) override;

	/// Returns the title of this object.
	virtual QString objectTitle() override { return tr("NetCDF"); }

	/// \brief Returns the user-defined mapping between data columns in the input file and
	///        the internal particle properties.
	const InputColumnMapping& customColumnMapping() const { return _customColumnMapping; }

	/// \brief Sets the user-defined mapping between data columns in the input file and
	///        the internal particle properties.
	void setCustomColumnMapping(const InputColumnMapping& mapping);

	/// Displays a dialog box that allows the user to edit the custom file column to particle
	/// property mapping.
	void showEditColumnMappingDialog(QWidget* parent = nullptr);

protected:

	/// The format-specific task object that is responsible for reading an input file in the background.
	class OVITO_NETCDF_EXPORT NetCDFImportTask : public ParticleImportTask
	{
	public:

		/// Normal constructor.
		NetCDFImportTask(const LinkedFileImporter::FrameSourceInformation& frame,
				bool useCustomColumnMapping, const InputColumnMapping& customColumnMapping)
			: ParticleImportTask(frame), _parseFileHeaderOnly(false), _useCustomColumnMapping(useCustomColumnMapping), _customColumnMapping(customColumnMapping), _ncIsOpen(false), _ncid(-1) {}

		/// Constructor used when reading only the file header information.
		NetCDFImportTask(const LinkedFileImporter::FrameSourceInformation& frame)
			: ParticleImportTask(frame), _parseFileHeaderOnly(true), _useCustomColumnMapping(false), _ncIsOpen(false), _ncid(-1) {}

		/// Returns the file column mapping used to load the file.
		const InputColumnMapping& columnMapping() const { return _customColumnMapping; }

		/// Is called in the background thread to perform the data file import.
		virtual void load(DataSetContainer& container, FutureInterfaceBase& futureInterface) override;

	protected:

		/// Parses the given input file and stores the data in this container object.
		virtual void parseFile(FutureInterfaceBase& futureInterface, CompressedTextParserStream& stream) override { };

	private:

		/// Is the NetCDF file open?
		bool _ncIsOpen;

		/// NetCDF ids.
		int _ncid, _frame_dim, _atom_dim, _spatial_dim;
		int _cell_spatial_dim, _cell_angular_dim;
		int _cell_origin_var, _cell_lengths_var, _cell_angles_var;
		int _shear_dx_var;

		/// Open NetCDF file, and load additional information
		void openNetCDF(const QString &filename);

		/// Close the current NetCDF file.
		void closeNetCDF();

		bool _parseFileHeaderOnly;
		bool _useCustomColumnMapping;
		InputColumnMapping _customColumnMapping;
	};

protected:

	/// \brief Saves the class' contents to the given stream.
	virtual void saveToStream(ObjectSaveStream& stream) override;

	/// \brief Loads the class' contents from the given stream.
	virtual void loadFromStream(ObjectLoadStream& stream) override;

	/// \brief Creates a copy of this object.
	virtual OORef<RefTarget> clone(bool deepCopy, CloneHelper& cloneHelper) override;

	/// \brief Creates an import task object to read the given frame.
	virtual ImportTaskPtr createImportTask(const FrameSourceInformation& frame) override {
		return std::make_shared<NetCDFImportTask>(frame, _useCustomColumnMapping, _customColumnMapping);
	}

	/// Retrieves the given file in the background and scans it for simulation timesteps.
	virtual void scanMultiTimestepFile(FutureInterface<QVector<LinkedFileImporter::FrameSourceInformation>>& futureInterface, const QUrl sourceUrl);

	/// \brief Guesses the mapping of input file columns to internal particle properties.
	static void mapVariableToColumn(InputColumnMapping &columnMapping, int column, QString name, int dataType);

private:

	/// Controls whether the mapping between input file columns and particle
	/// properties is done automatically or by the user.
	PropertyField<bool> _useCustomColumnMapping;

	/// Stores the user-defined mapping between data columns in the input file and
	/// the internal particle properties.
	InputColumnMapping _customColumnMapping;

	Q_OBJECT
	OVITO_OBJECT

	DECLARE_PROPERTY_FIELD(_useCustomColumnMapping);
};

/**
 * \brief A properties editor for the NetCDFImporter class.
 */
class OVITO_NETCDF_EXPORT NetCDFImporterEditor : public PropertiesEditor
{
public:

	/// Constructor.
	Q_INVOKABLE NetCDFImporterEditor() {}

protected:

	/// Creates the user interface controls for the editor.
	virtual void createUI(const RolloutInsertionParameters& rolloutParams) override;

protected Q_SLOTS:

	/// Is called when the user pressed the "Edit column mapping" button.
	void onEditColumnMapping();

private:

	Q_OBJECT
	OVITO_OBJECT
};


};

#endif // __OVITO_NETCDF_IMPORTER_H
