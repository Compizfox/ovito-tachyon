////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2017 Alexander Stukowski
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


#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include <ovito/crystalanalysis/objects/MicrostructurePhase.h>
#include <ovito/crystalanalysis/data/ClusterGraph.h>
#include <ovito/crystalanalysis/data/DislocationNetwork.h>
#include <ovito/particles/import/ParticleImporter.h>
#include <ovito/particles/import/ParticleFrameData.h>
#include <ovito/mesh/surface/HalfEdgeMesh.h>

namespace Ovito { namespace CrystalAnalysis {

/**
 * \brief Importer for output files generated by the Crystal Analysis Tool.
 */
class OVITO_CRYSTALANALYSIS_EXPORT CAImporter : public ParticleImporter
{
	/// Defines a metaclass specialization for this importer type.
	class OOMetaClass : public ParticleImporter::OOMetaClass
	{
	public:
		/// Inherit standard constructor from base meta class.
		using ParticleImporter::OOMetaClass ::OOMetaClass;

		/// Returns the file filter that specifies the files that can be imported by this service.
		virtual QString fileFilter() const override { return QStringLiteral("*"); }

		/// Returns the filter description that is displayed in the drop-down box of the file dialog.
		virtual QString fileFilterDescription() const override { return tr("Crystal Analysis files"); }

		/// Checks if the given file has format that can be read by this importer.
		virtual bool checkFileFormat(QFileDevice& input, const QUrl& sourceLocation) const override;
	};

	OVITO_CLASS_META(CAImporter, OOMetaClass)
	Q_OBJECT

public:

	/// \brief Constructs a new instance of this class.
	Q_INVOKABLE CAImporter(DataSet* dataset) : ParticleImporter(dataset) {}

	/// Returns the title of this object.
	virtual QString objectTitle() const override { return tr("CA File"); }

	/// Creates an asynchronous loader object that loads the data for the given frame from the external file.
	virtual std::shared_ptr<FileSourceImporter::FrameLoader> createFrameLoader(const Frame& frame, const QString& localFilename) override {
		return std::make_shared<FrameLoader>(frame, localFilename);
	}

	/// Creates an asynchronous frame discovery object that scans the input file for contained animation frames.
	virtual std::shared_ptr<FileSourceImporter::FrameFinder> createFrameFinder(const QUrl& sourceUrl, const QString& localFilename) override {
		return std::make_shared<FrameFinder>(sourceUrl, localFilename);
	}

protected:

	/// The format-specific data holder.
	class CrystalAnalysisFrameData : public ParticleFrameData
	{
	public:

		struct BurgersVectorFamilyInfo {
			int id = 0;
			QString name;
			Vector3 burgersVector = Vector3::Zero();
			Color color = Color(1,1,1);
		};

		struct PatternInfo {
			int id = 0;
			MicrostructurePhase::Dimensionality type = MicrostructurePhase::Dimensionality::Volumetric;
			MicrostructurePhase::CrystalSymmetryClass symmetryType = MicrostructurePhase::CrystalSymmetryClass::CubicSymmetry;
			QString shortName;
			QString longName;
			Color color = Color(1,1,1);
			QVector<BurgersVectorFamilyInfo> burgersVectorFamilies;
		};

	public:

		/// Inherit constructor from base class.
		using ParticleFrameData::ParticleFrameData;

		/// Inserts the loaded data into the provided pipeline state structure. This function is
		/// called by the system from the main thread after the asynchronous loading task has finished.
		virtual OORef<DataCollection> handOver(const DataCollection* existing, bool isNewFile, FileSource* fileSource) override;

		void addPattern(PatternInfo pattern) {
			_patterns.push_back(std::move(pattern));
		}

		Cluster* createCluster(int patternId) {
			return clusterGraph()->createCluster(patternId);
		}

		const std::shared_ptr<ClusterGraph>& clusterGraph() const {
			OVITO_ASSERT(_clusterGraph);
			return _clusterGraph;
		}

		const std::shared_ptr<DislocationNetwork>& dislocations() {
			if(!_dislocations) _dislocations = std::make_shared<DislocationNetwork>(clusterGraph());
			return _dislocations;
		}

		const std::unique_ptr<SurfaceMeshData>& defectSurface() const {
			return _defectSurface;
		}

		void setDefectSurface(std::unique_ptr<SurfaceMeshData> mesh) {
			_defectSurface = std::move(mesh);
		}

	protected:

		/// The structure pattern catalog.
		QVector<PatternInfo> _patterns;

		/// The crystal cluster list.
		ClusterGraphPtr _clusterGraph = std::make_shared<ClusterGraph>();

		/// The dislocation lines.
		std::shared_ptr<DislocationNetwork> _dislocations;

		/// The defect surface mesh.
		std::unique_ptr<SurfaceMeshData> _defectSurface;
	};

	/// The format-specific task object that is responsible for reading an input file in the background.
	class FrameLoader : public FileSourceImporter::FrameLoader
	{
	public:

		/// Inherit constructor from base class.
		using FileSourceImporter::FrameLoader::FrameLoader;

	protected:

		/// Loads the frame data from the given file.
		virtual FrameDataPtr loadFile(QFile& file) override;
	};

	/// The format-specific task object that is responsible for scanning the input file for animation frames.
	class FrameFinder : public FileSourceImporter::FrameFinder
	{
	public:

		/// Inherit constructor from base class.
		using FileSourceImporter::FrameFinder::FrameFinder;

	protected:

		/// Scans the given file for source frames.
		virtual void discoverFramesInFile(QFile& file, const QUrl& sourceUrl, QVector<FileSourceImporter::Frame>& frames) override;
	};
};

}	// End of namespace
}	// End of namespace
