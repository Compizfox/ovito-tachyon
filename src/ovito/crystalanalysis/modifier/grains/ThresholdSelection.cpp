////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2019 Alexander Stukowski
//  Copyright 2019 Peter Mahler Larsen
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

#include <ovito/crystalanalysis/CrystalAnalysis.h>
#include "GrainSegmentationEngine.h"

namespace Ovito { namespace CrystalAnalysis {

/******************************************************************************
* Calculate a threshold suggestion
******************************************************************************/
FloatType GrainSegmentationEngine::calculate_threshold_suggestion()
{
	// Temporarily sort dendrogram entries lexicographically, by (merge size, distance).
    // (we revert this at the end of the function).
	boost::sort(_dendrogram, [](const DendrogramNode& a, const DendrogramNode& b)
	{
		if (a.size == b.size)
			return a.distance < b.distance;
		else
			return a.size < b.size;
	});

	// Transform the data into x = log size, y = log median distance
	size_t i = 0;
	size_t start = 0, currentSize = _dendrogram[0].size;
	std::vector< std::tuple< FloatType, FloatType > > transformed; // median of y values for each dsize
	for(i=0;i<_dendrogram.size();i++) {
		DendrogramNode& node = _dendrogram[i];

		if (node.size != currentSize) {
			size_t count = i - start;
			FloatType median = _dendrogram[start + count / 2].distance;
			if (count % 2 == 0) {
				median += _dendrogram[start + count / 2 - 1].distance;
				median /= 2;
			}
			transformed.push_back(std::make_tuple(log(currentSize), log(median)));

			currentSize = node.size;
			start = i;
		}
	}

	size_t count = i - start;
	FloatType median = _dendrogram[i + count / 2].distance;
	if (count % 2 == 0) {
		median += _dendrogram[i + count / 2 - 1].distance;
		median /= 2;
	}
	transformed.push_back(std::make_tuple(log(currentSize), log(median)));

	// Use Theil-Sen estimator to perform a robust linear regression
	FloatType gradient, intercept;
	auto residuals = GrainSegmentationEngine::theil_sen_estimator(100000, transformed, gradient, intercept);
	for (size_t i=0;i<residuals.size();i++) {
		residuals[i] = fabs(residuals[i]);
	}

    // Calculate a robust estimate of the standard deviation of the residuals
    // c.f. https://en.wikipedia.org/wiki/Median_absolute_deviation#relation_to_standard_deviation
	FloatType sigma = 1.4826 * calculate_median(residuals);

    // Select the threshold as the inlier with the largest distance.
	FloatType minSuggestion = 0;
	for(DendrogramNode& node : _dendrogram) {
		FloatType x = log(node.size);
		FloatType y = log(node.distance);

		FloatType prediction = x * gradient + intercept;
		FloatType residual = y - prediction;
		if (residual < 2. * sigma) {
			minSuggestion = std::max(minSuggestion, log(node.distance));
		}
	}

	// Sort dendrogram entries by distance (undoing the lexicographic sorting performed above).
	boost::sort(_dendrogram, [](const DendrogramNode& a, const DendrogramNode& b) { return a.distance < b.distance; });

#if 0
	// Set a slightly higher threshold, ignoring small merges
	FloatType maxSuggestion = minSuggestion;
	for(DendrogramNode& node : _dendrogram) {
		if (log(node.distance) <= maxSuggestion) continue;
		if (node.size >= _minPlotSize) break;
		maxSuggestion = log(node.distance);
	}
#endif

	return minSuggestion;//(minSuggestion + maxSuggestion) / 2;
}

}	// End of namespace
}	// End of namespace
