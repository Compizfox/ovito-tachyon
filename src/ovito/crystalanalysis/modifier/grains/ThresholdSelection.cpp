////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2020 Alexander Stukowski
//  Copyright 2020 Peter Mahler Larsen
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

namespace {

FloatType calculate_median(std::vector< FloatType >& data)
{
	size_t n = data.size();
	std::sort(data.begin(), data.end());
	FloatType median = data[n / 2];
	if (n % 2 == 0) {
		median += data[n / 2 - 1];
		median /= 2;
	}

	return median;
}

std::vector< FloatType > theil_sen_estimator(size_t num_samples, std::vector< std::tuple< FloatType, FloatType> >& data,
											  FloatType& gradient, FloatType& intercept)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, data.size() - 1);

	// Fit gradient (median of random gradient samples)
	std::vector< FloatType > gradients;
	for (size_t it=0;it<num_samples;it++) {
		size_t i = dis(gen);
		size_t j = i;
		while (i == j) {
			j = dis(gen);
		}

		auto a = data[i];
		auto b = data[j];
		if (std::get<0>(b) < std::get<0>(a)) {
			std::swap(a, b);
		}

		FloatType dy = std::get<1>(b) - std::get<1>(a);
		FloatType dx = std::get<0>(b) - std::get<0>(a);
		gradients.push_back(dy / dx);
	}

	gradient = calculate_median(gradients);

	// Fit intercept (median of residuals)
	std::vector< FloatType > residuals;
	for (auto point: data) {
		FloatType x = std::get<0>(point);
		FloatType y = std::get<1>(point);
		residuals.push_back(y - gradient * x);
	}

	intercept = calculate_median(residuals);
	for (size_t i=0;i<residuals.size();i++) {
		residuals[i] -= intercept;
	}

	return residuals;
}

} // End of anonymous namespace

/******************************************************************************
* Calculate a threshold suggestion
******************************************************************************/
FloatType GrainSegmentationEngine1::calculate_threshold_suggestion()
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
	auto residuals = theil_sen_estimator(100000, transformed, gradient, intercept);
	for (size_t i=0;i<residuals.size();i++) {
		residuals[i] = fabs(residuals[i]);
	}

	FloatType mean_absolute_deviation = calculate_median(residuals);

    // Select the threshold as the inlier with the largest distance.
	FloatType minSuggestion = 0;
	for(DendrogramNode& node : _dendrogram) {
		FloatType x = log(node.size);
		FloatType y = log(node.distance);

		FloatType prediction = x * gradient + intercept;
		FloatType residual = y - prediction;
		if (residual < 3. * mean_absolute_deviation) {
			minSuggestion = std::max(minSuggestion, log(node.distance));
		}
	}

	// Sort dendrogram entries by distance (undoing the lexicographic sorting performed above).
	boost::sort(_dendrogram, [](const DendrogramNode& a, const DendrogramNode& b) { return a.distance < b.distance; });
	return minSuggestion;
}

}	// End of namespace
}	// End of namespace
