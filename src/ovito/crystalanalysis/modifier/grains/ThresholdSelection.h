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

std::vector< FloatType > theil_sen_estimator(size_t num_samples,
                                             std::vector< FloatType >& xs, std::vector< FloatType >& ys,
											 FloatType& gradient, FloatType& intercept)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, xs.size() - 1);

	// Fit gradient (median of random gradient samples)
	std::vector< FloatType > gradients;
	for (size_t it=0;it<num_samples;it++) {
		size_t i = dis(gen);
		size_t j = i;
		while (i == j) {
			j = dis(gen);
		}

		if (xs[j] < xs[i]) {
			std::swap(i, j);
		}

		FloatType dx = xs[j] - xs[i];
		FloatType dy = ys[j] - ys[i];
		gradients.push_back(dy / dx);
	}

	gradient = calculate_median(gradients);

	// Fit intercept (median of residuals)
	std::vector< FloatType > residuals;
    for (size_t i=0;i<xs.size();i++) {
		residuals.push_back(ys[i] - gradient * xs[i]);
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
namespace ThresholdSelection {

class Regressor
{
public:

    FloatType gradient = 0, intercept = 0;
    FloatType mean_absolute_deviation = 0;
    std::vector< FloatType > residuals;
    std::vector< FloatType > dsize;
    std::vector< FloatType > medianDistance;     // median of y values for each dsize

    Regressor(std::vector<GrainSegmentationEngine1::DendrogramNode>& dendrogram) {
        if (dendrogram.size() == 0)
            return;

	    // Transform the data into x = log size, y = log median distance
	    size_t start = 0, currentSize = dendrogram[0].size;
	    for(size_t i=0;i<dendrogram.size();i++) {
		    auto node = dendrogram[i];

		    if (node.size != currentSize) {
			    size_t count = i - start;
			    FloatType median = dendrogram[start + count / 2].distance;
			    if (count % 2 == 0) {
				    median += dendrogram[start + count / 2 - 1].distance;
				    median /= 2;
			    }
			    dsize.push_back(log(currentSize));
                medianDistance.push_back(log(median));

			    currentSize = node.size;
			    start = i;
		    }
	    }

	    size_t count = dendrogram.size() - start;
	    FloatType median = dendrogram[start + count / 2].distance;
	    if (count % 2 == 0) {
		    median += dendrogram[start + count / 2 - 1].distance;
		    median /= 2;
	    }
	    dsize.push_back(log(currentSize));
        medianDistance.push_back(log(median));

	    // Use Theil-Sen estimator to perform a robust linear regression
	    residuals = theil_sen_estimator(100000, dsize, medianDistance, gradient, intercept);
	    for (size_t i=0;i<residuals.size();i++) {
		    residuals[i] = fabs(residuals[i]);
	    }

	    mean_absolute_deviation = calculate_median(residuals);
    }

    FloatType calculate_threshold(std::vector<GrainSegmentationEngine1::DendrogramNode>& dendrogram, FloatType cutoff) {

        // Select the threshold as the inlier with the largest distance.
	    FloatType threshold = 0;
	    for(auto node : dendrogram) {
		    FloatType x = log(node.size);
		    FloatType y = log(node.distance);

		    FloatType prediction = x * gradient + intercept;
		    FloatType residual = y - prediction;

		    if (residual < cutoff * mean_absolute_deviation) {
			    threshold = std::max(threshold, y);
		    }
	    }

	    return threshold;
    }
};

}	// End of namespace
}	// End of namespace
}	// End of namespace
