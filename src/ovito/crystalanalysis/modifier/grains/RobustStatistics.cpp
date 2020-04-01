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
* Clustering using pair sampling algorithm.
******************************************************************************/
FloatType GrainSegmentationEngine::calculate_median(std::vector< FloatType >& data)
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

std::vector< FloatType > GrainSegmentationEngine::theil_sen_estimator(size_t num_samples, std::vector< std::tuple< FloatType, FloatType> >& data,
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

}	// End of namespace
}	// End of namespace
