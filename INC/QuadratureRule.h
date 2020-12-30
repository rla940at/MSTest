#pragma once

#include "MathVector.h"

class QuadratureRule
{
private:
	std::vector<MathVector> quadrature_points_;
	std::vector<double> quadrature_weights_;

public:
	explicit QuadratureRule(void) = default;

	explicit QuadratureRule(const std::vector<MathVector>& quadrature_points, const std::vector<double>& quadrature_weights)
		: quadrature_points_(quadrature_points), quadrature_weights_(quadrature_weights) {};

public:
	const std::vector<MathVector>& quadrature_point_set(void) const {
		return quadrature_points_; 
	};

	const std::vector<double>& quadrature_weight_set(void) const {
		return quadrature_weights_; 
	};
};
