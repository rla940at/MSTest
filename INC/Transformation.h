#pragma once

#include "FigureBase.h"
#include "VectorFunction.h"

class Transformation
{
private:
	VectorFunction<Polynomial> transformation_function_;

	JacobianMatrix<Polynomial> transformation_Jacobian_matrix_;

public:
	explicit Transformation(void) = default;

	explicit Transformation
	(
		const FigureType figure_type,
		const size_t transformation_order,
		const std::vector<MathVector>& original_point_set,
		const std::vector<const MathVector*>& transformed_point_set
	);

	explicit Transformation
	(
		const FigureType figure_type,
		const size_t transformation_order,
		const std::vector<MathVector>& original_point_set,
		const std::vector<MathVector>& transformed_point_set
	);
		

	MathVector operator()(const MathVector& point) const;
	
	std::vector<MathVector> operator()(const std::vector<MathVector>& points) const ;


	double calculate_Transformation_Scale(const FigureType figure_type, const MathVector& point) const;

	RowMajorMatrix calculate_Transformation_Normal_Matrix(const FigureType figure_type, const MathVector& point) const;

private:
	std::vector<Monomial> calculate_Transformation_Monomial_Set(const FigureType figure_type, const size_t transformation_order) const;	
};