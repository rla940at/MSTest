#pragma once

#include "Transformation.h"

#include <unordered_map>

class Figure;
namespace Math{
	std::vector<Polynomial> calculate_Initial_Basis_Function_Set(const Figure& figure, const size_t order);

	std::vector<Polynomial> calculate_Orthonormal_Basis_Function_Set(const Figure& figure, const size_t order);

	size_t combination(const size_t n, const size_t k);

	size_t combination_with_repetition(const size_t n, const size_t k);

	double integrate(const Polynomial& integrand, const QuadratureRule& quadrature_rule);

	double integrate(const Polynomial& integrand, const Figure& figure);

	double inner_product(const Polynomial& f1, const Polynomial& f2, const QuadratureRule& quadrature_rule);

	double inner_product(const Polynomial& f1, const Polynomial& f2, const Figure& figure);

	double L2_norm(const Polynomial& polynomial, const QuadratureRule& quadrature_rule);

	double L2_norm(const Polynomial& polynomial, const Figure& figure);

	std::vector<Polynomial> Gram_Schmidt_Process(const std::vector<Polynomial>& initial_polynomial_set, const QuadratureRule& quadrature_rule);

	std::vector<Polynomial> Gram_Schmidt_Process(const std::vector<Polynomial>& initial_polynomial_set, const Figure& figure);
}


class Figure : public FigureBase
{
private:
	size_t order_;

	FigureType type_;
	
	Transformation transformation_function_;

	std::vector<const MathVector*> consisting_node_set_;

public:
	explicit Figure(const FigureType type, const size_t order, std::vector<const MathVector*>&& consisting_node_set);


	MathVector calculate_Center_Node(void) const {
		return this->transformation_function_(this->calculate_Reference_Center_Node());
	};

	std::vector<std::vector<double>> calculate_Connecitivity(const size_t post_order, const size_t start_index) const;
	
	std::vector<std::vector<size_t>> calculate_Face_Consisting_Node_Index_Family(const std::vector<size_t>& figre_consisting_node_index_set) const;

	std::vector<std::vector<const MathVector*>> calculate_Face_Consisting_Node_Family(void) const;

	std::vector<FigureType> calculate_Face_Figure_Type_Set(void) const;

	MathVector calculate_Normal_Vector(const MathVector& point) const;

	std::vector<MathVector> calculate_Post_Point_Set(const size_t post_order) const;

	std::vector<double> calculate_Projection_Volume_Set(void) const;

	QuadratureRule calculate_Quadrature_Rule(const size_t integrand_order) const;
	
	std::vector<size_t> calculate_Vertex_Node_Index_Set(const std::vector<size_t>& figre_consisting_node_index_set) const;

	std::unordered_map<size_t, std::vector<size_t>> calculate_Vertex_Node_Index_To_Vertex_Simplex_Element_Consisting_Node_Index(const std::vector<size_t>& calculate_Vertex_MathVector_Index_Set) const;
		
	double calculate_Volume(void) const;

	size_t dimension(void) const;

	size_t order(void) const{
		return this->order_; 
	};

	FigureType type(void) const {
		return this->type_;
	};
	
	const Transformation& transformation_function(void) const {
		return this->transformation_function_;
	};

private:
	MathVector calculate_Reference_Center_Node(void) const;

	std::vector<std::vector<size_t>> calculate_Face_Consisting_Node_Index_Order_Family(void) const;

	MathVector calculate_Reference_Normal_Vector(void) const;

	FigureType calculate_Simplex_Figure_Type(void) const;

	std::vector<size_t> calculate_Vertex_Node_Index_Order_Set(void) const;
	
	std::vector<std::vector<size_t>> calculate_Vertex_Simplex_Element_Consisting_Node_Index_Order_Family(void) const;

	bool is_Simplex(void) const;
};


