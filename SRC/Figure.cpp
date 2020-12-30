#include "../INC/Figure.h"

namespace Math
{
	size_t Math::combination(const size_t n, const size_t k){ //calculate nCk			
		constexpr size_t allowable_maximum = 32;
		if (n > allowable_maximum)
			FATAL_ERROR("n(" << n << ") exceed allowable maximum(" << allowable_maximum << ")");
		if (n < k)
			FATAL_ERROR("n(" << n << ") should be bigger than k(" << k << ")");

		if (n == k || k == 0)
			return 1;
		else
			return combination(n - 1, k - 1) + combination(n - 1, k);
	}

	size_t Math::combination_with_repetition(const size_t n, const size_t k){ //calculate nHk		
		constexpr size_t allowable_maximum = 32;
		if (n + k - 1 > allowable_maximum)
			FATAL_ERROR("n+k-1(" << n + k - 1 << ") exceed allowable maximum(" << allowable_maximum << ")");

		return Math::combination(n + k - 1, k);
	}

	std::vector<Polynomial> calculate_Initial_Basis_Function_Set(const Figure& figure, const size_t order){
		const auto figure_dimension = figure.dimension();
		const auto center_node = figure.calculate_Center_Node();

		const auto num_basis = Math::combination_with_repetition(1 + figure_dimension, order);
		std::vector<Polynomial> initial_basis_set;
		initial_basis_set.reserve(num_basis);

		if (figure_dimension == 2){
			// 1 (x - x_c) ( y - y_c )  ...
			for (size_t a = 0; a <= order; ++a)
			for (size_t b = 0; b <= a; ++b)
				initial_basis_set.emplace_back(Monomial{ a - b, b });

			for (auto& basis : initial_basis_set)
				Math::translate(basis, center_node);
		}
		else
			FATAL_TYPE_ERROR;

		return initial_basis_set;
	}

	std::vector<Polynomial> calculate_Orthonormal_Basis_Function_Set(const Figure& figure, const size_t order){
		const auto initial_basis_function_set = Math::calculate_Initial_Basis_Function_Set(figure, order);

		const auto integrand_order = 2 * order;
		const auto quadrature_rule = figure.calculate_Quadrature_Rule(integrand_order);

		return Math::Gram_Schmidt_Process(initial_basis_function_set, quadrature_rule);
	}

	double integrate(const Polynomial& integrand, const QuadratureRule& quadrature_rule)	{
		const auto& QP_set = quadrature_rule.quadrature_point_set();
		const auto& QW_set = quadrature_rule.quadrature_weight_set();

		double result = 0.0; 
		for (size_t i = 0; i < QP_set.size(); ++i)
			result += integrand(QP_set[i]) * QW_set[i];

		return result;
	}

	double integrate(const Polynomial& integrand, const Figure& figure){
		const auto quadrature_rule = figure.calculate_Quadrature_Rule(integrand.order());
		return Math::integrate(integrand, quadrature_rule);
	}

	double inner_product(const Polynomial& f1, const Polynomial& f2, const QuadratureRule& quadrature_rule){
		return Math::integrate(f1 * f2, quadrature_rule);
	}

	double inner_product(const Polynomial& f1, const Polynomial& f2, const Figure& geometry){
		const auto quadrature_rule = geometry.calculate_Quadrature_Rule(f1.order() + f2.order());
		return Math::inner_product(f1, f2, quadrature_rule);
	}

	double L2_norm(const Polynomial& function, const QuadratureRule& quadrature_rule){
		return std::sqrt(Math::inner_product(function, function, quadrature_rule));
	}

	double L2_norm(const Polynomial& polynomial, const Figure& geometry){
		const auto quadrature_rule = geometry.calculate_Quadrature_Rule(polynomial.order() * 2);
		return Math::L2_norm(polynomial, quadrature_rule);
	}

	std::vector<Polynomial> Gram_Schmidt_Process(const std::vector<Polynomial>& initial_polynomial_set, const QuadratureRule& quadrature_rule)	{
		auto normalized_polynomial_set = initial_polynomial_set;
		for (size_t i = 0; i < initial_polynomial_set.size(); ++i){
			for (size_t j = 0; j < i; ++j)
				normalized_polynomial_set[i] -= normalized_polynomial_set[j] * Math::inner_product(normalized_polynomial_set[i], normalized_polynomial_set[j], quadrature_rule);

			normalized_polynomial_set[i] *= 1.0 / Math::L2_norm(normalized_polynomial_set[i], quadrature_rule);
		}

		return normalized_polynomial_set;
	}

	std::vector<Polynomial> Gram_Schmidt_Process(const std::vector<Polynomial>& initial_polynomial_set, const Figure& geometry)	{
		std::vector<size_t> order_set;
		order_set.reserve(initial_polynomial_set.size());

		for (const auto& polynomial : initial_polynomial_set)
			order_set.emplace_back(polynomial.order());

		const auto maximum_order = *std::max_element(order_set.begin(), order_set.end());
		const auto integrand_order = maximum_order * 2;

		const auto quadrature_rule = geometry.calculate_Quadrature_Rule(integrand_order);
		return Math::Gram_Schmidt_Process(initial_polynomial_set, quadrature_rule);
	}
}


Figure::Figure(const FigureType type, const size_t order, std::vector<const MathVector*>&& consisting_node_set)
	: type_(type), order_(order), consisting_node_set_(std::move(consisting_node_set))
{
	const auto& reference_transformation_point_set = FigureBase::ReferenceTransfomrationPointSet::get(this->type_, this->order_);

	this->transformation_function_ = Transformation(this->type_, this->order_, reference_transformation_point_set, this->consisting_node_set_);
};


std::vector<std::vector<double>> Figure::calculate_Connecitivity(const size_t post_order, const size_t start_index) const{
	const auto& reference_connectivity = FigureBase::ReferenceConnectivity::get(this->type_, post_order);

	auto connectivity = reference_connectivity;

	for (auto& simplex_consisting_node_index_set : connectivity)
	for (auto& node_index : simplex_consisting_node_index_set)
		node_index += start_index;

	return connectivity;
}

std::vector<std::vector<size_t>> Figure::calculate_Face_Consisting_Node_Index_Family(const std::vector<size_t>& figure_consisting_node_index_set) const{
	const auto face_consisting_node_index_order_family = this->calculate_Face_Consisting_Node_Index_Order_Family();
	const auto num_face = face_consisting_node_index_order_family.size();

	std::vector<std::vector<size_t>> face_consisting_node_index_family(num_face);
	for (size_t i = 0; i < num_face; ++i)
	{
		const auto num_consisting_node = face_consisting_node_index_order_family[i].size();
		face_consisting_node_index_family[i].resize(num_consisting_node);

		for (size_t j = 0; j < num_consisting_node; ++j)
			face_consisting_node_index_family[i][j] = figure_consisting_node_index_set[face_consisting_node_index_order_family[i][j]];
	}

	return face_consisting_node_index_family;
}

std::vector<std::vector<const MathVector*>> Figure::calculate_Face_Consisting_Node_Family(void) const{
	const auto face_consisting_node_index_order_family = this->calculate_Face_Consisting_Node_Index_Order_Family();
	const auto num_face = face_consisting_node_index_order_family.size();

	std::vector<std::vector<const MathVector*>> face_consisting_node_family(num_face);	
	for (size_t i = 0; i < num_face; ++i) {
		const auto num_consisting_node = face_consisting_node_index_order_family[i].size();
		face_consisting_node_family.resize(num_consisting_node);

		for (size_t j = 0; j < num_consisting_node; ++j)
			face_consisting_node_family[i][j] = this->consisting_node_set_[face_consisting_node_index_order_family[i][j]];
	}

	return face_consisting_node_family;
}

MathVector Figure::calculate_Normal_Vector(const MathVector& position_vector) const{
	const auto reference_normal_vector = this->calculate_Reference_Normal_Vector();
	const auto transformation_normal_matrix = this->transformation_function_.calculate_Transformation_Normal_Matrix(this->type_, position_vector);
	return Math::normalize(transformation_normal_matrix * reference_normal_vector);
}

std::vector<MathVector> Figure::calculate_Post_Point_Set(const size_t post_order) const{
	const auto& reference_post_point_set = FigureBase::ReferencePostPointSet::get(this->type_, post_order);
	return this->transformation_function_(reference_post_point_set);
}

std::vector<double> Figure::calculate_Projection_Volume_Set(void) const{
	switch (this->dimension())
	{
	case 2:	{
		double projection_volume_on_x0_axis = 0.0;
		double projection_volume_on_x1_axis = 0.0;

		for (const auto& face_consisting_node_set : this->calculate_Face_Consisting_Node_Family())
		{
			const auto& start_node = *face_consisting_node_set[0];
			const auto& end_node = *face_consisting_node_set[1];
			const auto node_to_node = end_node - start_node;

			projection_volume_on_x0_axis += std::abs(node_to_node[0]);
			projection_volume_on_x1_axis += std::abs(node_to_node[1]);
		}

		return { 0.5 * projection_volume_on_x0_axis, 0.5 * projection_volume_on_x1_axis };
	}
	default:
		FATAL_TYPE_ERROR;
		return std::vector<double>();
	}
}

QuadratureRule Figure::calculate_Quadrature_Rule(const size_t integrand_order) const{
	const auto& reference_quadrature_rule = FigureBase::ReferenceQuadratureRule::get(this->type_, integrand_order);
	const auto& reference_QP = reference_quadrature_rule.quadrature_point_set();
	const auto& reference_QW = reference_quadrature_rule.quadrature_weight_set();

	const auto transformed_QP = this->transformation_function_(reference_QP);

	const auto num_QP = reference_QP.size();
	std::vector<double> transformed_QW(num_QP);
	for (size_t ipoint = 0; ipoint < num_QP; ++ipoint)
	{
		const auto& node = reference_QP[ipoint];
		const auto trasnformation_scale = this->transformation_function_.calculate_Transformation_Scale(this->type_, node);

		transformed_QW[ipoint] = reference_QW[ipoint] * trasnformation_scale;
	}

	return QuadratureRule(transformed_QP, transformed_QW);
}

std::vector<size_t> Figure::calculate_Vertex_Node_Index_Set(const std::vector<size_t>& figre_consisting_node_index_set) const{
	const auto vertex_node_index_order_set = this->calculate_Vertex_Node_Index_Order_Set();
	const auto num_vertex = vertex_node_index_order_set.size();

	std::vector<size_t> vertex_node_index_set(num_vertex);
	for (size_t i = 0; i < num_vertex; ++i)
		vertex_node_index_set[i] = figre_consisting_node_index_set[vertex_node_index_order_set[i]];

	return vertex_node_index_set;
}

std::unordered_map<size_t, std::vector<size_t>> Figure::calculate_Vertex_Node_Index_To_Vertex_Simplex_Element_Consisting_Node_Index(const std::vector<size_t>& vertex_node_index_set) const{
	const auto vertex_simplex_element_consisting_node_index_order_family = this->calculate_Vertex_Simplex_Element_Consisting_Node_Index_Order_Family();
	const auto num_vertex = vertex_node_index_set.size();

	std::unordered_map<size_t, std::vector<size_t>> calculate_vertex_node_index_to_vertex_simplex_element_consisting_node_index;
	calculate_vertex_node_index_to_vertex_simplex_element_consisting_node_index.reserve(num_vertex);
	for (size_t i = 0; i < num_vertex; ++i)
	{
		const auto vertex_node_index = vertex_node_index_set[i];
		const auto num_consisting_node = vertex_simplex_element_consisting_node_index_order_family[i].size();

		std::vector<size_t> simplex_element_consisting_node_index(num_consisting_node);
		for (size_t j = 0; j < num_consisting_node; ++j)
			simplex_element_consisting_node_index[j] = vertex_node_index_set[vertex_simplex_element_consisting_node_index_order_family[i][j]];

		calculate_vertex_node_index_to_vertex_simplex_element_consisting_node_index.emplace(vertex_node_index, std::move(simplex_element_consisting_node_index));
	}

	return calculate_vertex_node_index_to_vertex_simplex_element_consisting_node_index;
}

double Figure::calculate_Volume(void) const{
	const size_t integrand_order = 0;	// int(1) = volume

	const auto calculate_Quadrature_Rule = this->calculate_Quadrature_Rule(integrand_order);
	const auto& QW_set = calculate_Quadrature_Rule.quadrature_weight_set();

	double calculate_Volume = 0;
	for (const double QW : QW_set)
		calculate_Volume += QW;

	return calculate_Volume;
}

size_t Figure::dimension(void) const{
	switch (this->type_)
	{
	case FigureType::Line:				return 1;
	case FigureType::Triangle:
	case FigureType::Quadrilateral:		return 2;
	default:
		FATAL_TYPE_ERROR;
		return NULL;
	}
}

MathVector Figure::calculate_Reference_Center_Node(void) const{
	switch (this->type_){
	case FigureType::Line:				return { 0, 0, 0 };
	case FigureType::Triangle:			return { -1.0 / 3.0, -1.0 / 3.0, 0 };
	case FigureType::Quadrilateral:		return { 0, 0, 0 };
	default:
		FATAL_TYPE_ERROR;
		return MathVector();
	}
}

std::vector<std::vector<size_t>> Figure::calculate_Face_Consisting_Node_Index_Order_Family(void) const{
	// it tells index order of i - th face consisting node at cell consisting node index

	std::vector<std::vector<size_t>> face_node_index_order;
	switch (this->type_){
	case FigureType::Line:{
		// 0 式式式式 1

		constexpr size_t num_face = 2;

		const std::vector<size_t> first_face_node_index = { 0 };
		const std::vector<size_t> second_face_node_index = { 1 };

		face_node_index_order.resize(num_face);
		face_node_index_order[0] = first_face_node_index;
		face_node_index_order[1] = second_face_node_index;

		break;
	}
	case FigureType::Triangle:{

		//  2
		//  弛 \
		//	弛  \
		//  0式式式1

		constexpr size_t num_face = 3;

		const std::vector<size_t> first_face_node_index = { 0,1 };
		const std::vector<size_t> second_face_node_index = { 1,2 };
		const std::vector<size_t> third_face_node_index = { 2,0 };

		face_node_index_order.resize(num_face);
		face_node_index_order[0] = first_face_node_index;
		face_node_index_order[1] = second_face_node_index;
		face_node_index_order[2] = third_face_node_index;

		//if (element_order > 1)
		//{
		//	const size_t num_additional_point = element_order - 1;
		//	
		//	size_t index = num_face;
		//	for (size_t iface = 0; iface < num_face; ++iface)
		//		for (size_t ipoint = 0; ipoint < num_additional_point; ++ipoint)
		//			face_node_index_order[iface].push_back(index++);
		//}

		break;
	}
	case FigureType::Quadrilateral:{

		//  3式式式式式2
		//  弛     弛
		//  0式式式式式1

		const size_t num_face = 4;

		const std::vector<size_t> first_face_node_index = { 0,1 };
		const std::vector<size_t> second_face_node_index = { 1,2 };
		const std::vector<size_t> third_face_node_index = { 2,3 };
		const std::vector<size_t> forth_face_node_index = { 3,0 };

		face_node_index_order.resize(num_face);
		face_node_index_order[0] = first_face_node_index;
		face_node_index_order[1] = second_face_node_index;
		face_node_index_order[2] = third_face_node_index;
		face_node_index_order[3] = forth_face_node_index;

		//if (element_order > 1)
		//{
		//	const size_t num_additional_point = element_order - 1;

		//	size_t index = num_face;
		//	for (size_t iface = 0; iface < num_face; ++iface)
		//		for (size_t ipoint = 0; ipoint<num_additional_point; ++ipoint)
		//			face_node_index_order[iface].push_back(index++);
		//}

		break;
	}
	default:
		FATAL_TYPE_ERROR;
	}

	return face_node_index_order;
}

std::vector<FigureType> Figure::calculate_Face_Figure_Type_Set(void) const{
	std::vector<FigureType> face_figure_type_set;
	switch (this->type_){
	case FigureType::Line:
	{
		// 0 式式式式 1

		constexpr size_t num_face = 2;

		face_figure_type_set.resize(num_face);
		face_figure_type_set[0] = FigureType::Point;
		face_figure_type_set[1] = FigureType::Point;
		break;
	}
	case FigureType::Triangle:{
		//  弛\
		//  2 1
		//	弛  \
		//  戌式0式式

		constexpr size_t num_face = 3;

		face_figure_type_set.resize(num_face);
		face_figure_type_set[0] = FigureType::Line;
		face_figure_type_set[1] = FigureType::Line;
		face_figure_type_set[2] = FigureType::Line;
		break;
	}
	case FigureType::Quadrilateral:{
		//  忙式 2 式忖
		//  3     1
		//  戌式 0 式戎

		constexpr size_t num_face = 4;

		face_figure_type_set.resize(num_face);
		face_figure_type_set[0] = FigureType::Line;
		face_figure_type_set[1] = FigureType::Line;
		face_figure_type_set[2] = FigureType::Line;
		face_figure_type_set[3] = FigureType::Line;
		break;
	}
	default:
		FATAL_TYPE_ERROR;
	}

	return face_figure_type_set;
}

MathVector Figure::calculate_Reference_Normal_Vector(void) const{
	switch (this->type_){
	case FigureType::Line:				return { 0,1 };
	case FigureType::Triangle:
	case FigureType::Quadrilateral:		return { 0,0,1 };
	default:
		FATAL_TYPE_ERROR;
		return MathVector();
	}
}

FigureType Figure::calculate_Simplex_Figure_Type(void) const{
	if (this->dimension() == 2)
		return FigureType::Triangle;
	else{
		FATAL_TYPE_ERROR;
		return FigureType::NotInList;
	}
}

std::vector<size_t> Figure::calculate_Vertex_Node_Index_Order_Set(void) const{
	switch (this->type_){
	case FigureType::Line:{
		// 0 式式式式 1

		return { 0,1 };
	}
	case FigureType::Triangle: {
		//  2
		//  弛 \
		//	弛  \
		//  0式式式1

		return { 0,1,2 }; 
	}
	case FigureType::Quadrilateral:{
		//  3式式式式式2
		//  弛     弛
		//  0式式式式式1

		return { 0,1,2,3 };
	}

	default:
		FATAL_TYPE_ERROR;
		return { 0 };
	}
}

std::vector<std::vector<size_t>> Figure::calculate_Vertex_Simplex_Element_Consisting_Node_Index_Order_Family(void) const{
	//For hMLP_BD Limiter

	if (this->order_ > 1)
		FATAL_ERROR("hMLP_BD doesn't work on high order mesh!");

	std::vector<std::vector<size_t>> vertex_simplex_space_node_index_order_family;
	switch (this->type_){
	case FigureType::Quadrilateral:
	{
		//  3式式式式式2
		//  弛     弛
		//  0式式式式式1

		constexpr size_t num_simplex = 4;

		const std::vector<size_t> first_simplex_node_index = { 0,1,3 };
		const std::vector<size_t> second_simplex_node_index = { 1,2,0 };
		const std::vector<size_t> third_simplex_node_index = { 2,3,1 };
		const std::vector<size_t> fourth_simplex_node_index = { 3,0,2 };

		vertex_simplex_space_node_index_order_family.reserve(num_simplex);
		vertex_simplex_space_node_index_order_family.emplace_back(first_simplex_node_index);
		vertex_simplex_space_node_index_order_family.emplace_back(second_simplex_node_index);
		vertex_simplex_space_node_index_order_family.emplace_back(third_simplex_node_index);
		vertex_simplex_space_node_index_order_family.emplace_back(fourth_simplex_node_index);

		break;
	}
	default:
		FATAL_TYPE_ERROR;
	}

	return vertex_simplex_space_node_index_order_family;
}

bool Figure::is_Simplex(void) const{
	if (this->type_ == FigureType::Triangle) // (figure_type == ReferenceFigureType::Triangle ||figure_type == ReferenceFigureType::Tetrahedral)
		return true;
	else
		return false;
}