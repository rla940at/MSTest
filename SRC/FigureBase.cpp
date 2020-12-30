#include "../INC/FigureBase.h"

namespace Editor {
	std::string to_String(const FigureType figure_type) {
		switch (figure_type){
		case FigureType::Line:				
			return "Line";
		case FigureType::Triangle:			
			return "Triangle";
		case FigureType::Quadrilateral:		
			return "Quadrilateral";
		default:
			FATAL_TYPE_ERROR;
			return "NotInList";
		}
	}
}


std::map<std::pair<FigureType, size_t>, std::vector<MathVector>> FigureBase::ReferenceTransfomrationPointSet::key_to_value_;

const std::vector<MathVector>& FigureBase::ReferenceTransfomrationPointSet::get(const FigureType figure_type, const size_t transformation_order){
	const auto key = std::make_pair(figure_type, transformation_order);
	if (FigureBase::ReferenceTransfomrationPointSet::key_to_value_.find(key) == FigureBase::ReferenceTransfomrationPointSet::key_to_value_.end())
		FigureBase::ReferenceTransfomrationPointSet::key_to_value_.emplace(key, add(figure_type, transformation_order));

	return FigureBase::ReferenceTransfomrationPointSet::key_to_value_.at(key);
}

std::vector<MathVector> FigureBase::ReferenceTransfomrationPointSet::add(const FigureType figure_type, const size_t transformation_order){
	const auto figure_type_name = Editor::to_String(figure_type);

	if (transformation_order > FigureBase::ReferenceTransfomrationPointSet::support_element_order(figure_type))
		FATAL_ERROR(figure_type_name << " figure doesn't support P" << transformation_order << " transformation");

	std::string read_file_path;
	read_file_path << "RSC/ReferenceFigures/" << figure_type_name << "_P" << transformation_order << ".msh";

	Text reference_transformation_point_set_text(read_file_path);
	Editor::remove(reference_transformation_point_set_text, "");

	return FigureBase::ReferenceTransfomrationPointSet::Text_To_Point_Set(reference_transformation_point_set_text);
}

size_t FigureBase::ReferenceTransfomrationPointSet::support_element_order(const FigureType figure_type){
	switch (figure_type){
	case FigureType::Line: 
		return 6;
	case FigureType::Triangle: 
		return 5;
	case FigureType::Quadrilateral: 
		return 6;
	default:
		FATAL_TYPE_ERROR;
		return NULL;
	}
}

std::vector<MathVector> FigureBase::ReferenceTransfomrationPointSet::Text_To_Point_Set(Text& reference_transformation_point_set_text){
	const auto num_nodes = StringEditor::toValue<size_t>(reference_transformation_point_set_text.front());
	reference_transformation_point_set_text.erase(reference_transformation_point_set_text.begin());

	std::vector<MathVector> point_set;
	point_set.reserve(num_nodes);
	for (const auto& sentence : reference_transformation_point_set_text){
		const char delimiter = ' ';		
		const auto parsed_str = StringEditor::parse(sentence, delimiter);
		auto point_value = Text(parsed_str).to_Value_Set<double>();

		point_value.erase(point_value.begin());
		point_set.emplace_back(std::move(point_value));
	}

	return point_set;
}


std::map<std::pair<FigureType, size_t>, QuadratureRule> FigureBase::ReferenceQuadratureRule::key_to_value_;

const QuadratureRule& FigureBase::ReferenceQuadratureRule::get(const FigureType figure_type, const size_t transformation_order){
	const auto key = std::make_pair(figure_type, transformation_order);
	if (FigureBase::ReferenceQuadratureRule::key_to_value_.find(key) == FigureBase::ReferenceQuadratureRule::key_to_value_.end())
		FigureBase::ReferenceQuadratureRule::key_to_value_.emplace(key, FigureBase::ReferenceQuadratureRule::add(figure_type, transformation_order));
	
	return FigureBase::ReferenceQuadratureRule::key_to_value_.at(key);
}

QuadratureRule FigureBase::ReferenceQuadratureRule::add(const FigureType figure_type, const size_t integrand_order){
	const auto required_order = FigureBase::ReferenceQuadratureRule::calculate_Required_Order(figure_type, integrand_order);
	const auto num_required_point = FigureBase::ReferenceQuadratureRule::calculate_Num_Required_Point(figure_type, required_order);

	std::string read_file_path;
	read_file_path << "RSC/Quadrature/Standard/" << Editor::to_String(figure_type) << "/P" << required_order << "_n" << num_required_point << ".txt";

	Text quadrature_text(read_file_path);
	return FigureBase::ReferenceQuadratureRule::Text_To_Quadrature_Rule(quadrature_text);
}

size_t FigureBase::ReferenceQuadratureRule::calculate_Num_Required_Point(const FigureType figure_type, const size_t required_order){
	switch (figure_type){
	case FigureType::Line:				
		return static_cast<size_t>((required_order + 1) * 0.5);
	case FigureType::Triangle:			
		return static_cast<size_t>((required_order * 0.5 + 1) * (required_order * 0.5 + 1));
	case FigureType::Quadrilateral:		
		return static_cast<size_t>((required_order + 1) * (required_order + 1) * 0.25);
	default:
		FATAL_TYPE_ERROR;
		return NULL;
	}
}

size_t FigureBase::ReferenceQuadratureRule::calculate_Required_Order(const FigureType figure_type, const size_t integrand_order){
	const auto remain = integrand_order % 2;

	switch (figure_type){
	case FigureType::Line:
	case FigureType::Quadrilateral:{
		if (remain == 0)
			return integrand_order + 1;
		else
			return integrand_order;
	}
	case FigureType::Triangle:{
		if (remain == 0)
			return integrand_order;
		else
			return integrand_order + 1;
	}
	default:
		FATAL_TYPE_ERROR;
		return NULL;
	}
}

QuadratureRule FigureBase::ReferenceQuadratureRule::Text_To_Quadrature_Rule(const Text& quadrature_text){
	const auto num_quadrature = quadrature_text.size();

	std::vector<MathVector> quadrature_point_set;
	std::vector<double> quadrature_weight_set;
	quadrature_point_set.reserve(num_quadrature);
	quadrature_weight_set.reserve(num_quadrature);

	for (const auto& sentence : quadrature_text)	{
		const char denominator = ' ';		
		const auto parsed_sentence = StringEditor::parse(sentence, denominator);

		auto quadrature_value_set = Text(parsed_sentence).to_Value_Set<double>();

		quadrature_weight_set.emplace_back(quadrature_value_set.back());
		quadrature_value_set.pop_back();
		quadrature_point_set.emplace_back(std::move(quadrature_value_set));
	}

	return QuadratureRule(quadrature_point_set, quadrature_weight_set);
}


std::map<std::pair<FigureType, size_t>, std::vector<MathVector>> FigureBase::ReferencePostPointSet::key_to_value_;

const std::vector<MathVector>& FigureBase::ReferencePostPointSet::get(const FigureType figure_type, const size_t post_order){
	const auto key = std::make_pair(figure_type, post_order);
	if (FigureBase::ReferencePostPointSet::key_to_value_.find(key) == FigureBase::ReferencePostPointSet::key_to_value_.end())
		FigureBase::ReferencePostPointSet::key_to_value_.emplace(key, add(figure_type, post_order));

	return FigureBase::ReferencePostPointSet::key_to_value_.at(key);
}

std::vector<MathVector> FigureBase::ReferencePostPointSet::add(const FigureType figure_type, const size_t post_order){
	std::vector<MathVector> reference_post_point_set;
	switch (figure_type){
	case FigureType::Triangle:{
		const size_t num_reference_post_point = static_cast<size_t>((post_order + 1) * (post_order + 2) * 0.5);
		reference_post_point_set.reserve(num_reference_post_point);

		const double delta = 2.0 / post_order;

		const double X0_start_coord = -1.0;
		const double X1_start_coord = -1.0;

		for (size_t j = 0; j <= post_order; ++j)
		for (size_t i = 0; i <= post_order - j; ++i)
		{
			const double X0_coord = X0_start_coord + delta * i;
			const double X1_coord = X1_start_coord + delta * j;

			reference_post_point_set.emplace_back(MathVector{ X0_coord, X1_coord });
		}
		break;
	}
	case FigureType::Quadrilateral:{
		const size_t num_reference_post_point = (post_order + 1) * (post_order + 1);
		reference_post_point_set.reserve(num_reference_post_point);

		const double delta = 2.0 / post_order;

		const double X0_start_coord = -1.0;
		const double X1_start_coord = -1.0;

		for (size_t j = 0; j <= post_order; ++j)
		for (size_t i = 0; i <= post_order; ++i)
		{
			const double X0_coord = X0_start_coord + delta * i;
			const double X1_coord = X1_start_coord + delta * j;

			reference_post_point_set.emplace_back(MathVector{ X0_coord, X1_coord });
		}
		break;
	}
	default:
		FATAL_TYPE_ERROR;
		break;
	}

	return reference_post_point_set;
}


std::map<std::pair<FigureType, size_t>, std::vector<std::vector<double>>> FigureBase::ReferenceConnectivity::key_to_value_;

const std::vector<std::vector<double>>& FigureBase::ReferenceConnectivity::get(const FigureType figure_type, const size_t post_order){
	const auto key = std::make_pair(figure_type, post_order);
	if (FigureBase::ReferenceConnectivity::key_to_value_.find(key) == FigureBase::ReferenceConnectivity::key_to_value_.end())
		FigureBase::ReferenceConnectivity::key_to_value_.emplace(key, add(figure_type, post_order));

	return FigureBase::ReferenceConnectivity::key_to_value_.at(key);
}

std::vector<std::vector<double>> FigureBase::ReferenceConnectivity::add(const FigureType figure_type, const size_t post_order){
	std::vector<std::vector<double>> simplex_decomposed_element_consisting_node_index_set;

	switch (figure_type){
	case FigureType::Triangle:{
		const size_t num_simplex = post_order * post_order;
		simplex_decomposed_element_consisting_node_index_set.resize(num_simplex);

		constexpr size_t num_simplex_consisting_node = 3;

		size_t isimplex = 0;
		for (size_t j = 0; j < post_order; j++)
		for (size_t i = 0; i < post_order - j; i++)
		{
			//    b  --------- b+1
			//    |            |
			//    |            |  
			//    a  --------- a+1

			// recurrence relation
			// a(n+1) = a(n) + (order + 2) - n
			//=> a(n) = a(0) + (order + 2) * n - n * (n+1) * 0.5
			const double a_point_index = i + 1 + (post_order + 2) * j - j * (j + 1) * 0.5;
			const double b_point_index = i + 1 + (post_order + 2) * (j + 1) - (j + 1) * (j + 2) * 0.5;

			simplex_decomposed_element_consisting_node_index_set[isimplex].resize(num_simplex_consisting_node);
			simplex_decomposed_element_consisting_node_index_set[isimplex][0] = a_point_index;
			simplex_decomposed_element_consisting_node_index_set[isimplex][1] = a_point_index + 1;
			simplex_decomposed_element_consisting_node_index_set[isimplex][2] = b_point_index;
			isimplex++;

			if (i == post_order - j - 1) continue;

			simplex_decomposed_element_consisting_node_index_set[isimplex].resize(num_simplex_consisting_node);
			simplex_decomposed_element_consisting_node_index_set[isimplex][0] = a_point_index + 1;
			simplex_decomposed_element_consisting_node_index_set[isimplex][1] = b_point_index + 1;
			simplex_decomposed_element_consisting_node_index_set[isimplex][2] = b_point_index;
			isimplex++;
		}

		return simplex_decomposed_element_consisting_node_index_set;
	}

	case FigureType::Quadrilateral:	{
		const size_t num_simplex = 2 * post_order * post_order;
		simplex_decomposed_element_consisting_node_index_set.resize(num_simplex);

		const size_t num_consisting_node = 3;


		size_t isimplex = 0;
		for (size_t j = 0; j < post_order; j++)
		for (size_t i = 0; i < post_order; i++)
		{
			//    b  --------- b+1
			//    |            |
			//    |            |  
			//    a  --------- a+1

			// recurrence relation
			// a(n+1) = a(n) + order + 1
			//=> a(n) = a(0) + (order + 1) * n
			const double a_point_index = static_cast<double>(i + 1 + (post_order + 1) * j);
			const double b_point_index = static_cast<double>(i + 1 + (post_order + 1) * (j + 1));

			simplex_decomposed_element_consisting_node_index_set[isimplex].resize(num_consisting_node);
			simplex_decomposed_element_consisting_node_index_set[isimplex][0] = a_point_index;
			simplex_decomposed_element_consisting_node_index_set[isimplex][1] = a_point_index + 1;
			simplex_decomposed_element_consisting_node_index_set[isimplex][2] = b_point_index;
			isimplex++;

			simplex_decomposed_element_consisting_node_index_set[isimplex].resize(num_consisting_node);
			simplex_decomposed_element_consisting_node_index_set[isimplex][0] = a_point_index + 1;
			simplex_decomposed_element_consisting_node_index_set[isimplex][1] = b_point_index + 1;
			simplex_decomposed_element_consisting_node_index_set[isimplex][2] = b_point_index;
			isimplex++;
		}

		return simplex_decomposed_element_consisting_node_index_set;
	}

	default:
		FATAL_TYPE_ERROR;
		return simplex_decomposed_element_consisting_node_index_set;
	}
}