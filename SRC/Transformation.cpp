#include "../INC/Transformation.h"

Transformation::Transformation
(
	const FigureType figure_type,
	const size_t transformation_order,
	const std::vector<MathVector>& original_position_vector_set,
	const std::vector<const MathVector*>& transformed_position_vector_set
)
{
	const auto num_original_position_vector = original_position_vector_set.size();
	const auto num_transformed_position_vector = transformed_position_vector_set.size();
	if (num_original_position_vector != num_transformed_position_vector)
		FATAL_ERROR(std::string() << "number of original position_vector(" << num_original_position_vector << ") and number of transformed position_vector(" << num_transformed_position_vector << ") is not 1:1 correspondence");

	//	C = X * INV(M)							
	//	C : transformation coefficient matrix	
	//	X : transformed position_vectors matrix			
	//	M : transformation monomial matrix	
	const auto transformation_monomial_set = this->calculate_Transformation_Monomial_Set(figure_type, transformation_order);
	const auto num_transformation_monomial = transformation_monomial_set.size();

	RowMajorMatrix M(MatrixType::Full, num_transformation_monomial, num_original_position_vector);
	for (size_t i = 0; i < num_transformation_monomial; ++i)
		for (size_t j = 0; j < num_original_position_vector; ++j)
			M.at(i, j) = transformation_monomial_set[i](original_position_vector_set[j]);

	constexpr size_t num_coord = 3;
	RowMajorMatrix X(MatrixType::Full, num_coord, num_transformed_position_vector);
	for (size_t i = 0; i < num_transformed_position_vector; ++i)
		X.change_Column(i, *transformed_position_vector_set[i]);

	const auto C = X * Math::inverse(M);
	const auto first_coord_trasnformation_coeffcient_set = C.row(0);
	const auto second_coord_trasnformation_coeffcient_set = C.row(1);
	const auto third_coord_trasnformation_coeffcient_set = C.row(2);

	this->transformation_function_.emplace_back(Polynomial(first_coord_trasnformation_coeffcient_set, transformation_monomial_set));
	this->transformation_function_.emplace_back(Polynomial(second_coord_trasnformation_coeffcient_set, transformation_monomial_set));
	this->transformation_function_.emplace_back(Polynomial(third_coord_trasnformation_coeffcient_set, transformation_monomial_set));

	constexpr size_t num_variable = 3;
	this->transformation_Jacobian_matrix_ = Math::jacobian(transformation_function_, num_variable);
}

Transformation::Transformation
(
	const FigureType figure_type,
	const size_t transformation_order,
	const std::vector<MathVector>& original_position_vector_set,
	const std::vector<MathVector>& transformed_position_vector_set
)
{
	const auto num_original_position_vector = original_position_vector_set.size();
	const auto num_transformed_position_vector = transformed_position_vector_set.size();
	if (num_original_position_vector != num_transformed_position_vector)
		FATAL_ERROR(std::string() << "number of original position_vector(" << num_original_position_vector << ") and number of transformed position_vector(" << num_transformed_position_vector << ") is not 1:1 correspondence");

	//	C = X * INV(M)							
	//	C : transformation coefficient matrix	
	//	X : transformed position_vectors matrix			
	//	M : transformation monomial matrix	
	const auto transformation_monomial_set = this->calculate_Transformation_Monomial_Set(figure_type, transformation_order);
	const auto num_transformation_monomial = transformation_monomial_set.size();

	RowMajorMatrix M(MatrixType::Full, num_transformation_monomial, num_original_position_vector);
	for (size_t i = 0; i < num_transformation_monomial; ++i)
	for (size_t j = 0; j < num_original_position_vector; ++j)
		M.at(i, j) = transformation_monomial_set[i](original_position_vector_set[j]);
	
	constexpr size_t num_coord = 3;
	RowMajorMatrix X(MatrixType::Full, num_coord, num_transformed_position_vector);
	for (size_t i = 0; i < num_transformed_position_vector; ++i)
		X.change_Column(i,transformed_position_vector_set[i]);

	const auto C = X * Math::inverse(M);
	const auto first_coord_trasnformation_coeffcient_set = C.row(0);
	const auto second_coord_trasnformation_coeffcient_set = C.row(1);
	const auto third_coord_trasnformation_coeffcient_set = C.row(2);

	this->transformation_function_.emplace_back(Polynomial(first_coord_trasnformation_coeffcient_set, transformation_monomial_set));
	this->transformation_function_.emplace_back(Polynomial(second_coord_trasnformation_coeffcient_set, transformation_monomial_set));
	this->transformation_function_.emplace_back(Polynomial(third_coord_trasnformation_coeffcient_set, transformation_monomial_set));

	constexpr size_t num_variable = 3;
	this->transformation_Jacobian_matrix_ = Math::jacobian(transformation_function_, num_variable);
}

MathVector Transformation::operator()(const MathVector& position_vector) const
{
	return this->transformation_function_(position_vector);
}

std::vector<MathVector> Transformation::operator()(const std::vector<MathVector>& position_vector_set) const
{
	std::vector<MathVector> transformed_position_vector_set;
	transformed_position_vector_set.reserve(position_vector_set.size());

	for (const auto& position_vector : position_vector_set)
		transformed_position_vector_set.emplace_back((*this)(position_vector));

	return transformed_position_vector_set;
}

double Transformation::calculate_Transformation_Scale(const FigureType figure_type, const MathVector& position_vector) const{
	switch (figure_type){
	case FigureType::Line:{
		const auto Jacobian_matrix = this->transformation_Jacobian_matrix_(position_vector);
		return Math::Frobenius_Norm(Jacobian_matrix);
	}
	case FigureType::Triangle:
	case FigureType::Quadrilateral:{
		const auto Jacobian_matrix = this->transformation_Jacobian_matrix_(position_vector);
		const auto cofactor_matrix = Math::cofactor_matrix(Jacobian_matrix);
		return Math::Frobenius_Norm(cofactor_matrix);
	}
	default:
		FATAL_TYPE_ERROR;
		return NULL;
	}
}

RowMajorMatrix Transformation::calculate_Transformation_Normal_Matrix(const FigureType figure_type, const MathVector& position_vector) const{
	switch (figure_type){
	case FigureType::Line: {
		constexpr size_t restricted_dimension = 2;
		const auto Jacobian_matrix = this->transformation_Jacobian_matrix_(position_vector);
		const auto Jacobian_sub_matrix = Jacobian_matrix.part(restricted_dimension);
		return Math::cofactor_matrix(Jacobian_sub_matrix);
	}
	case FigureType::Triangle:
	case FigureType::Quadrilateral: {
		const auto Jacobian_matrix = this->transformation_Jacobian_matrix_(position_vector);
		return Math::cofactor_matrix(Jacobian_matrix);
	}
	default:
		FATAL_TYPE_ERROR;
		return RowMajorMatrix();
	}
}

std::vector<Monomial> Transformation::calculate_Transformation_Monomial_Set(const FigureType figure_type, const size_t transformation_order) const{
	std::vector<Monomial> transformation_monomials;

	switch (figure_type){
	case FigureType::Line:{
		const size_t num_transformation_monomial = transformation_order + 1;
		transformation_monomials.reserve(num_transformation_monomial);

		for (size_t a = 0; a <= transformation_order; ++a)
			transformation_monomials.emplace_back(a);

		return transformation_monomials;	// 1 r r^2 ...
	}
	case FigureType::Triangle:{
		const size_t num_transformation_monomial = static_cast<size_t>((transformation_order + 2) * (transformation_order + 1) * 0.5);
		transformation_monomials.reserve(num_transformation_monomial);

		for (size_t a = 0; a <= transformation_order; ++a)
		for (size_t b = 0; b <= a; ++b)
			transformation_monomials.emplace_back(a - b, b);

		return transformation_monomials;	// 1 r s r^2 rs s^2 ...
	}
	case FigureType::Quadrilateral:{
		const size_t num_transformation_monomial = static_cast<size_t>((transformation_order + 1) * (transformation_order + 1));
		transformation_monomials.reserve(num_transformation_monomial);

		for (size_t a = 0; a <= transformation_order; ++a)
		{
			for (size_t b = 0; b <= a; ++b)
				transformation_monomials.emplace_back(a, b);

			if (a == 0)	continue;

			size_t index = a;
			while (true){
				transformation_monomials.emplace_back(--index, a);
				if (index == 0)		break;
			}
		}

		return transformation_monomials;	// 1 r rs s r^2 r^2s r^2s^2 rs^2 s^2...
	}
	default:
		FATAL_TYPE_ERROR;
		return transformation_monomials;
	}
}