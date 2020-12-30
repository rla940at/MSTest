#include "../INC/Monomial.h"

namespace Editor {
	std::string to_String(const Monomial& monomial){
		std::string str;
		if (monomial.is_Constant())		{
			str << "(1)";
			return str;
		}

		str << "(";
		for (size_t i = 0; i < monomial.num_variable(); ++i)		{
			if (monomial.exponent_set_[i] == 0)
				continue;
			else
				str << " (x" << i << ")^(" << monomial.exponent_set_[i] << ") ";
		}
		str << ")";

		return str;
	}
}


std::ostream& operator<<(std::ostream& ostream, const Monomial& monomial) {
	return ostream << Editor::to_String(monomial);
}


Monomial Monomial::operator*(const Monomial& other) const{
	const auto num_variable_this = this->num_variable();
	const auto num_variable_other = other.num_variable();

	std::vector<size_t> exponent_set(std::max(num_variable_this, num_variable_other));
	for (size_t i = 0; i < num_variable_this; ++i)
		exponent_set[i] += this->exponent_set_[i];
	for (size_t i = 0; i < num_variable_other; ++i)
		exponent_set[i] += other.exponent_set_[i];

	return std::move(exponent_set);
}

double Monomial::operator()(const MathVector& variable_vector) const{
	double result = 1.0;
	for (size_t i = 0; i < this->num_variable(); ++i)
		result *= std::pow(variable_vector[i], this->exponent_set_[i]);

	return result;
}

bool Monomial::operator<(const Monomial& other) const{
	const auto this_order = this->order();
	const auto other_order = other.order();

	if (this_order == other_order){
		auto min_num_variable = std::min(this->num_variable(), other.num_variable());
		for (size_t i = 0; i < min_num_variable; ++i){
			if (this->exponent_set_[i] == other.exponent_set_[i])
				continue;

			return this->exponent_set_[i] < other.exponent_set_[i];
		}
		return false;
	}

	return this_order < other_order;
}

size_t Monomial::exponent(size_t variable_index) const{
	if (this->num_variable() <= variable_index)
		return 0;
	else
		return this->exponent_set_[variable_index];
}

size_t Monomial::order(void) const{
	size_t order = 0;
	for (const auto exponent : this->exponent_set_)
		order += exponent;
	return order;
}

void Monomial::reduce_Order(const size_t variable_index){
	if (this->exponent(variable_index) == 0)
		*this = Monomial();
	else
		this->exponent_set_[variable_index]--;
}

bool Monomial::is_Constant(void) const{
	for (const auto exponent : this->exponent_set_)	{
		if (exponent != 0)
			return false;
	}
	return true;
}






















//
//
//void Monomial::showValue(void) const
//{
//	if (this->is_zero_)
//	{
//		std::cout << "(0)";
//		return;
//	}
//
//	std::cout << "(" << this->first_exponent_ << "," << this->second_exponent_ << "," << this->third_exponent_ << ")";
//}






//namespace
//{
//	std::vector<Monomial> getMonomialSet1D(const size_t num_monomial);
//	std::vector<Monomial> getMonomialSet2D(const size_t num_monomial);
//}
//
//std::vector<Monomial> MonomialSet::getMonomialSet(const size_t dimension, const size_t num_monomial)
//{
//	switch (dimension)
//	{
//	case 1:		return getMonomialSet1D(num_monomial);
//	case 2:		return getMonomialSet2D(num_monomial);
//	default:
//		ERROR(TOSTRING(dimension) + " dimensional monomial set is not supported.", __FILE__, __LINE__);
//		return std::vector<Monomial>();
//	}
//}












//
//
//
//namespace
//{
//	std::vector<Monomial> getMonomialSet1D(const size_t num_monomial)
//	{
//		std::vector<Monomial> monomials;
//		monomials.reserve(num_monomial);
//		for (size_t a = 0; a < num_monomial; ++a)
//			monomials.emplace_back(a);
//		return monomials;
//	}
//
//	std::vector<Monomial> getMonomialSet2D(const size_t num_monomial)
//	{
//		std::vector<Monomial> monomials;
//		monomials.reserve(num_monomial);
//
//		size_t index = 0;
//		size_t order = 0;
//		while (true)
//		{
//			const size_t a = order;
//			for (size_t b = 0; b <= a; ++b)
//			{
//				monomials.emplace_back(a - b, b);
//				++index;
//				if (index == num_monomial)
//					return monomials;
//			}
//			++order;
//		}
//	}
//}


//MonomialFunction1D::MonomialFunction1D(const size_t monomial_order) : monomial_order_(monomial_order)
//{
//	const size_t dimension = 1;
//	num_monomial_ = Math::calCombinationWithRepetition(1 + dimension, monomial_order);
//};
//
//std::vector<double> MonomialFunction1D::operator() (const Node& node) const
//{
//	const double r = node.getFirstCoord();
//
//	std::vector<double> monomials(num_monomial_);
//	monomials[0] = 1.0;
//	for (size_t i = 1; i < monomial_order_ + 1; i++)
//		monomials[i] = monomials[i - 1] * r;
//	return monomials;	//1 r r^2 ...
//}





//DivergenceMonomialFunction1D::DivergenceMonomialFunction1D(const size_t monomial_order) : monomial_order_(monomial_order)
//{
//	const size_t dimension = 1;
//	num_monomial_ = Math::calCombinationWithRepetition(1 + dimension, monomial_order);
//};
//std::vector<double> DivergenceMonomialFunction1D::operator()(const Node& node) const
//{
//	const double r = node.getFirstCoord();
//
//	std::vector<double> monomials(num_monomial_);
//	if (this->monomial_order_ == 0)
//		return monomials;
//	monomials[1] = 1.0;
//	for (size_t i = 2; i < this->monomial_order_ + 1; i++)
//		monomials[i] = static_cast<double>(i) / static_cast<double>(i - 1) * monomials[i - 1] * r;
//	return monomials;	//0 1 2r 3r^2 ...
//}
//
//
//MonomialFunction2D::MonomialFunction2D(const size_t monomial_order) : monomial_order_(monomial_order)
//{
//	const size_t dimension = 2;
//	num_monomial_ = Math::calCombinationWithRepetition(1 + dimension, monomial_order);
//};
//std::vector<double> MonomialFunction2D::operator() (const Node& node) const
//{
//	const double r = node.getFirstCoord();
//	const double s = node.getSecondCoord();
//
//	std::vector<double> monomials_2D(num_monomial_);
//	monomials_2D[0] = 1.0;
//	size_t index = 1;
//	for (size_t i = 1; i <= monomial_order_; i++) {
//		monomials_2D[index] = monomials_2D[index - i] * r;
//		index++;
//		for (size_t j = 0; j < i; j++) {
//			monomials_2D[index] = monomials_2D[index - i - 1] * s;
//			index++;
//		}
//	}
//	return monomials_2D;	//1 r s r^2 rs s^2 ...
//}
//
//
//DivergenceMonomialFunction2D::DivergenceMonomialFunction2D(const size_t monomial_order) : monomial_order_(monomial_order)
//{
//	const size_t dimension = 2;
//	num_monomial_ = Math::calCombinationWithRepetition(1 + dimension, monomial_order);
//};
//RowMajorMatrix DivergenceMonomialFunction2D::operator() (const Node& node) const
//{
//	const std::vector<double> dr_monomial_2D = this->DrMonomial2D(node);
//	const std::vector<double> ds_monomial_2D = this->DsMonomial2D(node);
//
//	const size_t dimension = 2;
//	const size_t num_row = dimension;
//	const size_t num_column = num_monomial_;
//
//	std::vector<double> value;
//	value.reserve(num_row * num_column);
//	VectorEditor::merge(value, dr_monomial_2D);
//	VectorEditor::merge(value, ds_monomial_2D);
//
//	return RowMajorMatrix(MatrixType::General, num_row, num_column, value).transpose();
//}
//
//
//std::vector<double> DivergenceMonomialFunction2D::DrMonomial2D(const Node& node) const
//{
//	const double r = node.getFirstCoord();
//	const double s = node.getSecondCoord();
//
//	std::vector<double> dr_monomials_2D(num_monomial_);
//	const DivergenceMonomialFunction1D divergence_monomial_function_1D(monomial_order_);
//	const std::vector<double> dr_monomials_1D = divergence_monomial_function_1D(Node(r));
//	
//	size_t index = 0;
//	for (size_t i = 0; i <= monomial_order_; i++)
//	{
//		dr_monomials_2D[index] = dr_monomials_1D[i];
//		index++;
//		for (size_t j = 0; j < i; j++) {
//			dr_monomials_2D[index] = dr_monomials_2D[index - i - 1] * s;
//			index++;
//		}
//	}
//	return dr_monomials_2D;	// 0 1 0 2r s 0 ...
//}
//std::vector<double> DivergenceMonomialFunction2D::DsMonomial2D(const Node& node) const
//{
//	const double r = node.getFirstCoord();
//	const double s = node.getSecondCoord();
//
//	std::vector<double> ds_monomials_2D(num_monomial_);
//	const DivergenceMonomialFunction1D divergence_monomial_function_1D(monomial_order_);
//	const std::vector<double> dr_monomials_1D = divergence_monomial_function_1D(Node(s));
//	
//	size_t index = 0;
//	for (size_t i = 0; i <= monomial_order_; i++)
//	{
//		for (size_t j = 0; j < i; j++) {
//			ds_monomials_2D[index] = ds_monomials_2D[index - i] * r;
//			index++;
//		}
//		ds_monomials_2D[index] = dr_monomials_1D[i];
//		index++;
//	}
//	return ds_monomials_2D;	// 0 0 1 0 r 2s ...
//}
//
//
//
//namespace
//{
//	std::vector<double> calTransformationMonomialVectorTriangle(const size_t order, const Node& point_in_figure);
//	std::vector<double> calTransformationMonomialVectorQuadrilateral(const size_t order, const Node& point_in_figure);
//
//	RowMajorMatrix calTransformationMonomialMatrixTriangle(const size_t element_order, const std::vector<Node>& points_in_reference_triangle);
//	RowMajorMatrix calTransformationMonomialMatrixQuadrilateral(const size_t element_order, const std::vector<Node>& points_in_reference_quadrilateral);
//
//	RowMajorMatrix calTransformationMonomialVectorJacobianMatrixTriangle(const size_t element_order, const Node& reference_quadrature_point);
//	RowMajorMatrix calTransformationMonomialVectorJacobianMatrixQuadrilateral(const size_t element_order, const Node& reference_quadrature_point);
//}
//
//std::vector<double> Monomials::calTransformationMonomialVector(const ReferenceFigureType figure_type, const size_t transformation_order, const Node& point_in_reference_figure)
//{
//	switch (figure_type)
//	{
//	case ReferenceFigureType::Triangle:			return calTransformationMonomialVectorTriangle(transformation_order, point_in_reference_figure);
//	case ReferenceFigureType::Quadrilateral:	return calTransformationMonomialVectorQuadrilateral(transformation_order, point_in_reference_figure);
//	default:
//		ERROR("This type of figure doesn't supprote transformation monomial matrix", __FILE__, __LINE__);
//		return std::vector<double>();
//	}
//
//}
//
//RowMajorMatrix Monomials::calTransformationMonomialMatrix(const ReferenceFigureType figure_type, const size_t transformation_order, const std::vector<Node>& points_in_reference_figure)
//{
//	switch (figure_type)
//	{
//	case ReferenceFigureType::Triangle:			return calTransformationMonomialMatrixTriangle(transformation_order, points_in_reference_figure);			
//	case ReferenceFigureType::Quadrilateral:	return calTransformationMonomialMatrixQuadrilateral(transformation_order, points_in_reference_figure);		
//	default:
//		ERROR("This reference figure type doesn't support transformation monomial matrix", __FILE__, __LINE__);
//		return RowMajorMatrix();
//	}
//}
//RowMajorMatrix Monomials::calTransformationMonomialVectorJacobianMatrix(const ReferenceFigureType figure_type, const size_t transformation_order, const Node& point_in_reference_figure)
//{
//	
//	switch (figure_type)
//	{
//	case ReferenceFigureType::Triangle:			return calTransformationMonomialVectorJacobianMatrixTriangle(transformation_order, point_in_reference_figure);
//	case ReferenceFigureType::Quadrilateral:	return calTransformationMonomialVectorJacobianMatrixQuadrilateral(transformation_order, point_in_reference_figure);	
//	default:
//		ERROR("This reference figure type doesn't support transformation monomial vector jacobian matrix", __FILE__, __LINE__);
//		return RowMajorMatrix();
//	}
//}
//
//
//
//namespace
//{
//	size_t numMonomialLine(const size_t order);
//	size_t numMonomialTriangle(const size_t order);
//	size_t numMonomialQuadrilateral(const size_t order);
//
//
//	std::vector<double> monomialDrLine(const size_t order, const double r);
//	std::vector<double> monomialDrTriangle(const size_t order, const double r, const double s);
//	std::vector<double> monomialDsTriangle(const size_t order, const double r, const double s);
//	std::vector<double> monomialDrQuadrilateral(const size_t order, const double r, const double s);
//	std::vector<double> monomialDsQuadrilateral(const size_t order, const double r, const double s);
//}
//
//
//
//
//namespace
//{
//	RowMajorMatrix calTransformationMonomialMatrixTriangle(const size_t transformation_order, const std::vector<Node>& points_in_reference_figure)
//	{
//		const size_t num_row = points_in_reference_figure.size();
//		const size_t num_column = numMonomialTriangle(transformation_order);
//		const MatrixType matrix_type = MatrixType::General;
//		std::vector<double> monomial_value;
//		monomial_value.reserve(num_row * num_column);
//
//		for (const Node& node : points_in_reference_figure)
//			VectorEditor::merge(monomial_value, calTransformationMonomialVectorTriangle(transformation_order, node));
//
//		RowMajorMatrix transposed_triangle_monomial_matrix(matrix_type, num_row, num_column, monomial_value);
//		transposed_triangle_monomial_matrix.transpose();
//
//		return transposed_triangle_monomial_matrix;
//	}
//
//	RowMajorMatrix calTransformationMonomialMatrixQuadrilateral(const size_t transformation_order, const std::vector<Node>& points_in_reference_figure)
//	{
//		const size_t num_row = points_in_reference_figure.size();
//		const size_t num_column = numMonomialQuadrilateral(transformation_order);
//		const MatrixType matrix_type = MatrixType::General;
//		std::vector<double> monomial_value;
//		monomial_value.reserve(num_row * num_column);
//
//		for (const Node& node : points_in_reference_figure)
//			VectorEditor::merge(monomial_value, calTransformationMonomialVectorQuadrilateral(transformation_order, node));
//
//		RowMajorMatrix transposed_quadrilateral_monomial_matrix(matrix_type, num_row, num_column, monomial_value);
//		transposed_quadrilateral_monomial_matrix.transpose();
//
//		return transposed_quadrilateral_monomial_matrix;
//	}
//}
//
//namespace
//{
//	RowMajorMatrix calTransformationMonomialVectorJacobianMatrixTriangle(const size_t transformation_order, const Node& reference_quadrature_point)
//	{
//		const size_t num_row = 2;
//		const size_t num_column = numMonomialTriangle(transformation_order);
//		const MatrixType matrix_type = MatrixType::General;
//
//		const double r1_coord = reference_quadrature_point.getFirstCoord();
//		const double r2_coord = reference_quadrature_point.getSecondCoord();
//
//		std::vector<double> monomial_dr_vector = monomialDrTriangle(transformation_order, r1_coord, r2_coord);
//		const std::vector<double> monomial_ds_vector = monomialDsTriangle(transformation_order, r1_coord, r2_coord);
//		VectorEditor::merge(monomial_dr_vector, monomial_ds_vector);
//
//		RowMajorMatrix transposed_triangle_monomial_jacobian_matrix(matrix_type, num_row, num_column, monomial_dr_vector);
//		transposed_triangle_monomial_jacobian_matrix.transpose();
//
//		return transposed_triangle_monomial_jacobian_matrix;
//	}
//
//	RowMajorMatrix calTransformationMonomialVectorJacobianMatrixQuadrilateral(const size_t transformation_order, const Node& reference_quadrature_point)
//	{
//		const size_t num_row = 2;
//		const size_t num_column = numMonomialQuadrilateral(transformation_order);
//		const MatrixType matrix_type = MatrixType::General;
//
//		const double r1_coord = reference_quadrature_point.getFirstCoord();
//		const double r2_coord = reference_quadrature_point.getSecondCoord();
//
//		std::vector<double> monomial_dr_vector = monomialDrQuadrilateral(transformation_order, r1_coord, r2_coord);
//		const std::vector<double> monomial_ds_vector =monomialDsQuadrilateral(transformation_order, r1_coord, r2_coord);
//		VectorEditor::merge(monomial_dr_vector, monomial_ds_vector);
//
//		RowMajorMatrix transposed_quadrilateral_monomial_jacobian_matrix(matrix_type, num_row, num_column, monomial_dr_vector);
//		transposed_quadrilateral_monomial_jacobian_matrix.transpose();
//		return transposed_quadrilateral_monomial_jacobian_matrix;
//	}
//}
//
//namespace
//{
//	size_t numMonomialLine(const size_t order)
//	{
//		return order + 1;
//	}
//	size_t numMonomialTriangle(const size_t order)
//	{
//		return (order + 1) * (order + 2) / 2;
//	}
//	size_t numMonomialQuadrilateral(const size_t order)
//	{
//		return (order + 1) * (order + 1);
//	}
//}
//
//namespace
//{
//	std::vector<double> monomialLine(const size_t order, const double r)
//	{
//		std::vector<double> results(numMonomialLine(order));
//		results[0] = 1.0;
//		for (size_t i = 1; i < order + 1; i++)
//			results[i] = results[i - 1] * r;
//		return results;	//1 r r^2 ...
//	}
//	std::vector<double> calTransformationMonomialVectorTriangle(const size_t transformation_order, const Node& point_in_figure)
//	{
//		const double r = point_in_figure.getFirstCoord();
//		const double s = point_in_figure.getSecondCoord();
//
//		std::vector<double> results(numMonomialTriangle(transformation_order));
//		results[0] = 1.0;
//		size_t index = 1;
//		for (size_t i = 1; i <= transformation_order; i++) {
//			results[index] = results[index - i] * r;
//			index++;
//			for (size_t j = 0; j < i; j++) {
//				results[index] = results[index - i - 1] * s;
//				index++;
//			}
//		}
//		return results;	//1 r s r^2 rs s^2 ...
//	}
//	std::vector<double> calTransformationMonomialVectorQuadrilateral(const size_t order, const Node& point_in_figure)
//	{
		//const double r = point_in_figure.getFirstCoord();
		//const double s = point_in_figure.getSecondCoord();

		//std::vector<double> results(numMonomialQuadrilateral(order));
		//results[0] = 1.0;
		//size_t index = 1;

		//for (size_t i = 1; i <= order; i++) {
		//	results[index] = results[index - 2 * i + 1] * r;
		//	index++;
		//	for (size_t j = 0; j < i; j++) {
		//		results[index] = results[index - 1] * s;
		//		index++;
		//	}
		//	for (size_t j = 0; j < i; j++) {
		//		results[index] = results[index - 2 * i - 1] * s;
		//		index++;
		//	}
		//}
		//return results;// 1 r rs s r^2 r^2s r^2s^2 rs^2 s^2 ...
//	}
//
//
//
//
//	std::vector<double> monomialDrLine(const size_t order, const double r)
//	{
//		std::vector<double> results(numMonomialLine(order));
//		results[0] = 0.0;
//		if (order == 0)
//			return results;
//		results[1] = 1.0;
//		for (size_t i = 2; i < order + 1; i++)
//			results[i] = double(i) / double(i - 1) * results[i - 1] * r;
//		return results;
//	}
//
//	std::vector<double> monomialDrTriangle(const size_t order, const double r, const double s)
//	{
//		std::vector<double> results(numMonomialTriangle(order));
//		const std::vector<double> monomials = monomialDrLine(order, r);
//		size_t index = 0;
//		for (size_t i = 0; i <= order; i++) {
//			results[index] = monomials[i];
//			index++;
//			for (size_t j = 0; j < i; j++) {
//				results[index] = results[index - i - 1] * s;
//				index++;
//			}
//		}
//		return results;
//	}
//	std::vector<double> monomialDsTriangle(const size_t order, const double r, const double s)
//	{
//		std::vector<double> results(numMonomialTriangle(order));
//		const std::vector<double> monomials = monomialDrLine(order, s);
//		size_t index = 0;
//		for (size_t i = 0; i <= order; i++) {
//			for (size_t j = 0; j < i; j++) {
//				results[index] = results[index - i] * r;
//				index++;
//			}
//			results[index] = monomials[i];
//			index++;
//		}
//		return results;
//	}
//	std::vector<double> monomialDrQuadrilateral(const size_t order, const double r, const double s)
//	{
//		std::vector<double> results(numMonomialQuadrilateral(order));
//		const std::vector<double> monomials = monomialDrLine(order, r);
//
//		size_t index = 0;
//		for (size_t i = 0; i <= order; i++) {
//			results[index] = monomials[i];
//			index++;
//			for (size_t j = 0; j < i; j++) {
//				results[index] = results[index - 1] * s;
//				index++;
//			}
//			for (size_t j = 0; j < i; j++) {
//				results[index] = results[index - 2 * i - 1] * s;
//				index++;
//			}
//		}
//		return results;
//	}
//	std::vector<double> monomialDsQuadrilateral(const size_t order, const double r, const double s)
//	{
//		std::vector<double> results(numMonomialQuadrilateral(order));
//		const std::vector<double> monomials = monomialDrLine(order, s);
//
//		results[0] = 0.0;
//		size_t index = 1;
//		for (size_t i = 1; i <= order; i++) {
//			for (size_t j = 0; j < i; j++) {
//				results[index] = results[index - 2 * i + 1] * r;
//				index++;
//			}
//			index = (i + 1) * (i + 1) - 1;
//			results[index] = monomials[i];
//			index--;
//			for (size_t j = 0; j < i; j++) {
//				results[index] = results[index + 1] * r;
//				index--;
//			}
//			index = (i + 1) * (i + 1);
//		}
//		return results;
//	}
//}






//size_t Monomial::numMonomialLine(const size_t order)
//{
//	return order + 1;
//}
//size_t Monomial::numMonomialTriangle(const size_t order)
//{
//	return (order + 1)*(order + 2) / 2;
//}
//size_t Monomial::numMonomialQuadrilateral(const size_t order)
//{
//	return (order + 1)*(order + 1);
//}
//size_t Monomial::numMonomialTets(const size_t order)
//{
//	return (order + 1)*(order + 2)*(order + 3) / 6;
//}
//size_t Monomial::numMonomialHexa(const size_t order)
//{
//	return (order + 1)*(order + 1)*(order + 1);
//}
//size_t Monomial::numMonomialPris(const size_t order)
//{
//	return (order + 1)*(order + 1)*(order + 2) / 2;
//}
//size_t Monomial::numMonomialPyra(const size_t order)
//{
//	return (order + 1)*(order + 2)*(2 * order + 3) / 6;
//}

//std::vector<double> Monomial::monomialLine(const size_t order, const double r)
//{
//	std::vector<double> results(numMonomialLine(order));
//	results[0] = 1.0;
//	for (size_t i = 1; i < order + 1; i++)
//		results[i] = results[i - 1] * r;
//	return results;	//1 r r^2 ...
//}
//std::vector<double> Monomial::monomialTriangle(const size_t order, const Node& point_in_figure)
//{
//	const double r = point_in_figure.getFirstCoord();
//	const double s = point_in_figure.getSecondCoord();
//
//	std::vector<double> results(numMonomialTriangle(order));
//	results[0] = 1.0;
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		results[index] = results[index - i] * r;
//		index++;
//		for (size_t j = 0; j < i; j++) {
//			results[index] = results[index - i - 1] * s;
//			index++;
//		}
//	}
//	return results;	//1 r s r^2 rs s^2 ...
//}
//std::vector<double> Monomial::monomialQuadrilateral(const size_t order, const Node& point_in_figure)
//{
//	const double r = point_in_figure.getFirstCoord();
//	const double s = point_in_figure.getSecondCoord();
//
//	std::vector<double> results(numMonomialQuadrilateral(order));
//	results[0] = 1.0;
//	size_t index = 1;
//
//	for (size_t i = 1; i <= order; i++) {
//		results[index] = results[index - 2 * i + 1] * r;
//		index++;
//		for (size_t j = 0; j < i; j++) {
//			results[index] = results[index - 1] * s;
//			index++;
//		}
//		for (size_t j = 0; j < i; j++) {
//			results[index] = results[index - 2 * i - 1] * s;
//			index++;
//		}
//	}
//	return results;// 1 r rs s r^2 r^2s r^2s^2 rs^2 s^2 ...
//}

//const std::vector<double> Monomial::monomialTets(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialTets(order));
//	results[0] = 1.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		results[index] = results[bef_index] * r;
//		index++;
//		def = index - bef_index;
//		bef_index = index - 1;
//		for (size_t j = 1; j <= i; j++) {
//			results[index] = results[index - def] * s;
//			index++; def++;
//			for (size_t k = 1; k <= j; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialHexa(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialHexa(order));
//	results[0] = 1.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		def = index - bef_index;
//		results[index] = results[bef_index] * r;
//		bef_index = index;
//		index++;
//		for (size_t j = 1; j <= i; j++) { // bottom
//			results[index] = results[index - 1] * s;
//			index++;
//		}
//		def += 2;
//		for (size_t j = 1; j <= i; j++) {
//			results[index] = results[index - def] * s;
//			index++;
//		}
//		def = 2 * i + 1;
//		for (size_t j = 1; j < i; j++) { // wall
//			for (size_t k = 0; k < def; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		// top
//		def = i*(3 * i + 1);
//		for (size_t j = 0; j < i; j++) {
//			for (size_t k = 0; k < 2 * j + 1; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		def = (i + 1)*(i + 1);
//		for (size_t j = 0; j < 2 * i + 1; j++) {
//			results[index] = results[index - def] * t;
//			index++;
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialPris(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPris(order));
//	results[0] = 1.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		results[index] = results[bef_index] * r;
//		index++;
//		def = index - bef_index;
//		bef_index = index - 1;
//		for (size_t j = 1; j <= i; j++) { // bottom
//			results[index] = results[index - def] * s;
//			index++;
//		}
//		def = i + 1;
//		for (size_t j = 1; j < i; j++) { // wall
//			for (size_t k = 0; k < def; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		// top
//		def = i*(i + 1) * 3 / 2;
//		for (size_t j = 1; j <= i; j++) {
//			for (size_t k = 1; k <= j; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		def = (i + 1)*(i + 2) / 2;
//		for (size_t j = 0; j < i + 1; j++) {
//			results[index] = results[index - def] * t;
//			index++;
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialPyra(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPyra(order));
//	results[0] = 1.0;
//
//	if (order == 0) return results;
//
//	double mix = 0.0;
//	if (std::abs(1.0 - t) > 1E-10)
//		mix = r*s / (1.0 - t);
//	const std::vector<double> mix_line = monomialLine(order, mix);
//	const std::vector<double> rst_tets = monomialTets(order, r, s, t);
//	const std::vector<double> rs_tris = monomialTriangle(order - 1, r, s);
//
//	std::vector<size_t> tets_index(order + 2);
//	std::vector<size_t> tris_index(order + 2);
//	for (size_t i = 0; i <= order + 1; i++) {
//		tets_index[i] = i*(i + 1)*(i + 2) / 6;
//		tris_index[i] = i*(i + 1) / 2;
//	}
//
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		for (size_t j = tets_index[i]; j < tets_index[i + 1]; j++)
//			results[index++] = rst_tets[j];
//		for (size_t j = 0; j < i; j++)
//			for (size_t k = tris_index[i - j - 1]; k < tris_index[i - j]; k++)
//				results[index++] = mix_line[j + 1] * rs_tris[k];
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialPyra2(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPyra(order));
//	results[0] = 1.0;
//
//	if (order == 0) return results;
//
//	const std::vector<double> rs_quad = monomialQuadrilateral(order, r, s);
//	const std::vector<double> t_line = monomialLine(order, t);
//
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		for (size_t j = 0; j <= i; j++) {
//			for (size_t k = j*j; k < (j + 1)*(j + 1); k++) {
//				results[index++] = rs_quad[k] * t_line[i - j];
//			}
//		}
//	}
//	return results;
//}
////const std::vector<double> Monomial::monomialDrLine(const size_t order, const double r)
////{
////	std::vector<double> results(numMonomialLine(order));
////	results[0] = 0.0;
////	if (order == 0)
////		return results;
////	results[1] = 1.0;
////	for (size_t i = 2; i < order + 1; i++)
////		results[i] = double(i) / double(i - 1)*results[i - 1] * r;
////	return results;
////}
////const std::vector<double> Monomial::monomialDrTriangle(const size_t order, const double r, const double s)
////{
////	std::vector<double> results(numMonomialTriangle(order));
////	const std::vector<double> monomials = monomialDrLine(order, r);
////	size_t index = 0;
////	for (size_t i = 0; i <= order; i++) {
////		results[index] = monomials[i];
////		index++;
////		for (size_t j = 0; j < i; j++) {
////			results[index] = results[index - i - 1] * s;
////			index++;
////		}
////	}
////	return results;
////}
////const std::vector<double> Monomial::monomialDsTriangle(const size_t order, const double r, const double s)
////{
////	std::vector<double> results(numMonomialTriangle(order));
////	const std::vector<double> monomials = monomialDrLine(order, s);
////	size_t index = 0;
////	for (size_t i = 0; i <= order; i++) {
////		for (size_t j = 0; j < i; j++) {
////			results[index] = results[index - i] * r;
////			index++;
////		}
////		results[index] = monomials[i];
////		index++;
////	}
////	return results;
////}
////const std::vector<double> Monomial::monomialDrQuadrilateral(const size_t order, const double r, const double s)
////{
////	std::vector<double> results(numMonomialQuadrilateral(order));
////	const std::vector<double> monomials = monomialDrLine(order, r);
////
////	size_t index = 0;
////	for (size_t i = 0; i <= order; i++) {
////		results[index] = monomials[i];
////		index++;
////		for (size_t j = 0; j < i; j++) {
////			results[index] = results[index - 1] * s;
////			index++;
////		}
////		for (size_t j = 0; j < i; j++) {
////			results[index] = results[index - 2 * i - 1] * s;
////			index++;
////		}
////	}
////	return results;
////}
////const std::vector<double> Monomial::monomialDsQuadrilateral(const size_t order, const double r, const double s)
////{
////	std::vector<double> results(numMonomialQuadrilateral(order));
////	const std::vector<double> monomials = monomialDrLine(order, s);
////
////	results[0] = 0.0;
////	size_t index = 1;
////	for (size_t i = 1; i <= order; i++) {
////		for (size_t j = 0; j < i; j++) {
////			results[index] = results[index - 2 * i + 1] * r;
////			index++;
////		}
////		index = (i + 1)*(i + 1) - 1;
////		results[index] = monomials[i];
////		index--;
////		for (size_t j = 0; j < i; j++) {
////			results[index] = results[index + 1] * r;
////			index--;
////		}
////		index = (i + 1)*(i + 1);
////	}
////	return results;
////}
//const std::vector<double> Monomial::monomialDrTets(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialTets(order));
//	const std::vector<double> monomials = monomialDrLine(order, r);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		results[index] = monomials[i];
//		index++;
//		def = index - bef_index;
//		bef_index = index - 1;
//		for (size_t j = 1; j <= i; j++) {
//			results[index] = results[index - def] * s;
//			index++; def++;
//			for (size_t k = 1; k <= j; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDsTets(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialTets(order));
//	const std::vector<double> monomials = monomialDrLine(order, s);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		def = index - bef_index;
//		bef_index = index;
//		for (size_t j = 0; j < i; j++) {
//			for (size_t k = 0; k <= j; k++) {
//				results[index] = results[index - def] * r;
//				index++;
//			}
//		}
//		results[index] = monomials[i];
//		index++;
//		def += (i + 1);
//		for (size_t k = 1; k <= i; k++) {
//			results[index] = results[index - def] * t;
//			index++;
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDtTets(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialTets(order));
//	const std::vector<double> monomials = monomialDrLine(order, t);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		def = index - bef_index;
//		bef_index = index;
//		for (size_t j = 0; j < i; j++) {
//			for (size_t k = 0; k <= j; k++) {
//				results[index] = results[index - def] * r;
//				index++;
//			}
//		}
//		def += i;
//		for (size_t k = 0; k < i; k++) {
//			results[index] = results[index - def] * s;
//			index++;
//		}
//		results[index] = monomials[i];
//		index++;
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDrHexa(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialHexa(order));
//	const std::vector<double> monomials = monomialDrLine(order, r);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		results[index] = monomials[i];
//		def = index - bef_index;
//		bef_index = index;
//		index++;
//		for (size_t j = 1; j <= i; j++) { // bottom
//			results[index] = results[index - 1] * s;
//			index++;
//		}
//		def += 2;
//		for (size_t j = 1; j < i; j++) {
//			results[index] = results[index - def] * s;
//			index++;
//		}
//		results[index++] = 0.0;
//		def = 2 * i + 1;
//		for (size_t j = 1; j < i; j++) { // wall
//			for (size_t k = 1; k < def; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//			results[index++] = 0.0;
//		}
//		// top
//		def = i*(3 * i + 1);
//		for (size_t j = 0; j < i; j++) {
//			for (size_t k = 1; k < 2 * j + 1; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//			results[index++] = 0.0;
//		}
//		def = (i + 1)*(i + 1);
//		for (size_t j = 1; j < 2 * i + 1; j++) {
//			results[index] = results[index - def] * t;
//			index++;
//		}
//		results[index++] = 0.0;
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDsHexa(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialHexa(order));
//	const std::vector<double> monomials = monomialDrLine(order, s);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		def = index - bef_index;
//		bef_index = index;
//		results[index++] = 0.0;
//		for (size_t j = 1; j < i; j++) {// bottom
//			results[index] = results[index - def] * r;
//			index++;
//		}
//		index += i;
//		results[index] = monomials[i];
//		for (size_t j = 1; j <= i; j++) {
//			index--;
//			results[index] = results[index + 1] * r;
//		}
//		index += (i + 1);
//		def = 2 * i + 1;
//		for (size_t j = 1; j < i; j++) { // wall
//			results[index++] = 0.0;
//			for (size_t k = 1; k < def; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		// top
//		def = i*(3 * i + 1);
//		for (size_t j = 0; j < i; j++) {
//			results[index++] = 0.0;
//			for (size_t k = 1; k < 2 * j + 1; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		def = (i + 1)*(i + 1);
//		results[index++] = 0.0;
//		for (size_t j = 1; j < 2 * i + 1; j++) {
//			results[index] = results[index - def] * t;
//			index++;
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDtHexa(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialHexa(order));
//	const std::vector<double> monomials = monomialDrLine(order, t);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		def = index - bef_index + 2;
//		bef_index = index;
//		for (size_t j = 0; j < 2 * i + 1; j++) { // bottom
//			results[index++] = 0.0;
//		}
//		for (size_t j = 1; j < i; j++) { // wall
//			if (j == i - 1) {
//				def -= (i - 1)*(i - 1);
//			}
//			results[index] = results[index - def] * r;
//			index++;
//			for (size_t k = 1; k <= i; k++) {
//				results[index] = results[index - 1] * s;
//				index++;
//			}
//			def += 2;
//			for (size_t k = 1; k <= i; k++) {
//				results[index] = results[index - def] * s;
//				index++;
//			}
//		}
//		// top
//		results[index++] = monomials[i];
//		def = 1;
//		for (size_t j = 1; j <= i; j++) {
//			results[index] = results[index - def] * r;
//			index++;
//			for (size_t k = 1; k <= j; k++) {
//				results[index] = results[index - 1] * s;
//				index++;
//			}
//			def += 2;
//			for (size_t k = 1; k <= j; k++) {
//				results[index] = results[index - def] * s;
//				index++;
//			}
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDrPris(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPris(order));
//	const std::vector<double> monomials = monomialDrLine(order, r);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		results[index] = monomials[i];
//		index++;
//		def = index - bef_index;
//		bef_index = index - 1;
//		for (size_t j = 1; j <= i; j++) { // bottom
//			results[index] = results[index - def] * s;
//			index++;
//		}
//		def = i + 1;
//		for (size_t j = 1; j < i; j++) { // wall
//			for (size_t k = 0; k < def; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		// top
//		def = i*(i + 1) * 3 / 2;
//		for (size_t j = 1; j <= i; j++) {
//			for (size_t k = 1; k <= j; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		def = (i + 1)*(i + 2) / 2;
//		for (size_t j = 0; j < i + 1; j++) {
//			results[index] = results[index - def] * t;
//			index++;
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDsPris(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPris(order));
//	const std::vector<double> monomials = monomialDrLine(order, s);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		def = index - bef_index;
//		bef_index = index;
//		results[index++] = 0.0;
//		for (size_t j = 1; j < i; j++) { // bottom
//			results[index] = results[index - def] * r;
//			index++;
//		}
//		results[index++] = monomials[i];
//		def = i + 1;
//		for (size_t j = 1; j < i; j++) { // wall
//			for (size_t k = 0; k < def; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		// top
//		def = i*(i + 1) * 3 / 2;
//		for (size_t j = 1; j <= i; j++) {
//			for (size_t k = 1; k <= j; k++) {
//				results[index] = results[index - def] * t;
//				index++;
//			}
//		}
//		def = (i + 1)*(i + 2) / 2;
//		for (size_t j = 0; j < i + 1; j++) {
//			results[index] = results[index - def] * t;
//			index++;
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDtPris(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPris(order));
//	const std::vector<double> monomials = monomialDrLine(order, t);
//
//	results[0] = 0.0;
//	size_t index = 1;
//	size_t bef_index = 0;
//	size_t def = 0;
//	for (size_t i = 1; i <= order; i++) {
//		def = index - bef_index;
//		bef_index = index;
//
//		for (size_t j = 1; j <= i + 1; j++) { // bottom
//			results[index++] = 0.0;
//		}
//
//		def++;
//		for (size_t j = 1; j < i; j++) { // wall
//			if (j == i - 1) {
//				def = i*(i + 1) - 1;
//			}
//			results[index] = results[index - def] * r;
//			index++;
//			def++;
//			for (size_t k = 1; k <= i; k++) {
//				results[index] = results[index - def] * s;
//				index++;
//			}
//		}
//
//		// top
//		results[index++] = monomials[i];
//		def = 1;
//		for (size_t j = 1; j <= i; j++) {
//			results[index] = results[index - def] * r;
//			index++;
//			def++;
//			for (size_t k = 1; k <= j; k++) {
//				results[index] = results[index - def] * s;
//				index++;
//			}
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDrPyra(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPyra(order));
//	results[0] = 0.0;
//
//	if (order == 0) return results;
//
//	double mix = 0.0;
//	double mix_dr = 0.0;
//	if (std::abs(1.0 - t) > 1E-10) {
//		mix = r*s / (1.0 - t);
//		mix_dr = s / (1.0 - t);
//	}
//
//	const std::vector<double> mix_line = monomialLine(order, mix);
//	const std::vector<double> rst_tets = monomialTets(order, r, s, t);
//	const std::vector<double> rs_tris = monomialTriangle(order - 1, r, s);
//	const std::vector<double> mix_line_dr = monomialDrLine(order, mix);
//	const std::vector<double> rst_tets_dr = monomialDrTets(order, r, s, t);
//	const std::vector<double> rs_tris_dr = monomialDrTriangle(order - 1, r, s);
//
//	std::vector<size_t> tets_index(order + 2);
//	std::vector<size_t> tris_index(order + 2);
//	for (size_t i = 0; i <= order + 1; i++) {
//		tets_index[i] = i*(i + 1)*(i + 2) / 6;
//		tris_index[i] = i*(i + 1) / 2;
//	}
//
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		for (size_t j = tets_index[i]; j < tets_index[i + 1]; j++)
//			results[index++] = rst_tets_dr[j];
//		for (size_t j = 0; j < i; j++)
//			for (size_t k = tris_index[i - j - 1]; k < tris_index[i - j]; k++)
//				results[index++] = mix_line_dr[j + 1] * mix_dr * rs_tris[k] + mix_line[j + 1] * rs_tris_dr[k];
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDsPyra(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPyra(order));
//	results[0] = 0.0;
//
//	if (order == 0) return results;
//
//	double mix = 0.0;
//	double mix_ds = 0.0;
//	if (std::abs(1.0 - t) > 1E-10) {
//		mix = r*s / (1.0 - t);
//		mix_ds = r / (1.0 - t);
//	}
//
//	const std::vector<double> mix_line = monomialLine(order, mix);
//	const std::vector<double> rst_tets = monomialTets(order, r, s, t);
//	const std::vector<double> rs_tris = monomialTriangle(order - 1, r, s);
//	const std::vector<double> mix_line_ds = monomialDrLine(order, mix);
//	const std::vector<double> rst_tets_ds = monomialDsTets(order, r, s, t);
//	const std::vector<double> rs_tris_ds = monomialDsTriangle(order - 1, r, s);
//
//	std::vector<size_t> tets_index(order + 2);
//	std::vector<size_t> tris_index(order + 2);
//	for (size_t i = 0; i <= order + 1; i++) {
//		tets_index[i] = i*(i + 1)*(i + 2) / 6;
//		tris_index[i] = i*(i + 1) / 2;
//	}
//
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		for (size_t j = tets_index[i]; j < tets_index[i + 1]; j++)
//			results[index++] = rst_tets_ds[j];
//		for (size_t j = 0; j < i; j++)
//			for (size_t k = tris_index[i - j - 1]; k < tris_index[i - j]; k++)
//				results[index++] = mix_line_ds[j + 1] * mix_ds * rs_tris[k] + mix_line[j + 1] * rs_tris_ds[k];
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDtPyra(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPyra(order));
//	results[0] = 0.0;
//
//	if (order == 0) return results;
//
//	double mix = 0.0;
//	double mix_dt = 0.0;
//	if (std::abs(1.0 - t) > 1E-10) {
//		mix = r*s / (1.0 - t);
//		mix_dt = mix / (1.0 - t);
//	}
//
//	const std::vector<double> mix_line = monomialLine(order, mix);
//	const std::vector<double> rst_tets = monomialTets(order, r, s, t);
//	const std::vector<double> rs_tris = monomialTriangle(order - 1, r, s);
//	const std::vector<double> mix_line_dt = monomialDrLine(order, mix);
//	const std::vector<double> rst_tets_dt = monomialDtTets(order, r, s, t);
//
//	std::vector<size_t> tets_index(order + 2);
//	std::vector<size_t> tris_index(order + 2);
//	for (size_t i = 0; i <= order + 1; i++) {
//		tets_index[i] = i*(i + 1)*(i + 2) / 6;
//		tris_index[i] = i*(i + 1) / 2;
//	}
//
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		for (size_t j = tets_index[i]; j < tets_index[i + 1]; j++)
//			results[index++] = rst_tets_dt[j];
//		for (size_t j = 0; j < i; j++)
//			for (size_t k = tris_index[i - j - 1]; k < tris_index[i - j]; k++)
//				results[index++] = mix_line_dt[j + 1] * mix_dt * rs_tris[k];
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDrPyra2(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPyra(order));
//	results[0] = 0.0;
//
//	if (order == 0) return results;
//
//	const std::vector<double> rs_quad = monomialDrQuadrilateral(order, r, s);
//	const std::vector<double> t_line = monomialLine(order, t);
//
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		for (size_t j = 0; j <= i; j++) {
//			for (size_t k = j*j; k < (j + 1)*(j + 1); k++) {
//				results[index++] = rs_quad[k] * t_line[i - j];
//			}
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDsPyra2(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPyra(order));
//	results[0] = 0.0;
//
//	if (order == 0) return results;
//
//	const std::vector<double> rs_quad = monomialDsQuadrilateral(order, r, s);
//	const std::vector<double> t_line = monomialLine(order, t);
//
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		for (size_t j = 0; j <= i; j++) {
//			for (size_t k = j*j; k < (j + 1)*(j + 1); k++) {
//				results[index++] = rs_quad[k] * t_line[i - j];
//			}
//		}
//	}
//	return results;
//}
//const std::vector<double> Monomial::monomialDtPyra2(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialPyra(order));
//	results[0] = 0.0;
//
//	if (order == 0) return results;
//
//	const std::vector<double> rs_quad = monomialQuadrilateral(order, r, s);
//	const std::vector<double> t_line = monomialDrLine(order, t);
//
//	size_t index = 1;
//	for (size_t i = 1; i <= order; i++) {
//		for (size_t j = 0; j <= i; j++) {
//			for (size_t k = j*j; k < (j + 1)*(j + 1); k++) {
//				results[index++] = rs_quad[k] * t_line[i - j];
//			}
//		}
//	}
//	return results;
//}
//
//const std::vector<double> Monomial::monomialTrisRpri(const size_t order, const double r, const double s)
//{
//	std::vector<double> results(numMonomialTriangle(order));
//	const std::vector<double> r_line = monomialLine(order, r);
//	const std::vector<double> s_line = monomialLine(order, s);
//	size_t index = 0;
//	for (size_t i = 0; i <= order; i++) {
//		for (size_t j = 0; j <= order - i; j++) {
//			results[index++] = s_line[i] * r_line[j];
//		}
//	}
//	return results;
//}
//
//const std::vector<double> Monomial::monomialQuadRpri(const size_t order, const double r, const double s)
//{
//	std::vector<double> results(numMonomialQuadrilateral(order));
//	const std::vector<double> r_line = monomialLine(order, r);
//	const std::vector<double> s_line = monomialLine(order, s);
//	size_t index = 0;
//	for (size_t i = 0; i <= order; i++) {
//		for (size_t j = 0; j <= order; j++) {
//			results[index++] = s_line[i] * r_line[j];
//		}
//	}
//	return results;
//}
//
//const std::vector<double> Monomial::monomialTetsRSpri(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialTets(order));
//	const std::vector<double> rs_tris = monomialTriangle(order, r, s);
//	const std::vector<double> t_line = monomialLine(order, t);
//	size_t index = 0;
//	for (size_t i = 0; i <= order; i++) {
//		for (size_t j = 0; j < numMonomialTriangle(order - i); j++) {
//			results[index++] = rs_tris[j] * t_line[i];
//		}
//	}
//	return results;
//}
//
//const std::vector<double> Monomial::monomialHexaRSpri(const size_t order, const double r, const double s, const double t)
//{
//	std::vector<double> results(numMonomialHexa(order));
//	const std::vector<double> rs_quad = monomialQuadrilateral(order, r, s);
//	const std::vector<double> t_line = monomialLine(order, t);
//	size_t index = 0;
//	for (size_t i = 0; i <= order; i++) {
//		for (size_t j = 0; j < numMonomialQuadrilateral(order); j++) {
//			results[index++] = rs_quad[j] * t_line[i];
//		}
//	}
//	return results;
//}
