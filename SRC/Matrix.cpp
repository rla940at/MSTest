
#include "../INC/Matrix.h"

namespace Math{
	RowMajorMatrix& abs(RowMajorMatrix& A) {
		Math::abs(A.value_);
		return A;
	}

	RowMajorMatrix cofactor_matrix(const RowMajorMatrix& A) { //	C_i,j = (-1)^(i+j) * (M_i,j)
		const auto [row, column] = A.size();
		RowMajorMatrix cofactor_matrix(MatrixType::Full, row, column);
		for (size_t i = 0; i < row; ++i)
			for (size_t j = 0; j < column; ++j)
				cofactor_matrix.at(i, j) = pow(-1.0, i + j) * Math::minor(A, i, j);

		return cofactor_matrix;
	}

	double determinant(RowMajorMatrix&& A) {
		if (!A.is_Square_Matrix())
			FATAL_TYPE_ERROR;

		A.update_Type();
		if (A.type_ == MatrixType::Full)
			return A.determinant_Full_Matrix();
		else
			return A.determinant_Triangle_Matrix();
	}	

	double determinant(const RowMajorMatrix& A) {
		return Math::determinant(RowMajorMatrix(A));
	}

	double Frobenius_Norm(const RowMajorMatrix& A){
		return Math::L2_norm(A.value_);
	}

	RowMajorMatrix& inverse(RowMajorMatrix& A) {
		if (A.type_ == MatrixType::Full)
			return A.inverse_Full_Matrix();
		else if (A.type_ == MatrixType::FullLowerTriangle || A.type_ == MatrixType::FullUpperTriangle)
			return A.inverse_Full_Triangle_Matrix();
		else 
			return A.inverse_Packed_Triangle_Matrix();
	}

	RowMajorMatrix inverse(const RowMajorMatrix& A) {
		RowMajorMatrix result = A;
		return Math::inverse(result);
	}

	double minor(const RowMajorMatrix& A, const size_t deleting_row_index, const size_t deleting_column_index){ //	M_i,j := Minor => determinant of eliminate i row j column
		auto minor_matrix = A.part(deleting_row_index, deleting_column_index);
		return Math::determinant(std::move(minor_matrix));
	}

	void Moore_Penrose_Inverse(RowMajorMatrix& A){
		//https://software.intel.com/content/www/us/en/develop/articles/implement-pseudoinverse-of-a-matrix-by-intel-mkl.html
		//https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/lapacke_dgesvd_row.c.htm
		//https://bskyvision.com/256

		// {U , S , tr(V)}
		auto decomposed_matrix_set = Math::singular_value_decomposition(A);

		auto& trU = Math::transpose(decomposed_matrix_set[0]);
		auto& S = decomposed_matrix_set[1];

		const auto [S_row, S_column] = S.size();
		const auto num_diagonal = std::min(S_row, S_column);
		for (size_t i = 0; i < num_diagonal; ++i){
			const auto singular_value = S.at(i, i);
			if (singular_value > 1.0E-9)
				S.at(i, i) = 1.0 / singular_value;
		}

		auto& PinvS = Math::transpose(S);
		auto& V = Math::transpose(decomposed_matrix_set[2]);

		//	Pinv(A) = V * Pinv(S) * tr(U)
		A = V * PinvS * trU;
	}
		
	std::vector<int> PLU_decomposition(RowMajorMatrix& A) { // A=PLU, A is overwriiten by L + U. The unit diagonal element of L are not stored		
		if (A.type_ != MatrixType::Full)
			FATAL_TYPE_ERROR;

		if (A.is_transposed())
			A.transpose_Value();

		const int matrix_layout = LAPACK_ROW_MAJOR;
		const lapack_int m = static_cast<int>(A.row_);
		const lapack_int n = static_cast<int>(A.column_);
		const lapack_int lda = n;
		std::vector<int> ipiv(std::min(m, n));
		lapack_int info = LAPACKE_dgetrf(matrix_layout, m, n, A.data(), lda, ipiv.data());

		if (info < 0)
			FATAL_ERROR("Fail to decompose");

		return ipiv;
	}

	std::array<RowMajorMatrix, 3> singular_value_decomposition(RowMajorMatrix& A) {
		//	singular value decomposition
		//	A = U * S * tr(V)	
		//	A := m * n matrix
		//	U := m * m orthogonal matrix
		//	S := m * n matrix which is zero excep for its min(m,n) diagonal elements
		//	tr(V) := transpoed n * n orthogonal matrix which is zero except for its min(m,n) diagonal elements
		//	return {U,S,trV}

		if (!A.is_Full_Storage_Type())
			A.to_Full_Matrix();

		if (A.is_transposed())
			A.transpose_Value();

		const auto [row, column] = A.size();

		const int matrix_layout = LAPACK_ROW_MAJOR;
		const char jobu = 'A';
		const char jobvt = 'A';
		const lapack_int m = static_cast<int>(row);
		const lapack_int n = static_cast<int>(column);
		const lapack_int lda = static_cast<int>(A.Leading_Dimension());
		const lapack_int ldu = m;
		const lapack_int ldvt = n;

		const auto num_singular_value = std::min(m, n);
		const auto maximum_num_unconverged_diagonal = num_singular_value - 1;

		RowMajorMatrix U(MatrixType::Full, m, m);
		RowMajorMatrix S(MatrixType::Full, m, n);
		RowMajorMatrix trV(MatrixType::Full, n, n);
		std::vector<double> singular_value(num_singular_value);
		std::vector<double> superb(maximum_num_unconverged_diagonal);

		const auto info = LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, A.data(), lda, singular_value.data(), U.data(), ldu, trV.data(), ldvt, superb.data());

		for (size_t i = 0; i < num_singular_value; ++i)
			S.at(i, i) = singular_value[i];

		if (info < 0)
			FATAL_ERROR(info << "th input of LAPACKE_dgesvd had an illegal value");

		return { U,S,trV };
	}

	RowMajorMatrix& transpose(RowMajorMatrix& A) {
		std::swap(A.row_, A.column_);

		if (A.is_transposed())
			A.transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
		else
			A.transpose_type_ = CBLAS_TRANSPOSE::CblasTrans;

		switch (A.type_) {
		case MatrixType::FullLowerTriangle:
			A.type_ = MatrixType::FullUpperTriangle;
			return A;
		case MatrixType::PackedLowerTriangle:
			A.type_ = MatrixType::PackedUpperTriangle;
			return A;
		case MatrixType::FullUpperTriangle:
			A.type_ = MatrixType::FullLowerTriangle;
			return A;
		case MatrixType::PackedUpperTriangle:
			A.type_ = MatrixType::PackedLowerTriangle;
			return A;
		default:
			return A;
		}
	}

	RowMajorMatrix transpose(const RowMajorMatrix& A) {
		RowMajorMatrix result = A;
		return Math::transpose(result);
	}

	RowMajorMatrix transpose(RowMajorMatrix&& A) {
		RowMajorMatrix result = std::move(A);
		return Math::transpose(result);
	}
}


namespace Editor {
	std::string to_String(const RowMajorMatrix& A) {
		std::string str;
		switch (A.type_) {
		case MatrixType::Full: {
			str << "[Full Matrix]\n";
			for (size_t i = 0; i < A.row_; ++i)	{
				for (size_t j = 0; j < A.column_; ++j)
					str << A.at(i, j) << "\t";
				str << "\n\n";
			}
			return str;
		}
		case MatrixType::FullLowerTriangle:
		case MatrixType::PackedLowerTriangle: {
			str << "[Lower Triangle Matrix]\n";
			for (size_t i = 0; i < A.row_; ++i) {
				for (size_t j = 0; j <= i; ++j)
					str << A.at(i, j) << "\t";
				str << "\n\n";
			}
			return str;
		}
		case MatrixType::FullUpperTriangle:
		case MatrixType::PackedUpperTriangle: {
			str << "[Upper Triangle Matrix]\n";
			for (size_t i = 0; i < A.row_; ++i) {
				for (size_t j = 0; j < i; ++j)
					str << "\t";
				for (size_t j = i; j < A.column_; ++j)
					str << A.at(i, j) << "\t";
				str << "\n\n";
			}
			return str;
		}
		default:
			FATAL_TYPE_ERROR;
			return str;
		}
	}
}


std::ostream& operator<<(std::ostream& os, const RowMajorMatrix& A) {
	return os << Editor::to_String(A);
}

RowMajorMatrix operator*(const double scalar, const RowMajorMatrix& A) {
	return A * scalar;
}


RowMajorMatrix::RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column)
	:type_(matrix_type), row_(num_row), column_(num_column) {
	size_t num_value = 0;
	if (this->is_Full_Storage_Type())
		num_value = this->row_ * this->column_;
	else if (this->is_Square_Matrix())
		num_value = static_cast<size_t>((this->row_) * (this->row_ + 1) * 0.5);
	else
		FATAL_TYPE_ERROR;

	this->value_.resize(num_value);
}

RowMajorMatrix::RowMajorMatrix(const MatrixType matrix_type, const size_t matrix_order)
	: type_(matrix_type), row_(matrix_order), column_(matrix_order) {
	size_t num_value = 0;
	if (this->is_Full_Storage_Type())
		num_value = this->row_ * this->column_;
	else 
		num_value = static_cast<size_t>((this->row_) * (this->row_ + 1) * 0.5);

	this->value_.resize(num_value);

	for (size_t i = 0; i < matrix_order; ++i)
		this->at(i, i) = 1;
}

RowMajorMatrix::RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column, const MathVector& value)
	: type_(matrix_type), row_(num_row), column_(num_column), value_(value) {
	this->inspect_Size();
	this->inspect_Value();
}

RowMajorMatrix::RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column, MathVector&& value)
	: type_(matrix_type), row_(num_row), column_(num_column), value_(std::move(value)) {
	this->inspect_Size();
	this->inspect_Value();
}

RowMajorMatrix& RowMajorMatrix::operator+=(const RowMajorMatrix& B) {
	RowMajorMatrix& A = *this;

	if (A.size() != B.size())
		FATAL_SIZE_ERROR;

	if (!A.is_Full_Storage_Type())
		A.to_Full_Matrix();

	if (A.is_transposed())
		A.transpose_Value();

	if (!B.is_Full_Storage_Type()) {
		auto full_B = B.to_Full_Matrix();

		if (full_B.is_transposed())
			full_B.transpose_Value();
		
		A.value_ += full_B.value_;

		A.update_Type();
		return A;
	}
	else {
		if (B.is_transposed())
			A.value_ += B.transpose_Value().value_;
		else
			A.value_ += B.value_;

		A.update_Type();
		return A;
	}
}

RowMajorMatrix RowMajorMatrix::operator+(const RowMajorMatrix& B) const{
	RowMajorMatrix A = *this;
	return A += B;
}

RowMajorMatrix& RowMajorMatrix::operator-=(const RowMajorMatrix& B){
	RowMajorMatrix& A = *this;

	if (A.size() != B.size())
		FATAL_SIZE_ERROR;

	if (!A.is_Full_Storage_Type())
		A.to_Full_Matrix();

	if (A.is_transposed())
		A.transpose_Value();

	if (!B.is_Full_Storage_Type()) {
		auto full_B = B.to_Full_Matrix();

		if (full_B.is_transposed())
			full_B.transpose_Value();

		A.value_ -= full_B.value_;

		A.update_Type();
		return A;
	}
	else {
		if (B.is_transposed())
			A.value_ -= B.transpose_Value().value_;
		else
			A.value_ -= B.value_;

		A.update_Type();
		return A;
	}
}

RowMajorMatrix RowMajorMatrix::operator-(const RowMajorMatrix& B) const{
	RowMajorMatrix A = *this;
	return A -= B;
}

RowMajorMatrix& RowMajorMatrix::operator*=(const double scalar){
	this->value_ *= scalar;
	return *this;
}

RowMajorMatrix RowMajorMatrix::operator*(const double scalar) const {
	RowMajorMatrix A = *this;
	return A *= scalar;
}

MathVector RowMajorMatrix::operator*(const MathVector& x) const {
	switch (this->type_) {
	case MatrixType::Full:
		return this->multiply_Full_Matrix_Vector(x);
	case MatrixType::FullLowerTriangle:
	case MatrixType::FullUpperTriangle:
		return this->multiply_Full_Triangle_Matrix_Vector(x);
	case MatrixType::PackedLowerTriangle:
	case MatrixType::PackedUpperTriangle:
		return this->multiply_Packed_Triangle_Matrix_Vector(x);
	default:
		FATAL_TYPE_ERROR;
		return MathVector();
	}
}

RowMajorMatrix& RowMajorMatrix::operator*=(const RowMajorMatrix& B) {
	if (this->is_Full_Storage_Type()) {
		switch (B.type_) {
		case MatrixType::Full:
			return this->multiply_Full_Matrix_Full_Matrix(B);
		case MatrixType::FullLowerTriangle:
		case MatrixType::FullUpperTriangle:
			return this->multiply_Full_Matrix_Full_Triangle_Matrix(B);
		case MatrixType::PackedLowerTriangle:
		case MatrixType::PackedUpperTriangle:
			return this->multiply_Full_Matrix_Full_Triangle_Matrix(B.to_Full_Matrix());
		default:
			FATAL_TYPE_ERROR;
			return *this;
		}
	}
	else {
		this->to_Full_Matrix();
		switch (B.type_) {
		case MatrixType::Full:
			return this->multiply_Full_Matrix_Full_Matrix(B);
		case MatrixType::FullLowerTriangle:
		case MatrixType::FullUpperTriangle:
			return this->multiply_Full_Matrix_Full_Triangle_Matrix(B);
		case MatrixType::PackedLowerTriangle:
		case MatrixType::PackedUpperTriangle:
			return this->multiply_Full_Matrix_Full_Triangle_Matrix(B.to_Full_Matrix());
		default:
			FATAL_TYPE_ERROR;
			return *this;
		}
	}
}

RowMajorMatrix RowMajorMatrix::operator*(const RowMajorMatrix& B) const{
	RowMajorMatrix A = *this;
	return A *= B;
}

bool RowMajorMatrix::operator==(const RowMajorMatrix& B) const {
	if (this->transpose_type_ != B.transpose_type_ ||
		this->type_ != B.type_ ||
		this->row_ != B.row_ ||
		this->column_ != B.column_ ||
		this->value_ != B.value_)
		return false;
	else
		return true;
}

double& RowMajorMatrix::at(const size_t irow, const size_t jcolumn){
	this->inspect_Range(irow, jcolumn);
	if (this->is_transposed()) {
		switch (this->type_) {
		case MatrixType::Full:
		case MatrixType::FullLowerTriangle:
		case MatrixType::FullUpperTriangle:			
			return this->value_[jcolumn * this->row_ + irow];
		case MatrixType::PackedLowerTriangle: {
			if (this->is_Upper_Part(irow, jcolumn))
				FATAL_ERROR("In packed lower triangle matrix, value reference of upper part doesn't supported");

			size_t position = 0;
			for (size_t i = 0; i <= jcolumn; ++i)
				position += this->row_ - i;

			position += irow - jcolumn;

			return this->value_[position];
		}
		case MatrixType::PackedUpperTriangle: {
			if (this->is_Lower_Part(irow, jcolumn))
				FATAL_ERROR("In packed upper triangle matrix, value reference of lower part doesn't supported");

			size_t position = 0;
			for (size_t i = 0; i <= jcolumn; ++i)
				position += i;

			position += irow;

			return this->value_[position];
		}
		default:
			FATAL_TYPE_ERROR;
			return *this->value_.data();
		}
	}
	else {
		switch (this->type_) {
		case MatrixType::Full:
		case MatrixType::FullLowerTriangle:
		case MatrixType::FullUpperTriangle:
			return this->value_[irow * this->column_ + jcolumn];
		case MatrixType::PackedLowerTriangle: {
			if (this->is_Upper_Part(irow, jcolumn))
				FATAL_ERROR("In packed lower triangle matrix, value reference of upper part deosn't supported");

			size_t position = 0;
			position += static_cast<size_t>(irow * (irow + 1) * 0.5);

			position += jcolumn;

			return this->value_[position];
		}
		case MatrixType::PackedUpperTriangle: {
			if (this->is_Lower_Part(irow, jcolumn))
				FATAL_ERROR("In packed upper triangle matrix, value reference of lower part deosn't supported");

			size_t position = 0;
			for (size_t i = 0; i < irow; ++i)
				position += this->row_ - i;

			position += jcolumn - irow;

			return this->value_[position];
		}
		default:
			FATAL_TYPE_ERROR;
			return *this->value_.data();
		}
	}
}

double RowMajorMatrix::at(const size_t irow, const size_t jcolumn) const {
	this->inspect_Range(irow, jcolumn);
	if (this->is_transposed()) {
		switch (this->type_) {
		case MatrixType::Full:
		case MatrixType::FullLowerTriangle:
		case MatrixType::FullUpperTriangle:
			return this->value_[jcolumn * this->row_ + irow];
		case MatrixType::PackedLowerTriangle: {
			if (this->is_Upper_Part(irow, jcolumn))
				return NULL;

			size_t position = 0;
			for (size_t i = 0; i <= jcolumn; ++i)
				position += this->row_ - i;

			position += irow - jcolumn;

			return this->value_[position];
		}
		case MatrixType::PackedUpperTriangle: {
			if (this->is_Lower_Part(irow, jcolumn))
				return NULL;

			size_t position = 0;
			for (size_t i = 0; i <= jcolumn; ++i)
				position += i;

			position += irow;

			return this->value_[position];
		}
		default:
			FATAL_TYPE_ERROR;
			return *this->value_.data();
		}
	}
	else {
		switch (this->type_) {
		case MatrixType::Full:
		case MatrixType::FullLowerTriangle:
		case MatrixType::FullUpperTriangle:
			return this->value_[irow * this->column_ + jcolumn];
		case MatrixType::PackedLowerTriangle: {
			if (this->is_Upper_Part(irow, jcolumn))
				return NULL;

			size_t position = 0;
			position += static_cast<size_t>(irow * (irow + 1) * 0.5);

			position += jcolumn;

			return this->value_[position];
		}
		case MatrixType::PackedUpperTriangle: {
			if (this->is_Lower_Part(irow, jcolumn))
				return NULL;

			size_t position = 0;
			for (size_t i = 0; i <= irow; ++i)
				position += this->row_ - i;

			position += jcolumn - irow;

			return this->value_[position];
		}
		default:
			FATAL_TYPE_ERROR;
			return *this->value_.data();
		}
	}
}

void RowMajorMatrix::change_Column(const size_t column_index, const MathVector& vector) {
	if (this->row_ != vector.size())
		FATAL_SIZE_ERROR;

	for (size_t i = 0; i < this->row_; ++i)
		this->at(i, column_index) = vector[i];
}

void RowMajorMatrix::change_Row(const size_t row_index, const MathVector& vector) {
	if (this->column_ != vector.size())
		FATAL_SIZE_ERROR;

	for (size_t i = 0; i < this->column_; ++i)
		this->at(row_index, i) = vector[i];
}

RowMajorMatrix RowMajorMatrix::part(const size_t start_row, const size_t end_row, const size_t start_column, const size_t end_column) const {
	const auto maximum_row_index = end_row - 1;
	const auto maximum_column_index = end_column - 1;
	this->inspect_Range(maximum_row_index, maximum_column_index);

	const auto num_part_row = end_row - start_row;
	const auto num_part_column = end_column - start_column;
	RowMajorMatrix part_matrix(MatrixType::Full, num_part_row, num_part_column);
	for (size_t i = 0; i < num_part_row; ++i)
		for (size_t j = 0; j < num_part_column; ++j)
			part_matrix.at(i, j) = this->at(start_row + i, start_column + j);

	return part_matrix;
}

RowMajorMatrix RowMajorMatrix::part(const size_t part_matrix_order) const {
	const size_t start_row = 0;
	const size_t end_row = part_matrix_order - 1;
	const size_t start_column = 0;
	const size_t end_column = part_matrix_order - 1;

	return this->part(start_row, end_row, start_column, end_column);
}

RowMajorMatrix RowMajorMatrix::part(const size_t deleting_row_index, const size_t deleting_column_index) const {
	this->inspect_Range(deleting_row_index, deleting_column_index);

	const size_t num_part_row = this->row_ - 1;
	const size_t num_part_column = this->column_ - 1;
	MathVector value(num_part_row * num_part_column);

	size_t index = 0;
	for (size_t i = 0; i < this->row_; ++i) {
		if (i == deleting_row_index)
			continue;

		for (size_t j = 0; j < this->column_; ++j) {
			if (j == deleting_column_index)
				continue;

			value[index++] = this->at(i, j);
		}
	}

	RowMajorMatrix part_matrix(MatrixType::Full, num_part_row, num_part_column, std::move(value));
	return part_matrix;
}

MathVector RowMajorMatrix::row(const size_t row_index) const {
	this->inspect_Range(row_index, this->column_ - 1);

	MathVector row_vector(this->column_);
	for (size_t j = 0; j < this->column_; ++j)
		row_vector[j] = this->at(row_index, j);

	return row_vector;
}

RowMajorMatrix RowMajorMatrix::row(const size_t start_row, const size_t end_row) const{
	if (start_row > end_row)
		FATAL_SIZE_ERROR;
	
	const size_t start_column = 0;
	const size_t end_column = this->column_;
	return this->part(start_row, end_row, start_column, end_column);
}

MathVector RowMajorMatrix::column(const size_t column_index) const{
	this->inspect_Range(this->row_ - 1, column_index);

	MathVector column_vector(this->row_);
	for (size_t i = 0; i < this->row_; ++i)
		column_vector[i] = this->at(i, column_index);

	return column_vector;
}

RowMajorMatrix RowMajorMatrix::column(const size_t start_column, const size_t end_column) const{
	if (start_column > end_column)
		FATAL_SIZE_ERROR;

	const size_t start_row = 0;
	const size_t end_row = this->row_;
	return this->part(start_row, end_row, start_column, end_column);
}

RowMajorMatrix& RowMajorMatrix::scale_Row(const size_t start_row, const size_t end_row, const double scale_factor){
	if (start_row > end_row)
		FATAL_SIZE_ERROR;

	this->inspect_Range(end_row - 1, this->column_ - 1);

	for (size_t i = start_row; i < end_row; ++i)
		for (size_t j = 0; j < this->column_; ++j)
			this->at(i, j) *= scale_factor;

	return *this;
}

RowMajorMatrix& RowMajorMatrix::scale_Column(const size_t start_column, const size_t end_column, const double scale_factor){
	if (start_column > end_column)
		FATAL_SIZE_ERROR;

	this->inspect_Range(this->row_ - 1, end_column - 1);

	for (size_t i = 0; i < this->row_; ++i)
		for (size_t j = start_column; j < end_column; ++j)
			this->at(i, j) *= scale_factor;

	return *this;
}

CBLAS_DIAG RowMajorMatrix::CBALS_Diagonal_Type(void) const {
	if (this->is_Unit_Diagonal())
		return CBLAS_DIAG::CblasUnit;
	else
		return CBLAS_DIAG::CblasNonUnit;
}

CBLAS_UPLO RowMajorMatrix::CBLAS_Triangle_Type(void) const {
	if (this->is_transposed()) {
		// MKL routine needs original triangle matrix type to calculate tranposed matrix multiplication
		if(this->is_Upper_Triangle_Type())
			return CBLAS_UPLO::CblasLower;
		else
			return CBLAS_UPLO::CblasUpper;
	}
	else {
		if (this->is_Upper_Triangle_Type())
			return CBLAS_UPLO::CblasUpper;
		else
			return CBLAS_UPLO::CblasLower;
	}
}

double RowMajorMatrix::determinant_Full_Matrix(void) {//https://en.wikipedia.org/wiki/LU_decomposition
	const auto pivot_indices = Math::PLU_decomposition(*this);

	double determinantU = 1.0;
	for (size_t i = 0; i < this->order(); ++i)
		determinantU *= this->at(i, i);

	double determinantP = 1;
	for (size_t indx = 0; indx < pivot_indices.size(); ++indx) {
		if (indx + 1 != pivot_indices[indx])		// it means row exchange is occured !
			determinantP *= -1;
	}

	return determinantP * determinantU;		//determinantL = 1
}

double RowMajorMatrix::determinant_Full_Matrix(void) const {
	RowMajorMatrix tmp = *this;
	return tmp.determinant_Full_Matrix();
}

double RowMajorMatrix::determinant_Triangle_Matrix(void) const {
	double determinant = 1.0;
	for (size_t i = 0; i < this->order(); ++i)
		determinant *= this->at(i, i);

	return determinant;
}

void RowMajorMatrix::inspect_Range(const size_t irow, const size_t jcolumn) const {
	if (this->row_ <= irow || this->column_ <= jcolumn)
		FATAL_ERROR("out of range");
}

void RowMajorMatrix::inspect_Size(void) const {
	switch (this->type_) {
	case MatrixType::Full: {
		if (!this->is_Full_Storage_Value())
			FATAL_SIZE_ERROR;
		return;
	}
	case MatrixType::FullLowerTriangle:
	case MatrixType::FullUpperTriangle: {
		if (!this->is_Full_Storage_Value() || !this->is_Square_Matrix())
			FATAL_SIZE_ERROR;
		return;
	}
	case MatrixType::PackedLowerTriangle:
	case MatrixType::PackedUpperTriangle: {
		if (!this->is_Packed_Storage_Value() || !this->is_Square_Matrix())
			FATAL_SIZE_ERROR;
		return;
	}
	default:
		FATAL_TYPE_ERROR;
		return;
	}
}

void RowMajorMatrix::inspect_Value(void) const {
	switch (this->type_)
	{
	case MatrixType::Full:
	case MatrixType::PackedLowerTriangle:
	case MatrixType::PackedUpperTriangle:
		return;
	case MatrixType::FullLowerTriangle: {
		if (!this->is_Lower_Triangle_Value())
			FATAL_ERROR("This matrix have non zero component in upper triangle part");
		return;
	}
	case MatrixType::FullUpperTriangle: {
		if (!this->is_Upper_Triangle_Value())
			FATAL_ERROR("This matrix have non zero component in lower triangle part");
		return;
	}
	default:
		FATAL_TYPE_ERROR;
		return;
	}
}

RowMajorMatrix& RowMajorMatrix::inverse_Full_Matrix(void) {
	const std::vector<int> ipiv = Math::PLU_decomposition(*this);

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const lapack_int n = static_cast<int>(this->order());
	const lapack_int lda = n;
	const lapack_int info = LAPACKE_dgetri(matrix_layout, n, this->data(), lda, ipiv.data());

	if (info > 0)
		FATAL_ERROR("U is singular matrix in L-U decomposition");
	else if (info < 0)
		FATAL_ERROR("fail to inverse the matrix");

	this->update_Type();
	return *this;
}

RowMajorMatrix& RowMajorMatrix::inverse_Full_Triangle_Matrix(void) {
	if (this->is_transposed())
		this->transpose_Value();

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const char uplo = this->LAPACK_Triangle_Type();
	const char diag = this->LAPACK_Diagonal_Type();
	const lapack_int n = static_cast<int>(this->order());
	const lapack_int lda = n;
		
	auto info = LAPACKE_dtrtri(matrix_layout, uplo, diag, n, this->data(), lda);

	if (info > 0)
		FATAL_ERROR("U is singular matrix in L-U decomposition");
	else if (info < 0)
		FATAL_ERROR("fail to inverse the matrix");

	this->update_Type();
	return *this;
}

RowMajorMatrix& RowMajorMatrix::inverse_Packed_Triangle_Matrix(void) {
	if (this->is_transposed())
		this->transpose_Value();

	const int matrix_layout = LAPACK_ROW_MAJOR;
	const char uplo = this->LAPACK_Triangle_Type();
	const char diag = this->LAPACK_Diagonal_Type();
	const lapack_int n = static_cast<int>(this->order());

	auto info = LAPACKE_dtptri(matrix_layout, uplo, diag, n, this->data());

	if (info > 0)
		FATAL_ERROR("U is singular matrix in L-U decomposition");
	else if (info < 0)
		FATAL_ERROR("fail to inverse the matrix");

	this->update_Type();
	return *this;
}

bool RowMajorMatrix::is_Full_Storage_Type(void) const {
	switch (this->type_) {
	case MatrixType::Full:
	case MatrixType::FullLowerTriangle:
	case MatrixType::FullUpperTriangle:
		return true;
	case MatrixType::PackedLowerTriangle:
	case MatrixType::PackedUpperTriangle:
		return false;
	default:
		FATAL_TYPE_ERROR;
		return false;
	}
}

bool RowMajorMatrix::is_Full_Storage_Value(void) const {
	if (this->row_ * this->column_ != this->value_.size())
		return false;
	else
		return true;
}

bool RowMajorMatrix::is_Upper_Triangle_Type(void) const {
	switch (this->type_) {
	case MatrixType::FullUpperTriangle:
	case MatrixType::PackedUpperTriangle:
		return true;
	case MatrixType::Full:
	case MatrixType::FullLowerTriangle:
	case MatrixType::PackedLowerTriangle:
		return false;
	default:
		FATAL_TYPE_ERROR;
		return false;
	}
}

bool RowMajorMatrix::is_Square_Matrix(void) const {
	if (this->row_ == this->column_)
		return true;
	else
		return false;
}

bool RowMajorMatrix::is_transposed(void) const {
	if (this->transpose_type_ == CBLAS_TRANSPOSE::CblasNoTrans)
		return false;
	else
		return true;
}

bool RowMajorMatrix::is_Upper_Triangle_Value(void) const {
	for (size_t i = 1; i < this->row_; ++i)
		for (size_t j = 0; j < i; ++j) {
			if (this->at(i, j) != 0)
				return false;
		}
	return true;
}

bool RowMajorMatrix::is_Lower_Triangle_Value(void) const {
	for (size_t i = 0; i < this->row_; ++i)
		for (size_t j = i + 1; j < this->column_; ++j) {
			if (this->at(i, j) != 0)
				return false;
		}
	return true;
}

bool RowMajorMatrix::is_Packed_Storage_Value(void) const {
	const size_t required_num_value = static_cast<size_t>((this->row_) * (this->row_ + 1) * 0.5);
	if (this->value_.size() == required_num_value)
		return true;
	else
		return false;
}

bool RowMajorMatrix::is_Unit_Diagonal(void) const {
	if (!this->is_Square_Matrix())
		return false;

	for (size_t i = 0; i < this->row_; ++i) {
		if (this->at(i, i) != 1)
			return false;
	}

	return true;
}

bool RowMajorMatrix::is_Lower_Part(const size_t irow, const size_t jcolumn) const {
	if (!this->is_Square_Matrix())
		FATAL_ERROR("Only Square matrix can have lower part of matrix");

	if (irow > jcolumn)
		return true;
	else
		return false;
}

bool RowMajorMatrix::is_Upper_Part(const size_t irow, const size_t jcolumn) const {
	if (!this->is_Square_Matrix())
		FATAL_ERROR("Only Square matrix can have upper part of matrix");

	if (irow < jcolumn)
		return true;
	else
		return false;
}

char RowMajorMatrix::LAPACK_Diagonal_Type(void) const {
	if (this->is_Unit_Diagonal())
		return 'U';
	else
		return 'N';
}

char RowMajorMatrix::LAPACK_Triangle_Type(void) const {
	if (this->is_Upper_Triangle_Type())
		return 'U';
	else
		return 'L';
}

size_t RowMajorMatrix::Leading_Dimension(void) const {
	// use when need original leading dimension
	if (this->is_transposed())
		return this->row_;
	else
		return this->column_;
}

RowMajorMatrix& RowMajorMatrix::multiply_Full_Matrix_Full_Matrix(const RowMajorMatrix& B) { // A = alpha * A * B + beta * A	=>	A *= B
	RowMajorMatrix& A = *this;

	if (A.column_ != B.row_)
		FATAL_SIZE_ERROR;

	const CBLAS_LAYOUT layout = CBLAS_LAYOUT::CblasRowMajor;
	const CBLAS_TRANSPOSE transA = A.transpose_type_;
	const CBLAS_TRANSPOSE transB = B.transpose_type_;
	const MKL_INT m = static_cast<int>(A.row_);
	const MKL_INT n = static_cast<int>(B.column_);
	const MKL_INT k = static_cast<int>(A.column_);
	const double alpha = 1;
	const MKL_INT lda = static_cast<int>(A.Leading_Dimension());
	const MKL_INT ldb = static_cast<int>(B.Leading_Dimension());
	const double beta = 0;
	const MKL_INT ldc = n;

	cblas_dgemm(layout, transA, transB, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, A.data(), ldc);

	A.transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	A.update_Type();
	return A;
}

RowMajorMatrix& RowMajorMatrix::multiply_Full_Matrix_Full_Triangle_Matrix(const RowMajorMatrix& B) { // A = alpha * A * B => A *= B
	// A : Full matrix
	// B : Full triangle matrix
	RowMajorMatrix& A = *this;

	if (A.column_ != B.row_)
		FATAL_SIZE_ERROR;

	if (A.is_transposed())
		A.transpose_Value();

	const CBLAS_LAYOUT layout = CBLAS_LAYOUT::CblasRowMajor;
	const CBLAS_SIDE side = CBLAS_SIDE::CblasRight;
	const CBLAS_UPLO uplo = B.CBLAS_Triangle_Type();
	const CBLAS_TRANSPOSE transB = B.transpose_type_;
	const CBLAS_DIAG diag = B.CBALS_Diagonal_Type();
	const MKL_INT m = static_cast<int>(A.row_);
	const MKL_INT n = static_cast<int>(A.column_);
	const double alpha = 1;
	const MKL_INT lda = n;
	const MKL_INT ldb = n;

	cblas_dtrmm(layout, side, uplo, transB, diag, m, n, alpha, B.data(), lda, A.data(), ldb);

	A.update_Type();
	return A;
}

MathVector RowMajorMatrix::multiply_Full_Matrix_Vector(const MathVector& x) const {// y = alpha Ax + beta y	
	if (this->type_ != MatrixType::Full)
		FATAL_TYPE_ERROR;
	if (this->column_ != x.size())
		FATAL_SIZE_ERROR;

	const auto layout = CBLAS_LAYOUT::CblasRowMajor;
	const auto trans = this->transpose_type_;
	MKL_INT m;
	MKL_INT n;
	if (trans == CBLAS_TRANSPOSE::CblasNoTrans) {
		m = static_cast<int>(this->row_);
		n = static_cast<int>(this->column_);
	}
	else {
		m = static_cast<int>(this->column_);
		n = static_cast<int>(this->row_);
	}
	const double alpha = 1;
	const MKL_INT lda = n;
	const MKL_INT incx = 1;
	const MKL_INT beta = 0;
	const MKL_INT incy = 1;

	MathVector y(this->row_);

	cblas_dgemv(layout, trans, m, n, alpha, this->data(), lda, x.data(), incx, beta, y.data(), incy);

	return y;
}

MathVector RowMajorMatrix::multiply_Full_Triangle_Matrix_Vector(const MathVector& vec) const {//x = Ax 
	if (this->type_ != MatrixType::FullLowerTriangle && this->type_ != MatrixType::FullUpperTriangle)
		FATAL_TYPE_ERROR;
	if (this->column_ != vec.size())
		FATAL_SIZE_ERROR;

	MathVector x(vec.size());

	const CBLAS_LAYOUT layout = CBLAS_LAYOUT::CblasRowMajor;
	const CBLAS_UPLO uplo = this->CBLAS_Triangle_Type();
	const CBLAS_TRANSPOSE trans = this->transpose_type_;
	const CBLAS_DIAG diag = this->CBALS_Diagonal_Type();
	const MKL_INT n = static_cast<int>(this->row_);
	const MKL_INT lda = n;
	const MKL_INT incx = 1;

	cblas_dtrmv(layout, uplo, trans, diag, n, this->data(), lda, x.data(), incx);

	return x;
}

MathVector RowMajorMatrix::multiply_Packed_Triangle_Matrix_Vector(const MathVector& vec) const {//x = Ax 
	if (this->is_Full_Storage_Type())
		FATAL_TYPE_ERROR;
	if (this->column_ != vec.size())
		FATAL_SIZE_ERROR;

	MathVector x(vec.size());

	const CBLAS_LAYOUT layout = CBLAS_LAYOUT::CblasRowMajor;
	const CBLAS_UPLO uplo = this->CBLAS_Triangle_Type();
	const CBLAS_TRANSPOSE trans = this->transpose_type_;
	const CBLAS_DIAG diag = this->CBALS_Diagonal_Type();
	const MKL_INT n = static_cast<int>(this->row_);
	const MKL_INT incx = 1;

	cblas_dtpmv(layout, uplo, trans, diag, n, this->data(), x.data(), incx);

	return x;
}

size_t RowMajorMatrix::order(void) const {
	if (this->is_Square_Matrix())
		return this->row_;
	else {
		FATAL_ERROR("It is not Sqaure Matrix");
		return NULL;
	}
}

RowMajorMatrix& RowMajorMatrix::transpose_Value(void) {
	if (!this->is_Full_Storage_Type())
		this->to_Full_Matrix();

	const char odering = 'R';
	const char trans = 'T';
	const size_t rows = this->row_;
	const size_t cols = this->column_;
	const double alpha = 1;
	const size_t lda = cols;
	const size_t ldb = rows;

	mkl_dimatcopy(odering, trans, rows, cols, alpha, this->data(), lda, ldb);
	this->transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;
	return *this;
}

RowMajorMatrix RowMajorMatrix::transpose_Value(void) const {
	RowMajorMatrix result(*this);
	return result.transpose_Value();
}

RowMajorMatrix& RowMajorMatrix::to_Full_Matrix(void) {
	switch (this->type_) {
	case MatrixType::PackedLowerTriangle: {
		RowMajorMatrix full_matrix(MatrixType::FullLowerTriangle, this->order());

		for (size_t i = 0; i < this->row_; ++i)
			for (size_t j = 0; j <= i; ++j)
				full_matrix.at(i, j) = this->at(i, j);

		return *this = std::move(full_matrix);
	}
	case MatrixType::PackedUpperTriangle: {
		RowMajorMatrix full_matrix(MatrixType::FullUpperTriangle, this->order());

		for (size_t i = 0; i < this->row_; ++i)
			for (size_t j = i; j < this->column_; ++j)
				full_matrix.at(i, j) = this->at(i, j);

		return *this = std::move(full_matrix);
	}
	default:
		FATAL_TYPE_ERROR;
		break;
	}

	return *this;
}

RowMajorMatrix RowMajorMatrix::to_Full_Matrix(void) const {
	RowMajorMatrix result = *this;
	return result.to_Full_Matrix();
}

const RowMajorMatrix& RowMajorMatrix::update_Type(void) const {
	if (!this->is_Full_Storage_Type())
		return *this;

	if (this->is_Square_Matrix()) {
		if (this->is_Lower_Triangle_Value())
			this->type_ = MatrixType::FullLowerTriangle;
		else if (this->is_Upper_Triangle_Value())
			this->type_ = MatrixType::FullUpperTriangle;
		else
			this->type_ = MatrixType::Full;
		return *this;
	}

	this->type_ = MatrixType::Full;
	return *this;
}