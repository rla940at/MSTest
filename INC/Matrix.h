#pragma once

#include "MathVector.h"

#include <array>

class RowMajorMatrix;
namespace Math {
	RowMajorMatrix& abs(RowMajorMatrix& A);

	RowMajorMatrix cofactor_matrix(const RowMajorMatrix& A);

	std::vector<int> PLU_decomposition(RowMajorMatrix& A);

	std::array<RowMajorMatrix, 3> singular_value_decomposition(RowMajorMatrix& A);

	double determinant(RowMajorMatrix&& A);

	double determinant(const RowMajorMatrix& A);

	double Frobenius_Norm(const RowMajorMatrix& A);

	RowMajorMatrix& inverse(RowMajorMatrix& A);

	RowMajorMatrix inverse(const RowMajorMatrix& A);

	double minor(const RowMajorMatrix& A, const size_t i, const size_t j);

	void Moore_Penrose_Inverse(RowMajorMatrix& A);

	RowMajorMatrix& transpose(RowMajorMatrix& A);

	RowMajorMatrix transpose(const RowMajorMatrix& A);

	RowMajorMatrix transpose(RowMajorMatrix&& A);
}


namespace Editor {
	std::string to_String(const RowMajorMatrix& A);
}


std::ostream& operator<<(std::ostream& os, const RowMajorMatrix& A);

RowMajorMatrix operator*(const double scalar, const RowMajorMatrix& A);


enum class MatrixType
{
	Full,
	FullUpperTriangle,
	FullLowerTriangle,
	PackedUpperTriangle,
	PackedLowerTriangle
};


class RowMajorMatrix
{
private:
	CBLAS_TRANSPOSE transpose_type_ = CBLAS_TRANSPOSE::CblasNoTrans;

	mutable MatrixType type_ = MatrixType::Full;

	size_t row_ = 0;

	size_t column_ = 0;

	MathVector value_;

private:
	friend RowMajorMatrix& Math::abs(RowMajorMatrix& A);

	friend double Math::determinant(RowMajorMatrix&& A);

	friend double Math::Frobenius_Norm(const RowMajorMatrix& A);

	friend std::vector<int> Math::PLU_decomposition(RowMajorMatrix& A);

	friend std::array<RowMajorMatrix, 3> Math::singular_value_decomposition(RowMajorMatrix& A);

	friend RowMajorMatrix& Math::inverse(RowMajorMatrix& A);

	friend void Math::Moore_Penrose_Inverse(RowMajorMatrix& A);

	friend RowMajorMatrix& Math::transpose(RowMajorMatrix& A);

	friend std::string Editor::to_String(const RowMajorMatrix& A);

public:
	explicit RowMajorMatrix() = default;

	explicit RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column);

	explicit RowMajorMatrix(const MatrixType matrix_type, const size_t matrix_order);

	RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column, const MathVector& value);

	RowMajorMatrix(const MatrixType matrix_type, const size_t num_row, const size_t num_column, MathVector&& value);

	RowMajorMatrix& operator+=(const RowMajorMatrix& other);

	RowMajorMatrix operator+(const RowMajorMatrix& other) const;

	RowMajorMatrix& operator-=(const RowMajorMatrix& B);

	RowMajorMatrix operator-(const RowMajorMatrix& other) const;

	RowMajorMatrix& operator*=(const double scalar);

	RowMajorMatrix operator*(const double scalar) const;
	
	MathVector operator*(const MathVector& vec) const;

	RowMajorMatrix& operator*=(const RowMajorMatrix& other);

	RowMajorMatrix operator*(const RowMajorMatrix& other) const;

	bool operator==(const RowMajorMatrix& B) const;

	bool operator!=(const RowMajorMatrix& B) const {
		return !(*this == B);
	};


	double& at(const size_t row, const size_t column);

	double at(const size_t row, const size_t column) const;

	void change_Column(const size_t column_index, const MathVector& vector);

	void change_Row(const size_t row_index, const MathVector& vector);

	RowMajorMatrix& change_Value(MathVector&& value) {
		this->value_ = std::move(value);
		return *this;
	};

	MathVector column(const size_t column_index) const;

	RowMajorMatrix column(const size_t start_column, const size_t end_column) const;

	RowMajorMatrix part(const size_t start_row, const size_t end_row, const size_t start_column, const size_t end_column) const;

	RowMajorMatrix part(const size_t part_matrix_order) const;

	RowMajorMatrix part(const size_t deleting_row_index, const size_t deleting_column_index) const;

	MathVector row(const size_t row_index) const;

	RowMajorMatrix row(const size_t start_row, const size_t end_row) const;

	std::pair<size_t, size_t> size(void) const {
		return { this->row_, this->column_ };
	};

	RowMajorMatrix& scale_Row(const size_t start_row, const size_t end_row, const double scale_factor);

	RowMajorMatrix& scale_Column(const size_t start_column, const size_t end_column, const double scale_factor);

private:
	CBLAS_DIAG CBALS_Diagonal_Type(void) const;

	CBLAS_UPLO CBLAS_Triangle_Type(void) const;

	double* data(void) {
		return this->value_.data();
	};

	const double* data(void) const {
		return this->value_.data();
	};

	double determinant_Full_Matrix(void);

	double determinant_Full_Matrix(void) const;

	double determinant_Triangle_Matrix(void) const;

	void inspect_Range(const size_t irow, const size_t jcolumn) const;

	void inspect_Size(void) const;

	void inspect_Value(void) const;

	RowMajorMatrix& inverse_Full_Matrix(void);

	RowMajorMatrix& inverse_Full_Triangle_Matrix(void);

	RowMajorMatrix& inverse_Packed_Triangle_Matrix(void);

	bool is_Full_Storage_Type(void) const;

	bool is_Full_Storage_Value(void) const;

	bool is_Lower_Triangle_Value(void) const;

	bool is_Lower_Part(const size_t irow, const size_t jcolumn) const;
		
	bool is_Packed_Storage_Value(void) const;

	bool is_Square_Matrix(void) const;

	bool is_transposed(void) const;

	bool is_Upper_Triangle_Value(void) const;

	bool is_Upper_Triangle_Type(void) const;

	bool is_Unit_Diagonal(void) const;
	
	bool is_Upper_Part(const size_t irow, const size_t jcolumn) const;

	char LAPACK_Diagonal_Type(void) const;

	char LAPACK_Triangle_Type(void) const;

	size_t Leading_Dimension(void) const;

	RowMajorMatrix& multiply_Full_Matrix_Full_Matrix(const RowMajorMatrix& B);

	RowMajorMatrix& multiply_Full_Matrix_Full_Triangle_Matrix(const RowMajorMatrix& B);

	MathVector multiply_Full_Matrix_Vector(const MathVector& vec) const;
	
	MathVector multiply_Full_Triangle_Matrix_Vector(const MathVector& vec) const;

	MathVector multiply_Packed_Triangle_Matrix_Vector(const MathVector& vec) const;

	size_t order(void) const;

	RowMajorMatrix& transpose_Value(void);

	RowMajorMatrix transpose_Value(void) const;

	RowMajorMatrix& to_Full_Matrix(void);

	RowMajorMatrix to_Full_Matrix(void) const;

	const RowMajorMatrix& update_Type(void) const;

	













};