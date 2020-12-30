#pragma once

#include "Matrix.h"
#include "Polynomial.h"

template <typename T>
class VectorFunction;
template <typename T>
class JacobianMatrix;

namespace Editor 
{
	template <typename T>
	std::string to_String(const JacobianMatrix<T>& jacobian_matrix){
		std::string str;

		const auto [num_row, num_column] = jacobian_matrix.size();
		for (size_t i = 0; i < num_row; ++i)
		for (size_t j = 0; j < num_column; ++j)
			str << "[" << i << "," << j << "]  :  " << jacobian_matrix.at(i, j) << "\n";

		return str;
	}
}

namespace Math {
	template <typename T>
	VectorFunction<T> gradient(const T& function, const size_t num_variable) {
		VectorFunction<T> gradient;
		for (size_t i = 0; i < num_variable; ++i)
			gradient.emplace_back(Math::differentiate(function, i));
		
		return gradient;
	}

	template <typename T>
	JacobianMatrix<T> jacobian(const VectorFunction<T>& vector_function, const size_t num_variable) {
		const auto num_function = vector_function.size();

		JacobianMatrix<T> jacobian_matrix(num_function, num_variable);
		for (size_t i = 0; i < num_function; ++i)
		for (size_t j = 0; j < num_variable; ++j)
			jacobian_matrix.at(i, j) = Math::differentiate(vector_function.at(i), j);

		return jacobian_matrix;
	}
}

template <typename T>
class VectorFunction
{
private:
	std::vector<T> function_set_;

public:
	VectorFunction(void) = default;

	VectorFunction(std::initializer_list<T> list) : function_set_(list) {};

	MathVector operator()(const MathVector& variable_vector) const {
		MathVector result;
		result.reserve(this->size());
		for (const auto& function : this->function_set_)
			result.emplace_back(function(variable_vector));
		return result;
	}

	const T& at(const size_t index) const {
		return this->function_set_[index];
	};

	template <typename ... Args>
	void emplace_back(Args&&... arguments) {
		this->function_set_.emplace_back(std::move(arguments...));
	}

	void emplace_back(T&& function) {
		this->function_set_.emplace_back(function);
	}

	size_t size(void) const { return function_set_.size(); };
};

template <typename T>
class JacobianMatrix
{
private:
	std::vector<std::vector<T>> function_set_;

public:
	JacobianMatrix(void) = default;

	JacobianMatrix(const size_t num_function, const size_t num_variable) {
		function_set_.resize(num_function);
		for (auto& set : function_set_)
			set.resize(num_variable);
	}


	RowMajorMatrix operator()(const MathVector& variable_vector) const {
		const auto [num_row, num_column] = this->size();

		RowMajorMatrix value(MatrixType::Full, num_row, num_column);
		for (size_t i = 0; i < num_row; ++i)
		for (size_t j = 0; j < num_column; ++j)
			value.at(i, j) = function_set_[i][j](variable_vector);

		return value;
	}

	T& at(const size_t i_index, const size_t j_index) {
		return function_set_[i_index][j_index];
	}

	const T& at(const size_t i_index, const size_t j_index) const {
		return function_set_[i_index][j_index];
	}

	std::pair<size_t, size_t> size(void) const {
		return { this->function_set_.size(), this->function_set_.front().size() }; 
	};
};