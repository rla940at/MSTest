#pragma once

#include "MathVector.h"

class Monomial;
namespace Editor {
	std::string to_String(const Monomial& monomial);
}


std::ostream& operator<<(std::ostream& ostream, const Monomial& monomial);


class Monomial
{
private:
	std::vector<size_t> exponent_set_;

private:
	friend std::string Editor::to_String(const Monomial& monomial);

public:
	Monomial(void) {
		this->exponent_set_.push_back(0);
	};

	template<typename T>
	Monomial(const T argument) { // Monomial(0) => x / Monomial(2) => z
		const auto variable_index = static_cast<size_t>(argument);
		this->exponent_set_.resize(variable_index + 1);
		this->exponent_set_.back() = 1;
	};

	template<typename ...Args>
	Monomial(const Args... arguments)
		: exponent_set_{ static_cast<size_t>(arguments)... } {};

	Monomial(std::initializer_list<size_t> list)
		: exponent_set_{ list } {}; // Monomial{0} => 1 / Monomial{2} = x^2

	Monomial(std::vector<size_t>&& exponent_set)
		: exponent_set_(std::move(exponent_set)) {};


	Monomial operator*(const Monomial& other) const;

	double operator()(const MathVector& variable_vector) const;

	bool operator<(const Monomial& other) const;


	size_t exponent(size_t variable_index) const;

	bool is_Constant(void) const;

	size_t num_variable(void) const {
		return this->exponent_set_.size();
	};

	size_t order(void) const;

	void reduce_Order(const size_t variable_index);
};