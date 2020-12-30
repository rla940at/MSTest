#pragma once

#include "Monomial.h"

#include <array>
#include <map>

class Polynomial;
namespace Math {
	Polynomial differentiate(const Monomial& monomial, const size_t differential_variable);

	Polynomial differentiate(const Polynomial& polynomial, const size_t differential_variable);

	Polynomial translate(const Monomial& monomial, const MathVector& translation_vector);

	Polynomial& translate(Polynomial& polynomial, const MathVector& translation_vector);

	Polynomial translate(const Polynomial& polynomial, const MathVector& translation_vector);
}


namespace Editor {
	std::string to_String(const Polynomial& polynomial);
}


std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial);

Polynomial operator*(const double scalar, const Polynomial& polynomial);


class Polynomial
{
private:
	std::map<Monomial, double> monomial_to_coefficient_;

private:
	friend std::string Editor::to_String(const Polynomial& polynomial);

	friend Polynomial Math::differentiate(const Polynomial& polynomial, const size_t differential_variable);

	friend Polynomial& Math::translate(Polynomial& polynomial, const MathVector& translation_vector);

public:
	explicit Polynomial(void)
		: monomial_to_coefficient_{ {Monomial(), 0} } {};

	explicit Polynomial(const double coefficient)
		: monomial_to_coefficient_{ {Monomial(), coefficient} } {};

	explicit Polynomial(const Monomial& monomial)
		: monomial_to_coefficient_{ {monomial, 1} } {};

	explicit Polynomial(const double coefficient, const Monomial& monomial)
		: monomial_to_coefficient_{ {monomial, coefficient} } {};

	template <typename coeffcient_container, typename monomial_container>
	explicit Polynomial(const coeffcient_container& coefficient_set, const monomial_container& monomial_set) {
		if (coefficient_set.size() != monomial_set.size())
			FATAL_SIZE_ERROR;

		for (size_t i = 0; i < coefficient_set.size(); ++i)
			this->insert(coefficient_set[i], monomial_set[i]);
	};

	Polynomial& operator*=(const double scalar);

	Polynomial operator*(const double scalar) const;

	Polynomial operator*(const Polynomial& other) const;

	Polynomial& operator+=(const Polynomial& other);

	Polynomial& operator-=(const Polynomial& other);

	Polynomial& operator*=(const Polynomial& other);

	double operator()(const MathVector& variable_vector) const;

	void insert(const double coefficient, const Monomial& monomial);

	void insert(const double coefficient, Monomial&& monomial);

	bool is_Zero(void) const;

	size_t order(void) const;
};