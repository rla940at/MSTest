
#include "../INC/Polynomial.h"

namespace Math{
	Polynomial translate(const Monomial& monomial, const MathVector& translation_vector){
		if (monomial.is_Constant())
			return Polynomial(monomial);

		Polynomial result(1.0);		
		for (size_t i = 0; i < monomial.num_variable(); ++i){
			const auto exponent = monomial.exponent(i);

			if (exponent == 0)
				continue;

			std::array<double, 2> coefficient_set = { 1, -translation_vector[i] };
			std::array<Monomial, 2> monomial_set = { Monomial(i),Monomial() };
			Polynomial translated_Xi(coefficient_set, monomial_set);

			for (size_t i = 0; i < exponent; ++i)
				result *= translated_Xi;
		}

		return result;
	}

	Polynomial& translate(Polynomial& polynomial, const MathVector& translation_vector) {
		Polynomial result;
		for (const auto& [monomial, coefficient] : polynomial.monomial_to_coefficient_)
			result += coefficient * Math::translate(monomial, translation_vector);
		polynomial = std::move(result);
		return polynomial;
	}

	Polynomial translate(const Polynomial& polynomial, const MathVector& translation_vector){
		Polynomial result(polynomial);
		return Math::translate(result, translation_vector);
	}

	Polynomial differentiate(const Monomial& monomial, const size_t differential_variable_index){
		if (monomial.is_Constant())
			return Polynomial();

		auto result_monomial = monomial;

		const auto coefficient = static_cast<double>(monomial.exponent(differential_variable_index));
		result_monomial.reduce_Order(differential_variable_index);

		return Polynomial(coefficient, result_monomial);
	}

	Polynomial differentiate(const Polynomial& polynomial, const size_t differential_variable){
		Polynomial result;

		for (const auto& [monomial, coefficient] : polynomial.monomial_to_coefficient_)
			result += coefficient * Math::differentiate(monomial, differential_variable);

		return result;
	}
}


namespace Editor {
	std::string to_String(const Polynomial& polynomial) {
		std::string str;

		str << "{ ";
		const auto start_iter = polynomial.monomial_to_coefficient_.rbegin();
		const auto end_iter = polynomial.monomial_to_coefficient_.rend();
		for (auto iter = start_iter; iter != end_iter; ++iter) {
			const auto& [monomial, coefficient] = *iter;
			str << coefficient << " X " << Editor::to_String(monomial) << " + ";
		}
		Editor::erase_back(str, 3);
		str << " }";

		return str;
	}
}


std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial) {
	return os << Editor::to_String(polynomial);
}

Polynomial operator*(const double scalar, const Polynomial& polynomial){
	return polynomial.operator*(scalar);
};


Polynomial& Polynomial::operator*=(const double scalar) {
	if (this->is_Zero())
		return *this;

	if (scalar == 0.0) {
		*this = Polynomial();
		return *this;
	}

	for (auto& [monomial, coefficient] : monomial_to_coefficient_)
		coefficient *= scalar;

	return *this;
}

Polynomial Polynomial::operator*(const double scalar) const{
	if (scalar == 0)	
		return Polynomial();

	Polynomial result(*this);
	return result *= scalar;
}

Polynomial& Polynomial::operator*=(const Polynomial& other) {
	if (this->is_Zero() || other.is_Zero())
		return *this = Polynomial();

	Polynomial result;
	for (const auto& [this_monomial, this_coefficient] : this->monomial_to_coefficient_)
		for (const auto& [other_monomial, other_coefficient] : other.monomial_to_coefficient_) {
			const auto result_coefficient = this_coefficient * other_coefficient;
			auto result_monomial = this_monomial * other_monomial;
			result.insert(result_coefficient, std::move(result_monomial));
		}
	return *this = std::move(result);
}

Polynomial Polynomial::operator*(const Polynomial& other) const{
	Polynomial result(*this);
	return result *= other;
}

Polynomial& Polynomial::operator+=(const Polynomial& other){
	for (const auto& [other_monomial, other_coefficient] : other.monomial_to_coefficient_)
		this->insert(other_coefficient, other_monomial);

	return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& other){
	for (const auto& [other_monomial, other_coefficient] : other.monomial_to_coefficient_)
		this->insert(-other_coefficient, other_monomial);

	return *this;
}

double Polynomial::operator()(const MathVector& variable_vector) const{
	double result = 0.0;
	for (const auto& [monomial, coefficient] : this->monomial_to_coefficient_)
		result += coefficient * monomial(variable_vector);

	return result;
}

void Polynomial::insert(const double coefficient, const Monomial& monomial){
	if (coefficient == 0) 
		return;

	auto iter = this->monomial_to_coefficient_.find(monomial);
	if (iter == this->monomial_to_coefficient_.end())
		this->monomial_to_coefficient_.emplace(monomial, coefficient);
	else{
		iter->second += coefficient;
		if (iter->second == 0)
			this->monomial_to_coefficient_.erase(iter);
	}
}

void Polynomial::insert(const double coefficient, Monomial&& monomial) {
	if (coefficient == 0)
		return;

	auto iter = this->monomial_to_coefficient_.find(monomial);
	if (iter == this->monomial_to_coefficient_.end())
		this->monomial_to_coefficient_.emplace(std::move(monomial), coefficient);
	else {
		iter->second += coefficient;
		if (iter->second == 0)
			this->monomial_to_coefficient_.erase(iter);
	}
}

bool Polynomial::is_Zero() const{
	if (this->monomial_to_coefficient_.size() != 1)
		return false;

	const auto& [monomial, coefficient] = *this->monomial_to_coefficient_.begin();

	if (coefficient == 0)
		return true;
	else
		return false;
}

size_t Polynomial::order(void) const{
	const auto& [highest_order_monomial, coefficient] = *this->monomial_to_coefficient_.rbegin();
	return highest_order_monomial.order();
}