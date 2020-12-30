#pragma once

#include "StringEditor.h"

#include <mkl.h>
#include <vector>

class MathVector;
namespace Math{
	MathVector& abs(MathVector& x);

	MathVector abs(const MathVector& x);
		
	MathVector abs(MathVector&& x);

	MathVector sqrt(const MathVector& x);

	MathVector& sqrt(MathVector& x);

	double L2_norm(const MathVector& x);

	double inner_product(const MathVector& x, const MathVector& y);

	MathVector& normalize(MathVector& x);

	MathVector normalize(const MathVector& x);

	MathVector normalize(MathVector&& x);
}


namespace Editor{	
	void merge(std::vector<double>& x, MathVector&& y);

	std::string to_String(const MathVector& x);
}


std::ostream& operator<<(std::ostream& os, const MathVector& x);


class MathVector	
{
	using Iter = std::vector<double>::iterator;
	using CIter = std::vector<double>::const_iterator;

private:
	std::vector<double> value_set_;

private:
	friend std::string Editor::to_String(const MathVector& x);

public:
	explicit MathVector(void) = default;

	explicit MathVector(const size_t num_value) {
		value_set_.resize(num_value); 
	};

	MathVector(std::initializer_list<double> list)
		: value_set_{ list } {};

	MathVector(const std::vector<double>& value_set)
		: value_set_(value_set) {};

	explicit MathVector(std::vector<double>&& value_set)
		: value_set_(std::move(value_set)) {};

	template <typename InputIterator>
	explicit MathVector(InputIterator first, InputIterator last)
		: value_set_(first, last) {};	


	MathVector operator-(const MathVector& y) const;

	MathVector operator+(const MathVector& y) const;

	MathVector& operator*=(const double scalar);

	MathVector operator*(const double scalar) const;

	MathVector& operator*=(const MathVector& y);

	MathVector operator*(const MathVector& y) const;

	MathVector& operator-=(const MathVector& y);

	MathVector& operator+=(const MathVector& y);

	bool operator==(const MathVector& y) const {
		return this->value_set_ == y.value_set_;		
	}

	bool operator!=(const MathVector& y) const {
		return this->value_set_ != y.value_set_;
	}

	double& operator[](const size_t index) {
		return value_set_[index]; 
	};

	const double& operator[](const size_t index) const {
		return value_set_[index]; 
	};

	MathVector& operator<<(const double scalar) {
		this->value_set_.push_back(scalar);
		return *this; 
	};


	Iter begin(void) {
		return value_set_.begin();
	};

	CIter begin(void) const {
		return value_set_.begin();
	};

	const double& back(void) const {
		return this->value_set_.back();
	};

	double* data(void) {
		return value_set_.data();
	};

	const double* data(void) const {
		return value_set_.data();
	};

	Iter end(void) {
		return value_set_.end();
	};

	CIter end(void) const {
		return value_set_.begin();
	};

	void emplace_back(const double val) {
		value_set_.emplace_back(val);
	};

	double& front(void) {
		return this->value_set_.front();
	};

	const double& front(void) const {
		return this->value_set_.front();
	};

	void pop_back(void) {
		value_set_.pop_back();
	};

	void push_back(const double val) {
		this->value_set_.push_back(val);
	};

	void reserve(const size_t required_memory) {
		value_set_.reserve(required_memory);
	};

	void resize(const size_t required_size) {
		value_set_.resize(required_size);
	};

	size_t size(void) const {
		return value_set_.size();
	};
			
	const std::vector<double>& get_Value(void) const {
		return value_set_;
	};
};