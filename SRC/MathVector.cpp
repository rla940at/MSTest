
#include "../INC/MathVector.h"

namespace Math{
	MathVector& abs(MathVector& x) {
		const MKL_INT n = static_cast<int>(x.size());
		vdAbs(n, x.data(), x.data());
		return x;
	}

	MathVector abs(const MathVector& x) {
		MathVector result = x;
		return Math::abs(result);
	}

	MathVector abs(MathVector&& x) {
		MathVector result = std::move(x);
		return Math::abs(result);
	}

	MathVector& sqrt(MathVector& x) {
		const MKL_INT n = static_cast<int>(x.size());
		vdSqrt(n, x.data(), x.data());
		return x;
	}

	MathVector sqrt(const MathVector& x) {
		MathVector result = x;
		return Math::sqrt(result);
	}

	double L2_norm(const MathVector& x) {
		const MKL_INT n = static_cast<int>(x.size());
		const MKL_INT incx = 1;
		return cblas_dnrm2(n, x.data(), incx);
	}

	double inner_product(const MathVector& x, const MathVector& y) {
		if (x.size() != y.size())
			FATAL_SIZE_ERROR;

		const MKL_INT n = static_cast<int>(x.size());
		const MKL_INT incx = 1;
		const MKL_INT incy = 1;
		return cblas_ddot(n, x.data(), incx, y.data(), incy);
	}

	MathVector& normalize(MathVector& x) {
		return x *= 1.0 / Math::L2_norm(x);
	}

	MathVector normalize(const MathVector& x) {
		MathVector result = x;
		return Math::normalize(result);
	}

	MathVector normalize(MathVector&& x) {
		MathVector result = std::move(x);
		return Math::normalize(result);
	}
}


namespace Editor
{
	void merge(std::vector<double>& x, MathVector&& y) {
		x.insert(x.end(), std::make_move_iterator(y.begin()), std::make_move_iterator(y.end()));
	}

	std::string to_String(const MathVector& x) {
		std::string result;
		for (const auto& value : x.value_set_)
			result << "\t" << value << "\n";
		return result;
	}
}


std::ostream& operator<<(std::ostream& os, const MathVector& x) {
	return os << Editor::to_String(x);
}


MathVector MathVector::operator-(const MathVector& y) const{
	MathVector x = *this;
	return x -= y;
}

MathVector MathVector::operator+(const MathVector& y) const{
	MathVector x = *this;
	return x += y;
}

MathVector& MathVector::operator+=(const MathVector& y) { //x = x + y
	MathVector& x = *this;

	if (x.size() != y.size())
		FATAL_SIZE_ERROR;

	const MKL_INT n = static_cast<int>(x.size());
	vdAdd(n, x.data(), y.data(), x.data());

	return x;
}

MathVector& MathVector::operator-=(const MathVector& y) { //x = x - y
	MathVector& x = *this;

	if (x.size() != y.size())
		FATAL_SIZE_ERROR;

	const MKL_INT n = static_cast<int>(x.size());
	vdSub(n, x.data(), y.data(), x.data());

	return x;
}

MathVector& MathVector::operator*=(const double scalar) { //x = a * x	
	MathVector& x = *this;

	const MKL_INT n = static_cast<int>(x.size());
	const double a = scalar;
	const MKL_INT incx = 1;
	cblas_dscal(n, a, x.data(), incx);

	return x;
}

MathVector MathVector::operator*(const double scalar) const {
	MathVector x = *this;
	return x *= scalar;
}

MathVector& MathVector::operator*=(const MathVector& y) { //x_i *= y_i value by value
	MathVector& x = *this;

	if (x.size() != y.size())
		FATAL_SIZE_ERROR;

	const MKL_INT n = static_cast<int>(x.size());
	vdMul(n, x.data(), y.data(), x.data());

	return x;
}

MathVector MathVector::operator*(const MathVector& y) const {
	MathVector x = *this;
	return x *= y;
}
