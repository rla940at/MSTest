#pragma once

#include "../INC/Matrix.h"

void MatrixTest1(void) {
	std::cout << " operator *= test \n\n";
	std::vector<double> v = { 1,0,3,4 };

	std::vector<double> v2 = { 1,3,4 };

	const RowMajorMatrix m1(MatrixType::Full, 2, 2, v);
	const auto m1t = Math::transpose(m1);

	const RowMajorMatrix m2(MatrixType::FullLowerTriangle, 2, 2, v);
	const auto m2t = Math::transpose(m2);

	const RowMajorMatrix m3(MatrixType::PackedLowerTriangle, 2, 2, v2);
	const auto m3t = Math::transpose(m3);

	auto ref1 = m1 * m1;
	auto ref2 = m1t * m1;
	auto ref3 = m1t * m1t;
	auto ref4 = m1 * m1t;

	std::cout << std::boolalpha;
	std::cout << "\nfull- full\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 *= m1)) << "\n";
		std::cout << (ref2 == (f2 *= m1)) << "\n";
		std::cout << (ref3 == (f3 *= m1t)) << "\n";
		std::cout << (ref4 == (f4 *= m1t)) << "\n";
	}
	std::cout << "\nfull - full triangle\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 *= m2)) << "\n";
		std::cout << (ref2 == (f2 *= m2)) << "\n";
		std::cout << (ref3 == (f3 *= m2t)) << "\n";
		std::cout << (ref4 == (f4 *= m2t)) << "\n";
	}
	std::cout << "\nfull - pack triangle\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 *= m3)) << "\n";
		std::cout << (ref2 == (f2 *= m3)) << "\n";
		std::cout << (ref3 == (f3 *= m3t)) << "\n";
		std::cout << (ref4 == (f4 *= m3t)) << "\n";
	}


	std::cout << "\nfull triangle - full\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 *= m1)) << "\n";
		std::cout << (ref2 == (f2 *= m1)) << "\n";
		std::cout << (ref3 == (f3 *= m1t)) << "\n";
		std::cout << (ref4 == (f4 *= m1t)) << "\n";
	}
	std::cout << "\nfull triangle - full triangle\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 *= m2)) << "\n";
		std::cout << (ref2 == (f2 *= m2)) << "\n";
		std::cout << (ref3 == (f3 *= m2t)) << "\n";
		std::cout << (ref4 == (f4 *= m2t)) << "\n";
	}
	std::cout << "\nfull triangle - pack triangle\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 *= m3)) << "\n";
		std::cout << (ref2 == (f2 *= m3)) << "\n";
		std::cout << (ref3 == (f3 *= m3t)) << "\n";
		std::cout << (ref4 == (f4 *= m3t)) << "\n";
	}


	std::cout << "\npack triangle - full\n"; {
		auto f1 = m3;
		auto f2 = m3t;
		auto f3 = m3t;
		auto f4 = m3;

		std::cout << (ref1 == (f1 *= m1)) << "\n";
		std::cout << (ref2 == (f2 *= m1)) << "\n";
		std::cout << (ref3 == (f3 *= m1t)) << "\n";
		std::cout << (ref4 == (f4 *= m1t)) << "\n";
	}
	std::cout << "\npack triangle - full triangle\n"; {
		auto f1 = m3;
		auto f2 = m3t;
		auto f3 = m3t;
		auto f4 = m3;

		std::cout << (ref1 == (f1 *= m2)) << "\n";
		std::cout << (ref2 == (f2 *= m2)) << "\n";
		std::cout << (ref3 == (f3 *= m2t)) << "\n";
		std::cout << (ref4 == (f4 *= m2t)) << "\n";
	}
	std::cout << "\npack triangle - pack triangle\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 *= m3)) << "\n";
		std::cout << (ref2 == (f2 *= m3)) << "\n";
		std::cout << (ref3 == (f3 *= m3t)) << "\n";
		std::cout << (ref4 == (f4 *= m3t)) << "\n";
	}
}

void MatrixTest2(void) {
	std::cout << "operator * test \n\n";
	std::vector<double> v = { 1,0,3,4 };

	std::vector<double> v2 = { 1,3,4 };

	const RowMajorMatrix m1(MatrixType::Full, 2, 2, v);
	const auto m1t = Math::transpose(m1);

	const RowMajorMatrix m2(MatrixType::FullLowerTriangle, 2, 2, v);
	const auto m2t = Math::transpose(m2);

	const RowMajorMatrix m3(MatrixType::PackedLowerTriangle, 2, 2, v2);
	const auto m3t = Math::transpose(m3);

	auto ref1 = m1 * m1;
	auto ref2 = m1t * m1;
	auto ref3 = m1t * m1t;
	auto ref4 = m1 * m1t;

	std::cout << std::boolalpha;

	std::cout << "\nfull- full\n";
	std::cout << (ref1 == m1 * m1) << "\n";
	std::cout << (ref2 == m1t * m1) << "\n";
	std::cout << (ref3 == m1t * m1t) << "\n";
	std::cout << (ref4 == m1 * m1t) << "\n";

	std::cout << "\nfull- full triangle\n";
	std::cout << (ref1 == m1 * m2) << "\n";
	std::cout << (ref2 == m1t * m2) << "\n";
	std::cout << (ref3 == m1t * m2t) << "\n";
	std::cout << (ref4 == m1 * m2t) << "\n";

	std::cout << "\nfull- packed triangle\n";
	std::cout << (ref1 == m1 * m3) << "\n";
	std::cout << (ref2 == m1t * m3) << "\n";
	std::cout << (ref3 == m1t * m3t) << "\n";
	std::cout << (ref4 == m1 * m3t) << "\n";


	std::cout << "\nfull triangle- full\n";
	std::cout << (ref1 == m2 * m1) << "\n";
	std::cout << (ref2 == m2t * m1) << "\n";
	std::cout << (ref3 == m2t * m1t) << "\n";
	std::cout << (ref4 == m2 * m1t) << "\n";

	std::cout << "\nfull triangle - full triangle\n";
	std::cout << (ref1 == m2 * m2) << "\n";
	std::cout << (ref2 == m2t * m2) << "\n";
	std::cout << (ref3 == m2t * m2t) << "\n";
	std::cout << (ref4 == m2 * m2t) << "\n";

	std::cout << "\nfull triangle- packed triangle\n";
	std::cout << (ref1 == m2 * m3) << "\n";
	std::cout << (ref2 == m2t * m3) << "\n";
	std::cout << (ref3 == m2t * m3t) << "\n";
	std::cout << (ref4 == m2 * m3t) << "\n";


	std::cout << "\npack triangle- full\n";
	std::cout << (ref1 == m3 * m1) << "\n";
	std::cout << (ref2 == m3t * m1) << "\n";
	std::cout << (ref3 == m3t * m1t) << "\n";
	std::cout << (ref4 == m3 * m1t) << "\n";

	std::cout << "\npack triangle - full triangle\n";
	std::cout << (ref1 == m3 * m2) << "\n";
	std::cout << (ref2 == m3t * m2) << "\n";
	std::cout << (ref3 == m3t * m2t) << "\n";
	std::cout << (ref4 == m3 * m2t) << "\n";

	std::cout << "\npack triangle- packed triangle\n";
	std::cout << (ref1 == m3 * m3) << "\n";
	std::cout << (ref2 == m3t * m3) << "\n";
	std::cout << (ref3 == m3t * m3t) << "\n";
	std::cout << (ref4 == m3 * m3t) << "\n";
}

void MatrixTest3(void) {
	std::cout << "test operator += \n\n";

	std::vector<double> v = { 1,0,3,4 };

	std::vector<double> v2 = { 1,3,4 };

	const RowMajorMatrix m1(MatrixType::Full, 2, 2, v);
	const auto m1t = Math::transpose(m1);

	const RowMajorMatrix m2(MatrixType::FullLowerTriangle, 2, 2, v);
	const auto m2t = Math::transpose(m2);

	const RowMajorMatrix m3(MatrixType::PackedLowerTriangle, 2, 2, v2);
	const auto m3t = Math::transpose(m3);

	auto ref1 = m1 + m1;
	auto ref2 = m1t + m1;
	auto ref3 = m1t + m1t;
	auto ref4 = m1 + m1t;

	std::cout << std::boolalpha;
	std::cout << "\nfull- full\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 += m1)) << "\n";
		std::cout << (ref2 == (f2 += m1)) << "\n";
		std::cout << (ref3 == (f3 += m1t)) << "\n";
		std::cout << (ref4 == (f4 += m1t)) << "\n";
	}
	std::cout << "\nfull - full triangle\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 += m2)) << "\n";
		std::cout << (ref2 == (f2 += m2)) << "\n";
		std::cout << (ref3 == (f3 += m2t)) << "\n";
		std::cout << (ref4 == (f4 += m2t)) << "\n";
	}
	std::cout << "\nfull - pack triangle\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 += m3)) << "\n";
		std::cout << (ref2 == (f2 += m3)) << "\n";
		std::cout << (ref3 == (f3 += m3t)) << "\n";
		std::cout << (ref4 == (f4 += m3t)) << "\n";
	}


	std::cout << "\nfull triangle - full\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 += m1)) << "\n";
		std::cout << (ref2 == (f2 += m1)) << "\n";
		std::cout << (ref3 == (f3 += m1t)) << "\n";
		std::cout << (ref4 == (f4 += m1t)) << "\n";
	}
	std::cout << "\nfull triangle - full triangle\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 += m2)) << "\n";
		std::cout << (ref2 == (f2 += m2)) << "\n";
		std::cout << (ref3 == (f3 += m2t)) << "\n";
		std::cout << (ref4 == (f4 += m2t)) << "\n";
	}
	std::cout << "\nfull triangle - pack triangle\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 += m3)) << "\n";
		std::cout << (ref2 == (f2 += m3)) << "\n";
		std::cout << (ref3 == (f3 += m3t)) << "\n";
		std::cout << (ref4 == (f4 += m3t)) << "\n";
	}


	std::cout << "\npack triangle - full\n"; {
		auto f1 = m3;
		auto f2 = m3t;
		auto f3 = m3t;
		auto f4 = m3;

		std::cout << (ref1 == (f1 += m1)) << "\n";
		std::cout << (ref2 == (f2 += m1)) << "\n";
		std::cout << (ref3 == (f3 += m1t)) << "\n";
		std::cout << (ref4 == (f4 += m1t)) << "\n";
	}
	std::cout << "\npack triangle - full triangle\n"; {
		auto f1 = m3;
		auto f2 = m3t;
		auto f3 = m3t;
		auto f4 = m3;

		std::cout << (ref1 == (f1 += m2)) << "\n";
		std::cout << (ref2 == (f2 += m2)) << "\n";
		std::cout << (ref3 == (f3 += m2t)) << "\n";
		std::cout << (ref4 == (f4 += m2t)) << "\n";
	}
	std::cout << "\npack triangle - pack triangle\n"; {
		auto f1 = m3;
		auto f2 = m3t;
		auto f3 = m3t;
		auto f4 = m3;

		std::cout << (ref1 == (f1 += m3)) << "\n";
		std::cout << (ref2 == (f2 += m3)) << "\n";
		std::cout << (ref3 == (f3 += m3t)) << "\n";
		std::cout << (ref4 == (f4 += m3t)) << "\n";
	}
}

void MatrixTest4(void) {
	std::cout << "test operator + \n\n";
	std::vector<double> v = { 1,0,3,4 };

	std::vector<double> v2 = { 1,3,4 };

	const RowMajorMatrix m1(MatrixType::Full, 2, 2, v);
	const auto m1t = Math::transpose(m1);

	const RowMajorMatrix m2(MatrixType::FullLowerTriangle, 2, 2, v);
	const auto m2t = Math::transpose(m2);

	const RowMajorMatrix m3(MatrixType::PackedLowerTriangle, 2, 2, v2);
	const auto m3t = Math::transpose(m3);

	auto ref1 = m1 + m1;
	auto ref2 = m1t + m1;
	auto ref3 = m1t + m1t;
	auto ref4 = m1 + m1t;

	std::cout << std::boolalpha;

	std::cout << "\nfull- full\n";
	std::cout << (ref1 == m1 + m1) << "\n";
	std::cout << (ref2 == m1t + m1) << "\n";
	std::cout << (ref3 == m1t + m1t) << "\n";
	std::cout << (ref4 == m1 + m1t) << "\n";

	std::cout << "\nfull- full triangle\n";
	std::cout << (ref1 == m1 + m2) << "\n";
	std::cout << (ref2 == m1t + m2) << "\n";
	std::cout << (ref3 == m1t + m2t) << "\n";
	std::cout << (ref4 == m1 + m2t) << "\n";

	std::cout << "\nfull- packed triangle\n";
	std::cout << (ref1 == m1 + m3) << "\n";
	std::cout << (ref2 == m1t + m3) << "\n";
	std::cout << (ref3 == m1t + m3t) << "\n";
	std::cout << (ref4 == m1 + m3t) << "\n";


	std::cout << "\nfull triangle- full\n";
	std::cout << (ref1 == m2 + m1) << "\n";
	std::cout << (ref2 == m2t + m1) << "\n";
	std::cout << (ref3 == m2t + m1t) << "\n";
	std::cout << (ref4 == m2 + m1t) << "\n";

	std::cout << "\nfull triangle - full triangle\n";
	std::cout << (ref1 == m2 + m2) << "\n";
	std::cout << (ref2 == m2t + m2) << "\n";
	std::cout << (ref3 == m2t + m2t) << "\n";
	std::cout << (ref4 == m2 + m2t) << "\n";

	std::cout << "\nfull triangle- packed triangle\n";
	std::cout << (ref1 == m2 + m3) << "\n";
	std::cout << (ref2 == m2t + m3) << "\n";
	std::cout << (ref3 == m2t + m3t) << "\n";
	std::cout << (ref4 == m2 + m3t) << "\n";


	std::cout << "\npack triangle- full\n";
	std::cout << (ref1 == m3 + m1) << "\n";
	std::cout << (ref2 == m3t + m1) << "\n";
	std::cout << (ref3 == m3t + m1t) << "\n";
	std::cout << (ref4 == m3 + m1t) << "\n";

	std::cout << "\npack triangle - full triangle\n";
	std::cout << (ref1 == m3 + m2) << "\n";
	std::cout << (ref2 == m3t + m2) << "\n";
	std::cout << (ref3 == m3t + m2t) << "\n";
	std::cout << (ref4 == m3 + m2t) << "\n";

	std::cout << "\npack triangle- packed triangle\n";
	std::cout << (ref1 == m3 + m3) << "\n";
	std::cout << (ref2 == m3t + m3) << "\n";
	std::cout << (ref3 == m3t + m3t) << "\n";
	std::cout << (ref4 == m3 + m3t) << "\n";
}

void MatrixTest5(void) {
	std::cout << "test operator -= \n\n";
	std::vector<double> v = { 1,0,3,4 };

	std::vector<double> v2 = { 1,3,4 };

	const RowMajorMatrix m1(MatrixType::Full, 2, 2, v);
	const auto m1t = Math::transpose(m1);

	const RowMajorMatrix m2(MatrixType::FullLowerTriangle, 2, 2, v);
	const auto m2t = Math::transpose(m2);

	const RowMajorMatrix m3(MatrixType::PackedLowerTriangle, 2, 2, v2);
	const auto m3t = Math::transpose(m3);

	auto ref1 = m1 - m1;
	auto ref2 = m1t - m1;
	auto ref3 = m1t - m1t;
	auto ref4 = m1 - m1t;

	std::cout << std::boolalpha;
	std::cout << "\nfull- full\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 -= m1)) << "\n";
		std::cout << (ref2 == (f2 -= m1)) << "\n";
		std::cout << (ref3 == (f3 -= m1t)) << "\n";
		std::cout << (ref4 == (f4 -= m1t)) << "\n";
	}
	std::cout << "\nfull - full triangle\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 -= m2)) << "\n";
		std::cout << (ref2 == (f2 -= m2)) << "\n";
		std::cout << (ref3 == (f3 -= m2t)) << "\n";
		std::cout << (ref4 == (f4 -= m2t)) << "\n";
	}
	std::cout << "\nfull - pack triangle\n"; {
		auto f1 = m1;
		auto f2 = m1t;
		auto f3 = m1t;
		auto f4 = m1;

		std::cout << (ref1 == (f1 -= m3)) << "\n";
		std::cout << (ref2 == (f2 -= m3)) << "\n";
		std::cout << (ref3 == (f3 -= m3t)) << "\n";
		std::cout << (ref4 == (f4 -= m3t)) << "\n";
	}


	std::cout << "\nfull triangle - full\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 -= m1)) << "\n";
		std::cout << (ref2 == (f2 -= m1)) << "\n";
		std::cout << (ref3 == (f3 -= m1t)) << "\n";
		std::cout << (ref4 == (f4 -= m1t)) << "\n";
	}
	std::cout << "\nfull triangle - full triangle\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 -= m2)) << "\n";
		std::cout << (ref2 == (f2 -= m2)) << "\n";
		std::cout << (ref3 == (f3 -= m2t)) << "\n";
		std::cout << (ref4 == (f4 -= m2t)) << "\n";
	}
	std::cout << "\nfull triangle - pack triangle\n"; {
		auto f1 = m2;
		auto f2 = m2t;
		auto f3 = m2t;
		auto f4 = m2;

		std::cout << (ref1 == (f1 -= m3)) << "\n";
		std::cout << (ref2 == (f2 -= m3)) << "\n";
		std::cout << (ref3 == (f3 -= m3t)) << "\n";
		std::cout << (ref4 == (f4 -= m3t)) << "\n";
	}


	std::cout << "\npack triangle - full\n"; {
		auto f1 = m3;
		auto f2 = m3t;
		auto f3 = m3t;
		auto f4 = m3;

		std::cout << (ref1 == (f1 -= m1)) << "\n";
		std::cout << (ref2 == (f2 -= m1)) << "\n";
		std::cout << (ref3 == (f3 -= m1t)) << "\n";
		std::cout << (ref4 == (f4 -= m1t)) << "\n";
	}
	std::cout << "\npack triangle - full triangle\n"; {
		auto f1 = m3;
		auto f2 = m3t;
		auto f3 = m3t;
		auto f4 = m3;

		std::cout << (ref1 == (f1 -= m2)) << "\n";
		std::cout << (ref2 == (f2 -= m2)) << "\n";
		std::cout << (ref3 == (f3 -= m2t)) << "\n";
		std::cout << (ref4 == (f4 -= m2t)) << "\n";
	}
	std::cout << "\npack triangle - pack triangle\n"; {
		auto f1 = m3;
		auto f2 = m3t;
		auto f3 = m3t;
		auto f4 = m3;

		std::cout << (ref1 == (f1 -= m3)) << "\n";
		std::cout << (ref2 == (f2 -= m3)) << "\n";
		std::cout << (ref3 == (f3 -= m3t)) << "\n";
		std::cout << (ref4 == (f4 -= m3t)) << "\n";
	}
}

void MatrixTest6(void){
	std::cout << "test operator - \n\n";
	std::vector<double> v = { 1,0,3,4 };

	std::vector<double> v2 = { 1,3,4 };

	const RowMajorMatrix m1(MatrixType::Full, 2, 2, v);
	const auto m1t = Math::transpose(m1);

	const RowMajorMatrix m2(MatrixType::FullLowerTriangle, 2, 2, v);
	const auto m2t = Math::transpose(m2);

	const RowMajorMatrix m3(MatrixType::PackedLowerTriangle, 2, 2, v2);
	const auto m3t = Math::transpose(m3);

	auto ref1 = m1 - m1;
	auto ref2 = m1t - m1;
	auto ref3 = m1t - m1t;
	auto ref4 = m1 - m1t;

	std::cout << std::boolalpha;

	std::cout << "\nfull- full\n";
	std::cout << (ref1 == m1 - m1) << "\n";
	std::cout << (ref2 == m1t - m1) << "\n";
	std::cout << (ref3 == m1t - m1t) << "\n";
	std::cout << (ref4 == m1 - m1t) << "\n";

	std::cout << "\nfull- full triangle\n";
	std::cout << (ref1 == m1 - m2) << "\n";
	std::cout << (ref2 == m1t - m2) << "\n";
	std::cout << (ref3 == m1t - m2t) << "\n";
	std::cout << (ref4 == m1 - m2t) << "\n";

	std::cout << "\nfull- packed triangle\n";
	std::cout << (ref1 == m1 - m3) << "\n";
	std::cout << (ref2 == m1t - m3) << "\n";
	std::cout << (ref3 == m1t - m3t) << "\n";
	std::cout << (ref4 == m1 - m3t) << "\n";


	std::cout << "\nfull triangle- full\n";
	std::cout << (ref1 == m2 - m1) << "\n";
	std::cout << (ref2 == m2t - m1) << "\n";
	std::cout << (ref3 == m2t - m1t) << "\n";
	std::cout << (ref4 == m2 - m1t) << "\n";

	std::cout << "\nfull triangle - full triangle\n";
	std::cout << (ref1 == m2 - m2) << "\n";
	std::cout << (ref2 == m2t - m2) << "\n";
	std::cout << (ref3 == m2t - m2t) << "\n";
	std::cout << (ref4 == m2 - m2t) << "\n";

	std::cout << "\nfull triangle- packed triangle\n";
	std::cout << (ref1 == m2 - m3) << "\n";
	std::cout << (ref2 == m2t - m3) << "\n";
	std::cout << (ref3 == m2t - m3t) << "\n";
	std::cout << (ref4 == m2 - m3t) << "\n";


	std::cout << "\npack triangle- full\n";
	std::cout << (ref1 == m3 - m1) << "\n";
	std::cout << (ref2 == m3t - m1) << "\n";
	std::cout << (ref3 == m3t - m1t) << "\n";
	std::cout << (ref4 == m3 - m1t) << "\n";

	std::cout << "\npack triangle - full triangle\n";
	std::cout << (ref1 == m3 - m2) << "\n";
	std::cout << (ref2 == m3t - m2) << "\n";
	std::cout << (ref3 == m3t - m2t) << "\n";
	std::cout << (ref4 == m3 - m2t) << "\n";

	std::cout << "\npack triangle- packed triangle\n";
	std::cout << (ref1 == m3 - m3) << "\n";
	std::cout << (ref2 == m3t - m3) << "\n";
	std::cout << (ref3 == m3t - m3t) << "\n";
	std::cout << (ref4 == m3 - m3t) << "\n";
}

void MatrixTest7(void) {
	std::cout << " test determinant \n\n";
	std::vector<double> v = { 1,0,3,4 };
	std::vector<double> v2 = { 1,3,4 };

	const RowMajorMatrix m1(MatrixType::Full, 2, 2, v);
	const auto m1t = Math::transpose(m1);

	const RowMajorMatrix m2(MatrixType::FullLowerTriangle, 2, 2, v);
	const auto m2t = Math::transpose(m2);

	const RowMajorMatrix m3(MatrixType::PackedLowerTriangle, 2, 2, v2);
	const auto m3t = Math::transpose(m3);

	auto ref1 = Math::determinant(m1);
	auto ref2 = Math::determinant(m1t);

	std::cout << std::boolalpha;
	std::cout << (ref1 == Math::determinant(m1)) << "\n";
	std::cout << (ref1 == Math::determinant(m2)) << "\n";
	std::cout << (ref1 == Math::determinant(m3)) << "\n";
	std::cout << (ref2 == Math::determinant(m1t)) << "\n";
	std::cout << (ref2 == Math::determinant(m2t)) << "\n";
	std::cout << (ref2 == Math::determinant(m3t)) << "\n";
}

void MatrixTest8(void) {
	std::cout << " test inverse \n\n";
	std::vector<double> v = { 1,0,3,4 };
	std::vector<double> v2 = { 1,3,4 };

	const RowMajorMatrix m1(MatrixType::Full, 2, 2, v);
	const auto m1t = Math::transpose(m1);

	const RowMajorMatrix m2(MatrixType::FullLowerTriangle, 2, 2, v);
	const auto m2t = Math::transpose(m2);

	const RowMajorMatrix m3(MatrixType::PackedLowerTriangle, 2, 2, v2);
	const auto m3t = Math::transpose(m3);

	auto ref1 = Math::inverse(m1);
	auto ref2 = Math::inverse(m1t);

	std::cout << std::boolalpha;
	std::cout << (ref1 == Math::inverse(m1)) << "\n";
	std::cout << (ref1 == Math::inverse(m2)) << "\n";
	std::cout << (ref1 == Math::inverse(m3)) << "\n";	//compare Full and Pack Isuue
	std::cout << (ref2 == Math::inverse(m1t)) << "\n";
	std::cout << (ref2 == Math::inverse(m2t)) << "\n";
	std::cout << (ref2 == Math::inverse(m3t)) << "\n";	//compare Full and Pack Isuue
}