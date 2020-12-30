#pragma once

#include "QuadratureRule.h"
#include "Text.h"

#include <map>

enum class FigureType;
namespace Editor {
	std::string to_String(const FigureType figure_type);
}


class FigureBase
{
protected:
	class ReferenceTransfomrationPointSet
	{
	private:
		static std::map<std::pair<FigureType, size_t>, std::vector<MathVector>> key_to_value_;

	public:
		static const std::vector<MathVector>& get(const FigureType figure_type, const size_t transformation_order);

	private:
		static std::vector<MathVector> add(const FigureType figure_type, const size_t transformation_order);

		static size_t support_element_order(const FigureType figure_type);

		static std::vector<MathVector> Text_To_Point_Set(Text& reference_transformation_point_set_text);
	};


	class ReferenceQuadratureRule
	{
	private:
		static std::map<std::pair<FigureType, size_t>, QuadratureRule> key_to_value_;

	public:
		static const QuadratureRule& get(const FigureType figure_type, const size_t transformation_order);

	private:
		static QuadratureRule add(const FigureType figure_type, const size_t transformation_order);

		static size_t calculate_Num_Required_Point(const FigureType figure_type, const size_t required_order);

		static size_t calculate_Required_Order(const FigureType figure_type, const size_t integrand_order);

		static QuadratureRule Text_To_Quadrature_Rule(const Text& quadrature_text);
	};


	class ReferencePostPointSet
	{
	private:
		static std::map<std::pair<FigureType, size_t>, std::vector<MathVector>> key_to_value_;

	public:
		static const std::vector<MathVector>& get(const FigureType figure_type, const size_t post_order);

	private:
		static std::vector<MathVector> add(const FigureType figure_type, const size_t post_order);
	};


	class ReferenceConnectivity
	{
	private:
		static std::map<std::pair<FigureType, size_t>, std::vector<std::vector<double>>> key_to_value_;

	public:
		static const std::vector<std::vector<double>>& get(const FigureType figure_type, const size_t post_order);

	private:
		static std::vector<std::vector<double>> add(const FigureType figure_type, const size_t post_order);
	};
};


enum class FigureType
{
	Point,
	Line,
	Triangle, Quadrilateral,
	NotInList
};