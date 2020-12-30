#pragma once

#include "FatalError.h"

#include <algorithm>	//std::transform
#include <iomanip>		//std::setprecision & std::noshowpoint
#include <sstream>		//std::ostringstream
#include <type_traits>	//std::is_arthmatic
#include <typeinfo>		//typeid
#include <vector>
#include <cctype>

namespace Tool {
	template <typename T>
	size_t find_First_Position(const std::string& str, const T& value_to_find, const size_t start_position = 0);

	template <typename T>
	std::vector<size_t> find_Position_Set(const std::string& str, const T& value_to_find);
}


namespace Editor {
	void erase_back(std::string& str, const size_t length);

	template <typename T1, typename T2>
	void replace(std::string& str, const T1& old_value, const T2& new_value);

	template <typename T1, typename T2>
	void replace(std::string& str, const std::vector<T1>& old_value_set, const T2& new_value);

	template <typename T>
	std::string& remove(std::string& str, const T& value_to_remove);

	template <typename T>
	std::string& remove(std::string& str, const std::vector<T>& value_set_to_remove);

	template <typename T>
	std::string& remove(std::string& str, std::initializer_list<T> remove_initialize_list);

	std::string to_String(const std::string& arg);

	template<typename T>
	std::string to_String(const T& arg);

	template<typename T>
	std::vector<std::string> to_String(const std::vector<T>& arg_set);
}


namespace StringEditor {
	bool is_same_without_Case_Sensitivity(const std::string& str1, const std::string& str2);

	std::vector<std::string> parse(const std::string& str, const size_t position);

	std::vector<std::string> parse(const std::string& str, const char delimiter);

	std::vector<std::string> parse(const std::string& str, const std::vector<char>& delimiters);

	std::string& remove_comment(std::string& str, const std::string& comment);

	template<typename ValueType>
	ValueType toValue(const std::string& str);

	std::string& UpperCase(std::string& str);

	std::string UpperCase(const std::string& str);

	std::string UpperCase(std::string&& str);
}


template <typename T>
std::string& operator<<(std::string& str, const T& value);

template <typename T>
std::string& operator<<(std::string&& str, const T& value);


namespace Tool {
	template <typename T>
	size_t find_First_Position(const std::string& str, const T& value_to_find, const size_t start_position) {
		const auto target = Editor::to_String(value_to_find);
		const auto target_size = target.size();

		if (target_size == 0 || str.size() < target_size)
			return std::string::npos;

		for (auto iter = str.begin() + start_position; iter != str.end() - target_size + 1; ++iter) {
			if (std::string(iter, iter + target_size) == target)
				return (iter - str.begin());
		}

		return std::string::npos;
	};

	template <typename T>
	std::vector<size_t> find_Position_Set(const std::string& str, const T& value_to_find) {
		const auto target = Editor::to_String(value_to_find);
		const auto target_size = target.size();

		if (target_size == 0)
			return std::vector<size_t>();

		std::vector<size_t> position_set;
		for (auto iter = str.begin(); iter != str.end() - target_size + 1; iter++)
		{
			if (std::string(iter, iter + target_size) == target)
				position_set.emplace_back(iter - str.begin());
		}

		return position_set;
	};
}


namespace Editor{
	template <typename T1, typename T2>
	void replace(std::string& str, const T1& old_value, const T2& new_value) {
		if (old_value == new_value)
			return;

		const auto old_string = Editor::to_String(old_value);
		const auto new_string = Editor::to_String(new_value);

		if (old_string.empty())
			return;

		size_t position = 0;
		while (true) {
			position = Tool::find_First_Position(str, old_string, position);

			if (position == std::string::npos)
				return;

			str.replace(position, old_string.size(), new_string);
			position += new_string.size();
		}
	};

	template <typename T1, typename T2>
	void replace(std::string& str, const std::vector<T1>& old_value_set, const T2& new_value) {
		const auto old_string_set = Editor::to_String(old_value_set);
		const auto new_string = Editor::to_String(new_value);

		const auto& reference_old_string = old_string_set.front();
		for (const auto& old_string : old_string_set)
			Editor::replace(str, old_string, reference_old_string);
	};

	template <typename T>
	std::string& remove(std::string& str, const T& value_to_remove) {
		const auto target = Editor::to_String(value_to_remove);

		if (str.empty() || target.empty())
			return str;

		while (true) {
			const auto position = Tool::find_First_Position(str, target);

			if (position == std::string::npos)
				return str;

			str.erase(position, target.size());
		}
	};

	template <typename T>
	std::string& remove(std::string& str, const std::vector<T>& value_set_to_remove) {
		const auto& reference_target = value_set_to_remove.front();
		Editor::replace(str, value_set_to_remove, reference_target);

		return Editor::remove(str, reference_target);
	}

	template <typename T>
	std::string& remove(std::string& str, std::initializer_list<T> remove_initialize_list) {
		std::vector<T> value_set_to_erase = remove_initialize_list;
		return Editor::remove(str, value_set_to_erase);
	}

	template<typename T>
	std::string to_String(const T& arg) {
		std::ostringstream out;
		out << std::setprecision(16) << std::noshowpoint << arg;
		return out.str();
	};

	template<typename T>
	std::vector<std::string> to_String(const std::vector<T>& arg_set) {
		std::vector<std::string> string_set;
		string_set.reserve(arg_set.size());

		for (const auto& arg : arg_set)
			string_set.emplace_back(Editor::to_String(arg));

		return string_set;
	}
}


namespace StringEditor {
	template<typename ValueType>
	ValueType toValue(const std::string& str) {
		ValueType value;

		if constexpr (std::is_same_v<ValueType, int>)								
			value = std::stoi(str);
		else if constexpr (std::is_same_v<ValueType, int*>)							
			*value = std::stoi(str);
		else if constexpr (std::is_same_v<ValueType, double>)						
			value = std::stod(str);
		else if constexpr (std::is_same_v<ValueType, double*>)						
			*value = std::stod(str);
		else if constexpr (std::is_same_v<ValueType, unsigned long long>)			
			value = std::stoull(str);
		else if constexpr (std::is_same_v<ValueType, std::string>)					
			value = str;
		else
			FATAL_TYPE_ERROR;

		return value;
	};
}


template <typename T>
std::string& operator<<(std::string& str, const T& value) {
	return str += Editor::to_String(value);
};

template <typename T>
std::string& operator<<(std::string&& str, const T& value) {
	return str += Editor::to_String(value);
}