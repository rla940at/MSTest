
#include "../INC/StringEditor.h"

namespace Editor{
	std::string to_String(const std::string& arg) {
		return arg;
	}

	void erase_back(std::string& str, const size_t length) {
		if (str.size() <= length) {
			str.clear();
			return;
		}

		str.erase(str.size() - length, length);
	};
}

namespace StringEditor {
	bool is_same_without_Case_Sensitivity(const std::string& str1, const std::string& str2) {
		if (StringEditor::UpperCase(str1) == StringEditor::UpperCase(str2))
			return true;
		else
			return false;
	};

	std::vector<std::string> parse(const std::string& str, const size_t position) {
		if (str.empty())				
			return std::vector<std::string>();

		std::vector<std::string> parsed_string_set;

		const auto iter1 = str.begin();
		const auto iter2 = iter1 + position;

		parsed_string_set.emplace_back(iter1, iter2);
		parsed_string_set.emplace_back(iter2, str.end());

		return parsed_string_set;
	};

	std::vector<std::string> parse(const std::string& str, const char delimiter) {
		std::vector<std::string> parsed_string_set;

		if (str.empty())	
			return parsed_string_set;

		const auto position_set = Tool::find_Position_Set(str, delimiter);

		auto iter1 = str.begin();
		for (const auto& position : position_set){
			auto iter2 = str.begin() + position;

			if (iter2 == iter1)	{
				iter1 = iter2 + 1;
				continue;
			}

			parsed_string_set.emplace_back(iter1, iter2);

			iter1 = iter2 + 1;
		}

		if (iter1 != str.end())
			parsed_string_set.emplace_back(iter1, str.end());

		return parsed_string_set;
	};

	std::vector<std::string> parse(const std::string& str, const std::vector<char>& delimiter_set) {
		if (str.empty())				
			return std::vector<std::string>();
		if (delimiter_set.empty())		
			return std::vector<std::string>(1, str);

		auto temporary_string = str;

		const auto& reference_delimiter = delimiter_set.front();
		for (const auto& other_delimiter : delimiter_set)
		{
			if (other_delimiter != reference_delimiter)
				Editor::replace(temporary_string, other_delimiter, reference_delimiter);
		}

		return StringEditor::parse(temporary_string, reference_delimiter);
	};

	std::string& remove_comment(std::string& str, const std::string& comment) {
		const auto position = Tool::find_First_Position(str, comment);

		return str.erase(position, std::string::npos);
	}


	std::string& UpperCase(std::string& str) {
		std::transform(str.begin(), str.end(), str.begin(), std::toupper);
		return str;
	};

	std::string UpperCase(const std::string& str) {
		auto result = str;
		return StringEditor::UpperCase(result);
	};

	std::string UpperCase(std::string&& str) {
		auto result = std::move(str);
		return StringEditor::UpperCase(result);
	};
}