#pragma once

#include "StringEditor.h"

#include <fstream>
#include <stdio.h> //std::rename
#include <filesystem>


namespace Tool {
	template <typename T>
	size_t find_First_Position(const std::ifstream& read_file, const T& value_to_find);
}


class Text;
namespace Editor{
	void merge(Text& text1, const Text& text2);

	void merge(Text& text1, Text&& text2);

	template <typename T>
	void remove(Text& text, const T& target);

	void replace(Text& text, const std::string& old_str, const std::string& new_str);
}

namespace FileEditor {
	std::ifstream& move_Line(std::ifstream& file_stream, const size_t num_move_line);

	void rename(const std::string& file_path, const std::string& old_name, const std::string& new_name);

	std::vector<std::string> read_Name(const std::string& path);
}

class Text
{
	using Iter = std::vector<std::string>::iterator;
	using CIter = std::vector<std::string>::const_iterator;

private:
	std::vector<std::string> sentence_set_;

private:
	friend void Editor::merge(Text& text1, const Text& text2);

	friend void Editor::merge(Text& text1, Text&& text2);

	template <typename T>
	friend void Editor::remove(Text& text, const T& target);

public:
	Text(void) = default;
	
	Text(const std::vector<std::string>& other_sentence_set)
		: sentence_set_(other_sentence_set) {};

	Text(std::initializer_list<std::string> list)
		: sentence_set_{ list } {};

	Text(const std::string& read_file);

	Text(const std::string& file_path, const std::streampos& start_position);

	Text(const std::string& file_path, const std::streampos& start_position, const size_t num_line);


	template <typename T>
	Text& operator<<(const T& value) {
		this->sentence_set_.emplace_back(Editor::to_String(value));
		return *this;
	};

	Text& operator<<(std::string&& sentence) {
		this->sentence_set_.emplace_back(std::move(sentence));
		return *this;
	};

	std::string& operator[](const size_t index) {
		return this->sentence_set_[index];
	};

	const std::string& operator[](const size_t index) const {
		return this->sentence_set_[index];
	};


	std::string& back(void) {
		return this->sentence_set_.back();
	};

	const std::string& back(void) const {
		return this->sentence_set_.back();
	};

	Iter begin(void) {
		return this->sentence_set_.begin(); 
	};

	CIter begin(void) const {
		return this->sentence_set_.begin();
	};

	Iter end(void) {
		return this->sentence_set_.end();
	};

	CIter end(void) const {
		return this->sentence_set_.end();
	};

	Iter erase(const Iter& iter) {
		return this->sentence_set_.erase(iter);
	};

	std::string& front(void) {
		return this->sentence_set_.front();
	};

	const std::string& front(void) const {
		return this->sentence_set_.front();
	};

	void pop_back(void) {
		this->sentence_set_.pop_back();
	};

	void reserve(const size_t required_capacity) {
		this->sentence_set_.reserve(required_capacity);
	};

	size_t size(void) const {
		return this->sentence_set_.size();
	};


	void add_Write(const std::string& file_path) const;

	template<typename ValueType>
	std::vector<ValueType> to_Value_Set(void) const {
		std::vector<ValueType> value_set;
		value_set.reserve(this->size());

		for (const auto& sentence : this->sentence_set_)
			value_set.push_back(StringEditor::toValue<ValueType>(sentence));

		return value_set;
	};

	void write(const std::string& file_path) const;
};


namespace Tool{
	template <typename T>
	size_t find_First_Position(const std::ifstream& read_file, const T& value_to_find) {
		const auto target = Editor::to_String(value_to_find);
		const auto target_size = target.size();

		if (target_size == 0)
			return std::string::npos;

		if (!read_file.is_open())
			FATAL_ERROR("Fail to open");

		std::streampos position = 0;
		std::string str;
		while (std::getline(read_file, str))
		{
			if (str.find(target) != std::string::npos)
				return position + Tool::find_First_Position(str, target);

			position = read_file.tellg();
		}

		return std::string::npos;
	}
}

namespace Editor{
	template <typename T>
	void remove(Text& text, const T& target) {
		text.sentence_set_.erase(std::remove(text.begin(), text.end(), Editor::to_String(target)), text.end());
	}
}
