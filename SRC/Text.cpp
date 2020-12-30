#include "../INC/Text.h"

namespace Editor{
	void merge(Text& text1, const Text& text2) {
		text1.sentence_set_.insert(text1.end(), text2.begin(), text2.end());
	}

	void merge(Text& text1, Text&& text2) {
		text1.sentence_set_.insert(text1.end(), std::make_move_iterator(text2.begin()), std::make_move_iterator(text2.end()));
	}

	void replace(Text& text, const std::string& old_str, const std::string& new_str){
		for (std::string& file_str : text)
			Editor::replace(file_str, old_str, new_str);
	}
}


namespace FileEditor {
	std::ifstream& move_Line(std::ifstream& file_stream, const size_t num_move_line) {
		size_t line_count = 0;
		std::string str;
		while (std::getline(file_stream, str)){
			if (++line_count == num_move_line)
				break;
		}

		return file_stream;
	}

	void rename(const std::string& file_path, const std::string& old_name, const std::string& new_name) {
		const std::string old_path = file_path + "\\" + old_name;
		const std::string new_path = file_path + "\\" + new_name;

		int result = std::rename(old_path.data(), new_path.data());

		if (result != 0)
			FATAL_ERROR("Fail renaming : " + old_path + " => " + new_path);
	}

	std::vector<std::string> read_Name(const std::string& path) {
		std::filesystem::path p(path);

		std::vector<std::string> name_set;
		for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator(p)) {
			std::string file_path = entry.path().string();
			Editor::remove(file_path, path + "\\");
			name_set.emplace_back(std::move(file_path));
		}

		return name_set;
	}
}


Text::Text(const std::string& file_path){
	std::ifstream read_file(file_path);
	if (!read_file.is_open())
		FATAL_ERROR("Fail to open" + file_path);

	std::string file_string;
	while (std::getline(read_file, file_string))
		this->sentence_set_.emplace_back(file_string);
}

Text::Text(const std::string& file_path, const std::streampos& start_position){
	std::ifstream read_file(file_path);
	if (!read_file.is_open())
		FATAL_ERROR("Fail to open" + file_path);

	if (!(start_position > 0))
		FATAL_ERROR("start position can't be negative.");

	read_file.seekg(start_position);

	std::string file_string;
	while (std::getline(read_file, file_string))
		this->sentence_set_.emplace_back(file_string);
}

Text::Text(const std::string& file_path, const std::streampos& start_position, const size_t num_line){
	std::ifstream read_file(file_path);
	if (!read_file.is_open())
		FATAL_ERROR("Fail to open" + file_path);
	if (!(start_position > 0))
		FATAL_ERROR("start position can't be negative.");

	read_file.seekg(start_position);

	std::string file_string;

	size_t line_count = 0;
	while (std::getline(read_file, file_string))
	{
		if (line_count++ == num_line)
			break;
		else
			this->sentence_set_.emplace_back(file_string);
	}
}

void Text::add_Write(const std::string& file_path) const{
	std::ofstream outfile(file_path, std::ios::app);

	if (!outfile.is_open())
		FATAL_ERROR("Fail to open" + file_path);

	const auto num_sentence = this->sentence_set_.size();
	for (size_t i = 0; i < num_sentence - 1; ++i)
		outfile << this->sentence_set_[i] << "\n";
	outfile << this->sentence_set_[num_sentence - 1];

	outfile.close();
}

void Text::write(const std::string& file_path) const{
	std::ofstream outfile(file_path);

	if (!outfile.is_open())
		FATAL_ERROR("Fail to open" + file_path);

	const auto num_sentence = this->sentence_set_.size();
	for (size_t i = 0; i < num_sentence - 1; ++i)
		outfile << this->sentence_set_[i] << "\n";
	outfile << this->sentence_set_[num_sentence - 1];

	outfile.close();
}