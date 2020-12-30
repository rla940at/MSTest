#pragma once

#include <iostream>
#include <string>

class FatalError
{
public:
	static void error(std::string file_name, const std::string& function_name, const int num_line, const std::string& comment = "") {
		file_name.erase(file_name.begin(), file_name.begin() + file_name.rfind("\\") + 1);

		std::cout << "\n\n============================ABNORMAL TERMINATION============================\n\n";
		std::cout << "Notice\t\t: " << comment << "\n";
		std::cout << "File\t\t: " << file_name << "\n";
		std::cout << "Function\t: " << function_name << "\n";
		std::cout << "Line\t\t: " << num_line << "\n";
		std::cout << "\n\n============================ABNORMAL TERMINATION============================\n\n";

		std::exit(99);
	};
};

#define FATAL_ERROR(comment)	FatalError::error(__FILE__, __FUNCTION__, __LINE__, std::string() << comment)
#define FATAL_SIZE_ERROR		FATAL_ERROR("Size Error")
#define FATAL_TYPE_ERROR		FATAL_ERROR("Type Error")