/*
ParameterFile.cpp - Parameter file IO function definitions
Copyright(C) 2024 State of New South Wales and Department of Climate Change, Energy, the Environment and Water (DCCEEW)
NSW DCCEEW Metrics and Forecasting Ecological Modelling Team 
ecological.modelling@environment.nsw.gov.au

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <https://www.gnu.org/licenses/>.
*/


#include <string>
#include <cstring>
#include <locale>
#include <fstream>
#include <sstream>
#include <iostream>
#include <ios>
#include <iomanip>
#include <vector>
#include <map>
#include "ParameterFile.h"

//sprintf depreciated
#pragma warning(disable : 4996)

/////////////////////////////////////////////////////////
//Legacy parameter functions using char arrays
/////////////////////////////////////////////////////////

//Get a text parameter in a section of a parameter file
//std::strcpy(value, is.c_str()) can overflow value! Caller needs to ensure value is large enough!
bool GetParam(const char * fileName, const char * section, const char * parameter, char * value)
{
	try {
		bool anySection = (strcmp("*", section) == 0);
		std::string is;
		std::ifstream fs(fileName);
		if (!fs.is_open()) return false;

		if (!anySection) {
			char ss[128];
			sprintf(ss, "[%s]", section);
			while (std::getline(fs, is))
				if (strcmp(is.c_str(), ss) == 0)
					break;
		}

		while (std::getline(fs, is)) {
			if (!anySection && is.front() == '[' && is.back() == ']') break;
			std::stringstream iss(is);
			std::getline(iss, is, '=');
			if (strcmp(is.c_str(), parameter) == 0) {
				std::getline(iss, is);
				strcpy(value, is.c_str());
				return true;
			}
		}
		return false;
	}
	catch (...) { return false; }
}

//Get a integer parameter in a section of a parameter file
bool GetParam(const char * fileName, const char * section, const char * parameter, int & value)
{
	return GetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Get a integer parameter in a section of a parameter file
bool GetParam(const char * fileName, const char * section, const char * parameter, unsigned int & value)
{
	return GetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Get a float parameter in a section of a parameter file
bool GetParam(const char * fileName, const char * section, const char * parameter, float & value)
{
	return GetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Get a double parameter in a section of a parameter file
bool GetParam(const char * fileName, const char * section, const char * parameter, double & value)
{
	return GetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Get a bool parameter in a section of a parameter file
bool GetParam(const char * fileName, const char * section, const char * parameter, bool & value)
{
	return GetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Set a text parameter in a section of a parameter file
bool SetParam(const char * fileName, const char * section, const char * parameter, const char * value)
{
	return SetParam(std::string(fileName), std::string(section), std::string(parameter), std::string(value));
}

//Set a integer parameter in a section of a parameter file
bool SetParam(const char * fileName, const char * section, const char * parameter, const int value)
{
	return SetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Set an unsigned integer parameter in a section of a parameter file
bool SetParam(const char * fileName, const char * section, const char * parameter, const unsigned int value)
{
	return SetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Set a float parameter in a section of a parameter file
bool SetParam(const char * fileName, const char * section, const char * parameter, const float value)
{
	return SetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Set a double parameter in a section of a parameter file
bool SetParam(const char * fileName, const char * section, const char * parameter, const double value)
{
	return SetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}

//Set a bool parameter in a section of a parameter file
bool SetParam(const char * fileName, const char * section, const char * parameter, const bool value)
{
	return SetParam(std::string(fileName), std::string(section), std::string(parameter), value);
}


/////////////////////////////////////////////////////////
//Parameter file utility functions
/////////////////////////////////////////////////////////

//Check whether a parameter section exists
bool ParamSectionExists(const std::string & fileName, const std::string & section)
{
	try {
		std::ifstream fs(fileName);
		if (!fs.is_open()) return false;
		std::string line;
		if (section.compare("*") != 0) {
			while (std::getline(fs, line))
				if (line.compare('[' + section + ']') == 0)
					return true;
		}
		return false;
	}
	catch (...) { return false; }
}


/////////////////////////////////////////////////////////
//Parameter functions using std::string
/////////////////////////////////////////////////////////

//Get a text parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::string & value)
{
	std::string tmp = value;
	try {
		std::ifstream fs(fileName);
		if (!fs.is_open()) return false;
		std::string line;
		if (section.compare("*") != 0) {
			while (std::getline(fs, line))
				if (line.compare('[' + section + ']') == 0)
					break;
		}
		while (std::getline(fs, line)) {
			if (line.front() == '#') continue;
			if (section.compare("*") != 0 && line.front() == '[' && line.back() == ']') break;
			if (line.find("=") > 0 && 
				line.find("=") < line.length() && 
				parameter.compare(line.substr(0,line.find("="))) == 0) 
			{
				value = line.substr(line.find("=") + 1);
				value.erase(0, value.find_first_not_of(" \t\n\r\f\v"));
				value.erase(value.find_last_not_of(" \t\n\r\f\v") + 1);
				return true;
			}
		}
		return false;
	}
	catch (...) { 
		value = tmp;
		return false;
	}
}

//Get an integer parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, int & value)
{
	int tmp = value;
	try {
		std::string str;
		if (!GetParam(fileName, section, parameter, str)) return false;
		value = std::stoi(str);
		return true;
	}
	catch (...) { 
		value = tmp;
		return false;
	}
}

//Get an unsigned integer parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, unsigned int & value)
{
	unsigned int tmp = value;
	try {
		std::string str;
		if (!GetParam(fileName, section, parameter, str)) return false;
		value = unsigned(std::stoi(str));
		return true;
	}
	catch (...) { 
		value = tmp;
		return false;
	}
}

//Get a float parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, float & value)
{
	float tmp = value;
	try {
		std::string str;
		if (!GetParam(fileName, section, parameter, str)) return false;
		value = std::stof(str);
		return true;
	}
	catch (...) { 
		value = tmp;
		return false;
	}
}

//Get a double parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, double & value)
{
	double tmp = value;
	try {
		std::string str;
		if (!GetParam(fileName, section, parameter, str)) return false;
		std::stringstream ss(str);
		ss >> std::setprecision(13) >> value;
		return true;
	}
	catch (...) { 
		value = tmp;
	return false; 
	}
}

//Get a bool parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, bool & value)
{
	try {
		std::string str;
		if (!GetParam(fileName, section, parameter, str)) return false;
		if (str.compare("true") == 0 || str.compare("TRUE") == 0 || str.compare("t") == 0 || str.compare("T") == 0 || str.compare("1") == 0) value = true; 
		else if (str.compare("false") == 0 || str.compare("FALSE") == 0 || str.compare("f") == 0 || str.compare("F") == 0 || str.compare("0") == 0) value = false; 
		else return false;
		return true;
	}
	catch (...) { return false; }
}

//Get a comma delimited string parameter in a section of a parameter file as a vector of type T

//Fix cuda_cba.so: undefined reference to `bool GetParam<double>
// template<typename T>
// bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<T> &values)
// {
// 	try {
// 		T value;
// 		std::string str, token;
// 		if (!GetParam(fileName, section, parameter, str)) return false;
// 		std::stringstream ss(str);
// 		while (std::getline(ss, token, ',')) {
// 			token.erase(0, token.find_first_not_of(" \t\n\r\f\v"));
// 			token.erase(token.find_last_not_of(" \t\n\r\f\v") + 1);
// 			if (std::stringstream(token) >> std::boolalpha >> value) values.push_back(value);
// 		}
// 		return true;
// 	}
// 	catch (...) { return false; }
// }

//Set a text parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const std::string & value)
{
	try {
		bool inSection = false;
		bool found = false;
		std::string line;
		std::vector<std::string> lines;
		std::fstream fs(fileName, std::ios::in);
		if (fs.is_open()) {
			while (std::getline(fs, line)) {
				if (line.front() == '[' && line.back() == ']') {
					if (!found && inSection) {
						if (!lines.empty() && lines.back().compare("") == 0) {//Added 20190806 swap parameter with previous empty line
							lines.back() = parameter + "=" + value;
							lines.push_back("");
						} else
							lines.push_back(parameter + "=" + value);
						found = true;
					}
					inSection = line.compare("[" + section + "]") == 0;
				}
				else if ((inSection || section.compare("*") == 0) && 
					line.find("=") > 0 && 
					line.find("=") < line.length() && 
					parameter.compare(line.substr(0, line.find("="))) == 0 ) 
				{
					line = parameter + "=" + value;
					found = true;
				}
				lines.push_back(line);
			}
			fs.close();
		}
		if (!found && !inSection && section.compare("*") != 0) {
			if (lines.size() > 0) lines.push_back("");
			lines.push_back("[" + section + "]");
		}
		if (!found) lines.push_back(parameter + "=" + value);

		fs.open(fileName, std::ios::out);
		if (!fs.is_open()) return false;
		for (std::string str : lines) fs << str << std::endl;
		fs.close();
		return true;
	}
	catch (...) { return false; }
}

//Set a char array parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const char * value)
{
	try { return SetParam(fileName, section, parameter, std::string(value)); }
	catch (...) { return false; }
}

//Set a integer parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const int value)
{
	try { return SetParam(fileName, section, parameter, std::to_string(value)); }
	catch (...) { return false; }
}

//Set an unsigned integer parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const unsigned int value)
{
	try { return SetParam(fileName, section, parameter, std::to_string(value)); }
	catch (...) { return false; }
}

//Set a float parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const float value)
{
	try { return SetParam(fileName, section, parameter, std::to_string(value)); }
	catch (...) { return false; }
}

//Set a double parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const double value)
{
	try {
		std::stringstream ss;
		ss << std::setprecision(13) << value;
		return SetParam(fileName, section, parameter, ss.str());
	}
	catch (...) { return false; }
}

//Set a bool parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const bool value)
{
	try { return SetParam(fileName, section, parameter, std::string(value ? "true" : "false")); }
	catch (...) { return false; }
}

/////////////////////////////////////////////////////////
//Parameter map functions
/////////////////////////////////////////////////////////

//Read parameter file into ParameterMap
bool ReadParamFile(const std::string & fileName, ParameterMap &paramMap) {
	try {
		std::string line, section, param, value;
		std::ifstream paramFS(fileName);
		int nComments = 0;
		if (!paramFS.is_open()) return false;

		while (std::getline(paramFS, line)) {
			if (line.empty()) continue;
			if (line.front() == '#') {
				line.erase(line.begin());
				paramMap["Comments"][std::string("Comment ").append(std::to_string(++nComments))] = line;
				continue;
			}
			if (line.front() == '[' && line.back() == ']') {
				section = line;
				section.erase(section.end() - 1);
				section.erase(section.begin());
				continue;
			}
			if (line.find('=') != std::string::npos) {
				std::stringstream lineStream(line);
				std::getline(lineStream, param, '=');
				std::getline(lineStream, value);
				paramMap[section][param] = value;
			}
		}
		return true;
	}
	catch (...) { return false; }
}

//Write ParameterMap to parameter file
bool WriteParamFile(const std::string & fileName, ParameterMap &paramMap) {
	try {
		std::ofstream paramFS(fileName);
		if (!paramFS.is_open()) return false;

		for (auto comment : paramMap["Comments"])
			paramFS << "#" << comment.second << std::endl;

		for (auto section : paramMap) {
			if (section.first == "Comments") continue;
			paramFS << std::endl << "[" << section.first << "]" << std::endl;
			for (auto parameter : section.second)
				paramFS << parameter.first << "=" << parameter.second << std::endl;
		}
		return true;
	}
	catch (...) { return false; }
}