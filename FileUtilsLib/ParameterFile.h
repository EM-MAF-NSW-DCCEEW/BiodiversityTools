/*
ParameterFile.h - Parameter file IO function declerations
Copyright(C) 2024 State of New South Wales and Department of Climate Change, Energy, the Environment and Water (DCCEEW)
Author: Jamie Love, Ecological Modelling, Science and Insights Division

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


#pragma once
#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H
#include <map>
#include <vector>

//////////////////////////////////////////
//Legacy Get Parameters using char *
//////////////////////////////////////////
//Get a text parameter in a section of a parameter file
//Can overflow value! Caller needs to ensure value is large enough!
bool GetParam(const char *fileName, const char *section, const char *parameter, char *value);
//Get a integer parameter in a section of a parameter file
bool GetParam(const char *fileName, const char *section, const char *parameter, int &value);
//Get an unsigned integer parameter in a section of a parameter file
bool GetParam(const char *fileName, const char *section, const char *parameter, unsigned int &value);
//Get a float parameter in a section of a parameter file
bool GetParam(const char *fileName, const char *section, const char *parameter, float &value);
//Get a double parameter in a section of a parameter file
bool GetParam(const char *fileName, const char *section, const char *parameter, double &value);
//Get a bool parameter in a section of a parameter file
bool GetParam(const char *fileName, const char *section, const char *parameter, bool &value);

//////////////////////////////////////////
//Legacy Set Parameters using char *
//////////////////////////////////////////
//Set a text parameter in a section of a parameter file
bool SetParam(const char *fileName, const char *section, const char *parameter, const char *value);
//Set a integer parameter in a section of a parameter file
bool SetParam(const char *fileName, const char *section, const char *parameter, const int value);
//Set an unsigned integer parameter in a section of a parameter file
bool SetParam(const char *fileName, const char *section, const char *parameter, const unsigned int value);
//Set a float parameter in a section of a parameter file
bool SetParam(const char *fileName, const char *section, const char *parameter, const float value);
//Set a double parameter in a section of a parameter file
bool SetParam(const char *fileName, const char *section, const char *parameter, const double value);
//Set a bool parameter in a section of a parameter file
bool SetParam(const char *fileName, const char *section, const char *parameter, const bool value);

//////////////////////////////////////////
//Parameter file utility functions
//////////////////////////////////////////
//Check whether a parameter section exists
bool ParamSectionExists(const std::string & fileName, const std::string & section);

//////////////////////////////////////////
//Get Parameters using std::string
//////////////////////////////////////////
//Get a text parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::string & value);
//Get a integer parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, int & value);
//Get a integer parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, unsigned int & value);
//Get a float parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, float & value);
//Get a double parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, double & value);
//Get a bool parameter in a section of a parameter file
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, bool & value);

//Get a comma delimited string parameter in a section of a parameter file as a vector of type T
//template <typename T>
//bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<T> &values);

//Fix cuda_cba.so: undefined reference to `bool GetParam<double>
#include <string>
#include <sstream>

template<typename T>
bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<T> &values)
{
	try {
		T value;
		std::string str, token;
		if (!GetParam(fileName, section, parameter, str)) return false;
		std::stringstream ss(str);
		while (std::getline(ss, token, ',')) {
			token.erase(0, token.find_first_not_of(" \t\n\r\f\v"));
			token.erase(token.find_last_not_of(" \t\n\r\f\v") + 1);
			if (std::stringstream(token) >> std::boolalpha >> value) values.push_back(value);
		}
		return true;
	}
	catch (...) { return false; }
}

// extern template bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<std::string> &values);
// extern template bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<int> &values);
// extern template bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<unsigned int> &values);
// extern template bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<float> &values);
// extern template bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<double> &values);
// extern template bool GetParam(const std::string & fileName, const std::string & section, const std::string & parameter, std::vector<bool> &values);

//////////////////////////////////////////
//Set Parameters using std::string
//////////////////////////////////////////
//Set a text parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const std::string & value);
//Set a char array parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const char * value);
//Set a integer parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const int value);
//Set an unsigned integer parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const unsigned int value);
//Set a float parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const float value);
//Set a double parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const double value);
//Set a bool parameter in a section of a parameter file
bool SetParam(const std::string & fileName, const std::string & section, const std::string & parameter, const bool value);

//ParameterMap
typedef std::map<std::string, std::map<std::string, std::string>> ParameterMap;
//Read a parameter file into a ParameterMap
bool ReadParamFile(const std::string & fileName, ParameterMap &paramMap);
//Write a ParameterMap to a parameter file
bool WriteParamFile(const std::string & fileName, ParameterMap &paramMap);

#endif