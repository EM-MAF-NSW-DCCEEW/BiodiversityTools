/*
FileUtils.h - Various file management functions
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

#pragma once
#ifndef FILE_UTILS_H
#define FILE_UTILS_H
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#ifdef _MSC_VER
#if _MSC_VER < 1910
#include <boost\filesystem.hpp>
namespace fileSys = boost::filesystem;
#elif _MSC_VER < 1920
#include <filesystem>
namespace fileSys = std::experimental::filesystem;
#else
#include <filesystem>
namespace fileSys = std::filesystem;
#endif
#else
#include <filesystem>
namespace fileSys = std::filesystem;
#endif

const std::string pathSep{ fileSys::path::preferred_separator };

//////////////////////////////////////////
//Legacy file utils using char arrays
//////////////////////////////////////////
//Copy a file using binary file streams
bool CopyFS(const char *srcFN, const char *dstFN, bool overwrite = false);
//Replace a files extension with a new extenstion
bool ReplaceFileExt(const char *inFN, char *outFN, const char *newExt);
//Remove a file's extention
bool RemoveFileExt(const char* inFN, char *outFN);
//Append appStr to file's name prior to its extension 
bool AppendToFileName(const char* inFN, char *outFN, const char *appStr);
//Get the fodler a file is in
bool GetFilesFolder(const char* inFN, char *outFolder);
//Returns whether or not a file exists and can be opened as an input stream
bool FileExists(const char *FN);

//////////////////////////////////////////
//New file utils using strings
//////////////////////////////////////////
//Copy a file using binary file streams checking for overwrite
bool CopyFS(const std::string & srcFN, const std::string dstFN, bool overwrite = false);
//Get a file's extension after last .
bool GetFileExt(const std::string & inFN, std::string & fileExt);
std::string GetFileExt(const std::string & inFN);
//Replace a file's extension with a new extenstion
bool ReplaceFileExt(const std::string & inFN, std::string & outFN, const std::string & newExt);
std::string ReplaceFileExt(const std::string & inFN, const std::string & newExt);
//Remove a file's extention
bool RemoveFileExt(const std::string & inFN, std::string & outFN);
std::string RemoveFileExt(const std::string & inFN);
//Append appStr to file's name prior to its extension if present 
bool AppendToFileName(const std::string & inFN, std::string & outFN, const std::string appStr);
//Get the folder that inFN is in
bool GetFilesFolder(const std::string & inFN, std::string & outFolder);
//Returns whether or not a file exists and can be opened as an input stream
bool FileExists(const std::string & FN);
#endif


// Template functions to write and read a 2D vector of T to a file
//TODO declarations of explicit specializations of function templates

template<typename T>
bool writeVectorToFile(const std::string& filename, const std::vector<T> & data) {
	std::ofstream outFile(filename, std::ios::binary);
	if (!outFile) return false;

	size_t dataSize = data.size();
	outFile.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));
	outFile.write(reinterpret_cast<const char*>(data.data()), dataSize * sizeof(T));

	outFile.close();
	return true;
}

template<typename T>
bool readVectorFromFile(const std::string& filename, std::vector<T> & data) {
	std::ifstream inFile(filename, std::ios::binary);
	if (!inFile) return false;

	size_t dataSize;
	inFile.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));
	data.resize(dataSize);
	inFile.read(reinterpret_cast<char*>(data.data()), dataSize * sizeof(T));

	inFile.close();
	return true;
}

template<typename T>
bool write2dVectorToFile(const std::string& filename, const std::vector<std::vector<T>> & data) {
	std::ofstream outFile(filename, std::ios::binary);
	if (!outFile) return false;

	// Write the size of the outer vector
	size_t outerSize = data.size();
	outFile.write(reinterpret_cast<const char*>(&outerSize), sizeof(outerSize));

	// Write each inner vector
	for (const auto& innerVec : data) {
		size_t innerSize = innerVec.size();
		outFile.write(reinterpret_cast<const char*>(&innerSize), sizeof(innerSize));
		outFile.write(reinterpret_cast<const char*>(innerVec.data()), innerSize * sizeof(T));
	}

	outFile.close();
	return true;
}

template<typename T>
bool read2dVectorFromFile(const std::string& filename, std::vector<std::vector<T>> & data) {
	std::ifstream inFile(filename, std::ios::binary);
	if (!inFile) return false;

	// Read the size of the outer vector
	size_t outerSize;
	inFile.read(reinterpret_cast<char*>(&outerSize), sizeof(outerSize));
	data.resize(outerSize);

	// Read each inner vector
	for (auto& innerVec : data) {
		size_t innerSize;
		inFile.read(reinterpret_cast<char*>(&innerSize), sizeof(innerSize));
		innerVec.resize(innerSize);
		inFile.read(reinterpret_cast<char*>(innerVec.data()), innerSize * sizeof(T));
	}

	inFile.close();
	return true;
}
