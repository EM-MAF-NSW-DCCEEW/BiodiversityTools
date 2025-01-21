/*
FileUtils.cpp - Various file management functions
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


#include <fstream> 
#include <string>
#include <cstring>
#include "FileUtils.h"
#include <iostream>

#pragma warning(disable : 4996)

//////////////////////////////////////////
//Legacy file utils using char arrays
//////////////////////////////////////////

//Copy a file using binary file streams checking for overwrite
bool CopyFS(const char *srcFN, const char *dstFN, bool overwrite)
{
	try { return CopyFS(std::string(srcFN), std::string(dstFN), overwrite); }
	catch (...) { return false; }
}

//Replace a file's extension with a new extenstion
bool ReplaceFileExt(const char *inFN, char *outFN, const char *newExt)
{
	try {
		std::string in(inFN), out;
		size_t lastDot = in.find_last_of('.');
		size_t lastSep = in.find_last_of("\\/");
		if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
			out = in.append(".").append(newExt);
		else
			out = in.substr(0, lastDot + 1).append(newExt);
		strcpy(outFN, out.c_str());
		return true;
	}
	catch (...) { return false; }
}

//Remove a file's extention
bool RemoveFileExt(const char* inFN, char *outFN)
{
	try {
		std::string in(inFN), out;
		size_t lastDot = in.find_last_of('.');
		size_t lastSep = in.find_last_of("\\/");
		if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
			out = std::string(inFN);
		else
			out = in.substr(0, lastDot);
		strcpy(outFN, out.c_str());
		return true;
	}
	catch (...) { return false; }
}

//Append appStr to file's name prior to its extension if present 
bool AppendToFileName(const char* inFN, char *outFN, const char *appStr)
{
	try {
		std::string in(inFN), out;
		size_t lastDot = in.find_last_of('.');
		size_t lastSep = in.find_last_of("\\/");
		if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
			out = in.append(appStr);
		else
			out = in.substr(0, lastDot).append(appStr).append(in.substr(lastDot));
		strcpy(outFN, out.c_str());
		return true;
	}
	catch (...) { return false; }
}

//Get the folder that inFN is in
bool GetFilesFolder(const char* inFN, char *outFolder)
{
	try {
		std::string in(inFN), out;
		size_t lastSep = in.find_last_of("\\/");
		if (lastSep == 0 || lastSep == std::string::npos) return false;
		out = in.substr(0, lastSep);
		strcpy(outFolder, out.c_str());
		return true;
	}
	catch (...) { return false; }
}

//Returns whether or not a file exists and can be opened as an input stream
bool FileExists(const char *FN)
{
	try { return FileExists(std::string(FN)); }
	catch (...) { return false; }
}

//////////////////////////////////////////
//New file utils using strings
//////////////////////////////////////////

//Copy a file using binary file streams checking for overwrite
bool CopyFS(const std::string & srcFN, const std::string dstFN, bool overwrite)
{
	try {
		if (!overwrite && FileExists(dstFN)) return false;
		std::ifstream ifs(srcFN, std::ios::binary);
		std::ofstream ofs(dstFN, std::ios::binary);
		if (!ifs.is_open() || !ofs.is_open()) return false;
		ofs << ifs.rdbuf();
		ofs.flush();
		//ifs.close();
		//ofs.close();
		return true;
	}
	catch (...) { return false; }
}

bool GetFileExt(const std::string & inFN, std::string & fileExt)
{
	try {

		size_t lastDot = inFN.find_last_of('.');
		size_t lastSep = inFN.find_last_of("\\/");
		if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
			return false;
		else
			fileExt = inFN.substr(lastDot + 1);
		return true;
	}
	catch (...) { return false; }

}
std::string GetFileExt(const std::string & inFN)
{
	std::string fileExt("");
	GetFileExt(inFN, fileExt);
	return fileExt;
}


//Replace a file's extension with a new extenstion
bool ReplaceFileExt(const std::string & inFN, std::string & outFN, const std::string & newExt)
{
	try {
		size_t lastDot = inFN.find_last_of('.');
		size_t lastSep = inFN.find_last_of("\\/");
		if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
			outFN = inFN + "." + newExt;
		else
			outFN = inFN.substr(0, lastDot + 1) + newExt;
		return true;
	}
	catch (...) { return false; }
}
std::string ReplaceFileExt(const std::string & inFN, const std::string & newExt)
{
	std::string outFN("");
	ReplaceFileExt(inFN, outFN, newExt);
	return outFN;
}

//Remove a file's extention
bool RemoveFileExt(const std::string & inFN, std::string & outFN)
{
	try {
		size_t lastDot = inFN.find_last_of('.');
		size_t lastSep = inFN.find_last_of("\\/");
		if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
			outFN = inFN;
		else
			outFN = inFN.substr(0, lastDot);
		return true;
	}
	catch (...) { return false; }
}
std::string RemoveFileExt(const std::string & inFN)
{
	std::string outFN("");
	RemoveFileExt(inFN, outFN);
	return outFN;
}


//Append appStr to file's name prior to its extension if present 
bool AppendToFileName(const std::string & inFN, std::string & outFN, const std::string appStr)
{
	try {
		size_t lastDot = inFN.find_last_of('.');
		size_t lastSep = inFN.find_last_of("\\/");
		if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
			outFN = inFN + appStr;
		else
			outFN = inFN.substr(0, lastDot) + appStr + inFN.substr(lastDot);
		return true;
	}
	catch (...) { return false; }
}

//Get the folder that inFN is in
bool GetFilesFolder(const std::string & inFN, std::string & outFolder)
{
	try {
		size_t lastSep = inFN.find_last_of("\\/");
		if (lastSep == 0 || lastSep == std::string::npos) return false;
		outFolder = inFN.substr(0, lastSep);
		return true;
	}
	catch (...) { return false; }
}

//Returns whether or not a file exists and can be opened as an input stream
bool FileExists(const std::string & FN)
{
	try {
		if (FN.empty()) return false;
		std::ifstream fs(FN);
		return fs.is_open() && fs.good();
	}
	catch (...) { return false; }
}

//Get the parent folder of inPath
bool GetParent(const std::string & inPath, std::string & parentPath)
{
	try {
		size_t lastSep = inPath.find_last_of("\\/", inPath.length() - 1);
		if (lastSep == 0 || lastSep == std::string::npos) return false;
		parentPath = inPath.substr(0, lastSep);
		return true;
	}
	catch (...) { return false; }
}
