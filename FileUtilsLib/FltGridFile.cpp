/*
FltGridFile.cpp - Various ESRI format floating point grid (.flt,.hdr) related functions
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
#include <iomanip>
#include <random>
#include <memory>

#include <cfloat>
#include <cstring>
//#ifndef VS_2013
//Added gdal2.2.4 built with vs2013
//Need to resolve for 32bit
#include "gdal_priv.h"
#include "gdalwarper.h"
//#endif
#include "FltGridFile.h"
#include "gstream.h"
#include "MessageUtils.h"
#include "FileUtils.h"

//Disable sprintf is depreciated
#pragma warning(disable : 4996)

//Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//TODO check strings for nullPtr

//Returns filename with .flt file extension
std::string GridFN(const std::string & inFN)
{
	if (inFN == "") return "";
	size_t lastDot = inFN.find_last_of('.');
	size_t lastSep = inFN.find_last_of("\\/");
	if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
		return inFN + ".flt";
	else
		return inFN.substr(0, lastDot) + ".flt";
}

//Returns filename with .hdr file extension
std::string HeaderFN(const std::string & inFN)
{
	if (inFN == "") return "";
	size_t lastDot = inFN.find_last_of('.');
	size_t lastSep = inFN.find_last_of("\\/");
	if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
		return inFN + ".hdr";
	else
		return inFN.substr(0, lastDot) + ".hdr";
}

//Returns filename with .prj file extension
std::string PrjFN(const std::string & inFN)
{
	if (inFN == "") return "";
	size_t lastDot = inFN.find_last_of('.');
	size_t lastSep = inFN.find_last_of("\\/");
	if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
		return inFN + ".prj";
	else
		return inFN.substr(0, lastDot) + ".prj";
}

//Returns filename with .tif file extension
//
std::string TifFN(const std::string & inFN)
{
	if (inFN == "") return "";
	size_t lastDot = inFN.find_last_of('.');
	size_t lastSep = inFN.find_last_of("\\/");
	if (lastDot == std::string::npos || (lastSep != std::string::npos && lastDot < lastSep))
		return inFN + ".tif";
	else
		return inFN.substr(0, lastDot) + ".tif";
}

//Get flt and hdr filename 
bool GetFltHdrFileNames(const char * inFN, char * gridFN, char * hdrFN)
{
	try {
		if (sprintf(gridFN, "%s.flt", RemoveFileExt(inFN).c_str()) < 0 ||
			sprintf(hdrFN, "%s.hdr", RemoveFileExt(inFN).c_str()) < 0) return false;
		return true;
	}
	catch (...) { return false; }
}

//Get flt and hdr filename 
bool GetFltHdrFileNames(const std::string &inFN, std::string &gridFN, std::string &hdrFN)
{
	try {
		gridFN = GridFN(inFN);
		hdrFN = HeaderFN(inFN);
		return true;
	}
	catch (...) { return false;	}
}

//Delete all grid files associated with inFN
bool DeleteGrid(const std::string &inFN)
{
	remove(GridFN(inFN).c_str());
	remove(HeaderFN(inFN).c_str());
	remove(PrjFN(inFN).c_str());
	remove((GridFN(inFN) + ".aux.xml").c_str());

	if (FileExists(GridFN(inFN)) || FileExists(HeaderFN(inFN)) || FileExists(PrjFN(inFN))) return false;
	return true;
}

//Get a string from a flt grid header file
bool GetHeader(const char *hdrFN, const char* keyword, char* value)
{
	try {
		char var[32];
		std::ifstream fs(hdrFN, std::ios::in);
		if (!fs.is_open()) return false;
		while (fs >> var >> value)
			if (strcmp(var, keyword) == 0)
				return true;
		return false;
	}
	catch (...) { return false; }
}

bool GetHeader(const std::string & hdrFN, const std::string & keyword, std::string & value)
{
	try {
		std::string var;
		std::ifstream fs(HeaderFN(hdrFN), std::ios::in);
		if (!fs.is_open()) return false;
		while (fs >> var >> value)
			if (var.compare(keyword) == 0)
				return true;
		return false;
	}
	catch (...) { return false; }
}

//Get a integer from a flt grid header file
bool GetHeader(const char *hdrFN, const char* keyword, int &value)
{
	return GetHeader(std::string(hdrFN), std::string(keyword), value);
}
bool GetHeader(const std::string & hdrFN, const std::string & keyword, int &value)
{
	try {
		std::string str;
		if (!GetHeader(hdrFN, keyword, str))
			return false;
		value = std::stoi(str);
		return true;
	}
	catch (...) { return false; }
}

//Get an unsigned integer from a flt grid header file
bool GetHeader(const char *hdrFN, const char* keyword, unsigned int &value)
{
	return GetHeader(std::string(hdrFN), std::string(keyword), value);
}
bool GetHeader(const std::string & hdrFN, const std::string & keyword, unsigned int &value)
{
	try {
		std::string str;
		if (!GetHeader(hdrFN, keyword, str))
			return false;
		value = unsigned(std::stoi(str));
		return true;
	}
	catch (...) { return false; }
}

//Get a float from a flt grid header file
bool GetHeader(const char *hdrFN, const char* keyword, float &value)
{
	return GetHeader(std::string(hdrFN), std::string(keyword), value);
}
bool GetHeader(const std::string & hdrFN, const std::string & keyword, float &value)
{
	try {
		std::string str;
		if (!GetHeader(hdrFN, keyword, str))
			return false;
		value = std::stof(str);
		return true;
	}
	catch (...) { return false; }
}

//Get a double from a flt grid header file
bool GetHeader(const char *hdrFN, const char* keyword, double &value)
{
	return GetHeader(std::string(hdrFN), std::string(keyword), value);
}
bool GetHeader(const std::string & hdrFN, const std::string & keyword, double &value)
{
	try {
		std::string str;
		if (!GetHeader(hdrFN, keyword, str))
			return false;
		value = std::stod(str);
		return true;
	}
	catch (...) { return false; }
}


//Read a flt grid header file
bool ReadGridHeader(const char *hdrFN, unsigned int &nCols, unsigned int &nRows, double &xll, double &yll, double &cellSize, float &noData) 
{
	return ReadGridHeader(std::string(hdrFN), nCols, nRows, xll, yll, cellSize, noData);
}
bool ReadGridHeader(const std::string & hdrFN, unsigned int &nCols, unsigned int &nRows, double &xll, double &yll, double &cellSize, float &noData)
{
	try {
		return (
			GetHeader(hdrFN, "ncols", nCols) &&
			GetHeader(hdrFN, "nrows", nRows) &&
			GetHeader(hdrFN, "xllcorner", xll) &&
			GetHeader(hdrFN, "yllcorner", yll) &&
			GetHeader(hdrFN, "cellsize", cellSize) &&
			GetHeader(hdrFN, "NODATA_value", noData));
	}
	catch (...) { return false; }
}

//Read a flt grid header file ignoring xll and yll
bool ReadGridHeader(const char *hdrFN, unsigned int &nCols, unsigned int &nRows, double &cellSize, float &noData)
{
	return ReadGridHeader(std::string(hdrFN), nCols, nRows, cellSize, noData);
}
bool ReadGridHeader(const std::string & hdrFN, unsigned int &nCols, unsigned int &nRows, double &cellSize, float &noData)
{
	try {
		return (
			GetHeader(hdrFN, "ncols", nCols) &&
			GetHeader(hdrFN, "nrows", nRows) &&
			GetHeader(hdrFN, "cellsize", cellSize) &&
			GetHeader(hdrFN, "NODATA_value", noData));
	}
	catch (...) { return false; }
}

//Write a flt grid header file
bool WriteGridHeader(const char *hdrFN, unsigned int nCols, unsigned int nRows, double xll, double yll, double cellSize, float noData = -9999.0f) 
{
	return WriteGridHeader(std::string(hdrFN), nCols, nRows, xll, yll, cellSize, noData);
}
bool WriteGridHeader(const std::string & hdrFN, unsigned int nCols, unsigned int nRows, double xll, double yll, double cellSize, float noData = -9999.0f)
{
	try {
		std::ofstream hdrFS;
		hdrFS.open(HeaderFN(hdrFN), std::ios::out);
		if (!hdrFS.is_open()) return false;
		hdrFS << "ncols\t\t" << nCols << std::endl;
		hdrFS << "nrows\t\t" << nRows << std::endl;
		hdrFS << "xllcorner\t\t" << std::setprecision(13) << xll << std::endl;
		hdrFS << "yllcorner\t\t" << std::setprecision(13) << yll << std::endl;
		hdrFS << "cellsize\t\t" << cellSize << std::endl;
		hdrFS << "NODATA_value\t\t" << noData << std::endl;
		hdrFS << "byteorder\t\t" << "LSBFIRST" << std::endl;
		hdrFS.close();
		return true;
	}
	catch (...) { return false; }
}

//Write ascii header to file
bool WriteAsciiHeader(const char *asciiFN, unsigned int nCols, unsigned int nRows, double xll, double yll, double cellSize, float noData = -9999.0f)
{
	return WriteAsciiHeader(std::string(asciiFN), nCols, nRows, xll, yll, cellSize, noData);
}
bool WriteAsciiHeader(const std::string & asciiFN, unsigned int nCols, unsigned int nRows, double xll, double yll, double cellSize, float noData = -9999.0f)
{
	try {
		std::ofstream hdrFS;
		hdrFS.open(asciiFN, std::ios::out);
		if (!hdrFS.is_open()) return false;
		hdrFS << "ncols\t\t" << nCols << std::endl;
		hdrFS << "nrows\t\t" << nRows << std::endl;
		hdrFS << "xllcorner\t\t" << std::setprecision(13) << xll << std::endl;
		hdrFS << "yllcorner\t\t" << std::setprecision(13) << yll << std::endl;
		hdrFS << "cellsize\t\t" << cellSize << std::endl;
		hdrFS << "NODATA_value\t\t" << noData << std::endl;
		hdrFS.close();
		return true;
	}
	catch (...) { return false; }
}

//Compare two flt grid header files
bool CompareGridHeaders(const char *hdrFNA, const char *hdrFNB)
{
	return CompareGridHeaders(std::string(hdrFNA), std::string(hdrFNB));
}
bool CompareGridHeaders(const std::string & hdrFNA, const std::string & hdrFNB)
{
	try {
		unsigned int nColsA, nColsB;
		unsigned int nRowsA, nRowsB;
		double xllA, xllB;
		double yllA, yllB;
		double cellSizeA, cellSizeB;
		float noDataA, noDataB;
		double xyExtTol = 0.0001;

		if (!ReadGridHeader(HeaderFN(hdrFNA), nColsA, nRowsA, xllA, yllA, cellSizeA, noDataA) ||
			!ReadGridHeader(HeaderFN(hdrFNB), nColsB, nRowsB, xllB, yllB, cellSizeB, noDataB))
			return false;

		return (
			nColsA == nColsB &&
			nRowsA == nRowsB &&
			//Fix for REMP Pu->Pi grid where habitat grid extents have high nDP 
			//Added telerance to comparison for when extent a xor b precision has been rounded
			std::abs(xllA - xllB) < xyExtTol &&
			std::abs(yllA - yllB) < xyExtTol &&
			cellSizeA == cellSizeB &&
			noDataA == noDataB);
	}
	catch (...) { return false; }
}

//Compare two flt grid header files without using extents
bool CompareGridHeadersNoExt(const std::string & hdrFNA, const std::string & hdrFNB)
{
	try {
		unsigned int nColsA, nColsB;
		unsigned int nRowsA, nRowsB;
		double cellSizeA, cellSizeB;
		float noDataA, noDataB;

		if (!ReadGridHeader(HeaderFN(hdrFNA), nColsA, nRowsA, cellSizeA, noDataA) ||
			!ReadGridHeader(HeaderFN(hdrFNB), nColsB, nRowsB, cellSizeB, noDataB))
			return false;

		return (
			nColsA == nColsB &&
			nRowsA == nRowsB &&
			cellSizeA == cellSizeB &&
			noDataA == noDataB);
	}
	catch (...) { return false; }
}

//Get the min value from a flt grid
bool GetGridMin(const char *FN, float &min) 
{ return GetGridMin(std::string(FN), min); }
bool GetGridMin(const std::string & FN, float &min)
{
	try {
		int nCols, nRows;
		float noData, imin = FLT_MAX;
		std::ifstream gridFS(GridFN(FN), std::ios::binary);
		if (!gridFS.is_open()) return false;
		if (!GetHeader(HeaderFN(FN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(FN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(FN), "NODATA_value", noData)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			gridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				imin = fminf(dataPtr[c], imin);
			}
		}
		min = imin;
		return true;
	}
	catch (...) { return false; }
}

template<typename T>
bool GetGridMin2(const std::string & FN, T &min)
{
	try {
		//int nCols, nRows;
		T imin = FLT_MAX;
		igstream gridFS(FN);
		if (!gridFS.is_open()) return false;
		T noData = (T) gridFS.noData(); // if (!GetHeader(HeaderFN(FN), "NODATA_value", noData)) return false;
		std::unique_ptr<T[]> dataPtr = std::make_unique<T[]>(gridFS.nCols());
		for (int r = 0; r < gridFS.nRows(); r++) {
			gridFS.read((char *)dataPtr.get(), gridFS.nCols() * sizeof(T));// sizeof(float));
			for (int c = 0; c < gridFS.nCols(); c++) {
				if (dataPtr[c] == noData) continue;
				imin = std::min(dataPtr[c], imin);
			}
		}
		min = imin;
		return true;
	}
	catch (...) { return false; }
}

//Get the max value from a flt grid
bool GetGridMax(const char *FN, float &max)
{ return GetGridMax(std::string(FN), max); }
bool GetGridMax(const std::string & FN, float &max)
{
	try {
		int nCols, nRows;
		float noData, imax = -FLT_MAX;
		std::ifstream gridFS(GridFN(FN), std::ios::binary);
		if (!gridFS.is_open()) return false;
		if (!GetHeader(HeaderFN(FN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(FN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(FN), "NODATA_value", noData)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			gridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				imax = fmaxf(dataPtr[c], imax);
			}
		}
		max = imax;
		return true;
	}
	catch (...) { return false; }
}

//Get the smallest positive power of 10 greater or equal to the highest grid value
bool GetGridMaxPow10(const std::string gridFN, float &maxPow10) {
	try {
		float gridMax;
		if (!GetGridMax(GridFN(gridFN), gridMax)) return false;
		if (gridMax <= 0.0f) return false;
		float limit = powf(10.0f, ceilf(log10(gridMax)));
		if (limit > 0.0f) {
			maxPow10 = limit;
			return true;
		}
		return false;
	}
	catch (...) { return false; }
}

//Get the min and max value from a flt grid
bool GetGridMinMax(const char *FN, float &min, float &max)
{ return GetGridMinMax(std::string(FN), min, max); }
bool GetGridMinMax(const std::string & FN, float &min, float &max)
{
	try {
		int nCols, nRows;
		float noData, imin = FLT_MAX, imax = -FLT_MAX;
		std::ifstream gridFS(GridFN(FN), std::ios::binary);
		if (!gridFS.is_open()) return false;
		if (!GetHeader(HeaderFN(FN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(FN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(FN), "NODATA_value", noData)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			gridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				imin = fmin(dataPtr[c], imin);
				imax = fmaxf(dataPtr[c], imax);
			}
		}
		min = imin;
		max = imax;
		return true;
	}
	catch (...) { return false; }
}

//Get the sum of values from a flt grid
bool GetGridSum(const char *FN, double &sum)
{
	return GetGridSum(std::string(FN), sum);
}
bool GetGridSum(const std::string & FN, double &sum)
{
	try {
		int nCols, nRows;
		float noData;
		double isum = 0.0;
		std::ifstream gridFS(GridFN(FN), std::ios::binary);
		if (!gridFS.is_open()) return false;
		if (!GetHeader(HeaderFN(FN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(FN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(FN), "NODATA_value", noData)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			gridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				if (dataPtr[c] > DBL_MAX - isum) return false;
				isum += double(dataPtr[c]);
			}
		}
		sum = isum;
		return true;
	}
	catch (...) { return false; }
}

//Get the average of values from a flt grid
bool GetGridAvg(const char *FN, double &avg)
{
	return GetGridAvg(std::string(FN), avg);
}
bool GetGridAvg(const std::string & FN, double &avg)
{
	try {
		int nCols, nRows;
		float noData;
		double mean = 0.0, delta = 0.0;
		ull count = 0;

		std::ifstream gridFS(GridFN(FN), std::ios::binary);
		if (!gridFS.is_open()) return false;
		if (!GetHeader(HeaderFN(FN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(FN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(FN), "NODATA_value", noData)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			gridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				if (ULLONG_MAX - count < 1ULL) return false;
				delta = dataPtr[c] - mean;
				mean += delta / ++count;
			}
		}
		avg = mean;
		return true;
	}
	catch (...) { return false; }
}

//Get the number of data cells in a flt grid
bool CountDataCells(const std::string & FN, ull & nDataCells)
{
	try {
		int nCols, nRows;
		float noData;
		ull count = 0;
		
		std::ifstream gridFS(GridFN(FN), std::ios::binary);
		if (!gridFS.is_open()) return false;
		if (!GetHeader(HeaderFN(FN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(FN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(FN), "NODATA_value", noData)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			gridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				if (ULLONG_MAX - count < 1ULL) return false;
				++count;
			}
		}
		nDataCells = count;
		return true;
	}
	catch (...) { return false; }
}

//Rescale flt grid values to between outMin and outMax
bool RescaleGrid(const char *inFN, const char *outFN, float outMin, float outMax)
{ return RescaleGrid(std::string(inFN), std::string(outFN), outMin, outMax); }
bool RescaleGrid(const std::string & inFN, const std::string & outFN, float inMin, float inMax, float outMin, float outMax)
{
	try {
		int nCols, nRows;
		float noData;
		float inRange = inMax - inMin;
		float outRange = outMax - outMin;
		if (inRange == 0.0f) inRange = 1.0f;//output = outMin if inRange == 0.0f
		std::ifstream inGridFS(GridFN(inFN), std::ios::binary);
		std::ofstream outGridFS(GridFN(outFN), std::ios::binary);
		if (!(inGridFS.is_open() && outGridFS.is_open())) return false;
		if (!GetHeader(HeaderFN(inFN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(inFN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(inFN), "NODATA_value", noData)) return false;
		if (!CopyFS(HeaderFN(inFN), HeaderFN(outFN), true)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			inGridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				dataPtr[c] = (dataPtr[c] - inMin) / inRange * outRange + outMin;
			}	
			outGridFS.write((const char *)dataPtr.get(), nCols * sizeof(float));
		}
		return true;
	}
	catch (...) { return false; }
}

bool RescaleGrid(const std::string & inFN, const std::string & outFN, float outMin, float outMax){
	float inMin = FLT_MAX;
	float inMax = -FLT_MAX;
	if (!GetGridMinMax(inFN, inMin, inMax)) return false;
	return RescaleGrid(inFN, outFN, inMin, inMax, outMin, outMax);
}

//Create a permeability flt grid from a habitat grid
bool CreatePrmGrid(const char *habFN, const char *prmFN, float distMin, float distMax)
{ return CreatePrmGrid(std::string(habFN), std::string(prmFN), distMin, distMax); }
bool CreatePrmGrid(const std::string & habFN, const std::string & prmFN, float distMin, float distMax)
{
	float cellSize;
	float inMin = FLT_MAX;
	float inMax = -FLT_MAX;
	if (!GetGridMinMax(habFN, inMin, inMax)) return false;
	if (!GetHeader(HeaderFN(habFN), "cellsize", cellSize)) return false;
	float outMin = expf(-cellSize / distMin);
	float outMax = expf(-cellSize / distMax);
	return RescaleGrid(habFN, prmFN, inMin, inMax, outMin, outMax);
}

//Create a circle flt grid area in ha
bool CreateCircleGrid(const char *outFN, double cellSize, double areaHA, float inVal, float outVal) 
{ return CreateCircleGrid(std::string(outFN), cellSize, areaHA, inVal, outVal); }
bool CreateCircleGrid(const std::string & outFN, double cellSize, double areaHA, float inVal, float outVal)
{
	try {
		//Create the circle
		double rMetres = sqrt(areaHA / M_PI) * 100;
		double rCells = rMetres / (cellSize > 1.0f ? cellSize : cellSize * 100000);
		int dim = int(ceil(rCells)) * 8 - 1;
		int cntr = dim / 2;

		//Open output grid file streams and Write output grid header
		std::ofstream outFS(GridFN(outFN), std::ios::binary);
		if (!outFS.is_open()) return false;
		if (!WriteGridHeader(HeaderFN(outFN), dim, dim, 0.0f, 0.0f, cellSize, -9999.0f)) return false;

		//Write output grid data
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(dim);
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				dataPtr[j] = (pow(abs(i - cntr), 2.0f) + pow(abs(j - cntr), 2.0f) <= pow(rCells, 2.0f)) ? inVal : outVal;
			}
			outFS.write((const char *)(dataPtr.get()), dim * sizeof(float));
		}
		return true;
	}
	catch (...) { return false; }
}

//Create a circle flt grid area in ha
bool CreateCircleGrid(const char *outFN, double cellSize, double areaHA, int dim, float inVal, float outVal) 
{ return CreateCircleGrid(std::string(outFN), cellSize, areaHA, dim, inVal, outVal); }
bool CreateCircleGrid(const std::string & outFN, double cellSize, double areaHA, int dim, float inVal, float outVal)
{
	try {
		//Create the circle
		double rMetres = sqrt(areaHA / M_PI ) * 100;
		double rCells = rMetres / (cellSize > 1.0f ? cellSize : cellSize * 100000);
		int cntr = dim / 2;

		//Open output grid file streams and Write output grid header
		std::ofstream outFS(GridFN(outFN), std::ios::binary);
		if (!outFS.is_open()) return false;
		if (!WriteGridHeader(HeaderFN(outFN), dim, dim, 0.0f, 0.0f, cellSize, -9999.0f)) return false;

		//Write output grid data
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(dim);
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				dataPtr[j] = (pow(abs(i - cntr), 2.0f) + pow(abs(j - cntr), 2.0f) <= pow(rCells, 2.0f)) ? inVal : outVal;
			}
			outFS.write((const char *)dataPtr.get(), dim * sizeof(float));
		}
		return true;
	}
	catch (...) { return false; }
}

//Create a random point grid using a habitat weight grid
//TODO needs comparison with previous Links method and improved output text
bool CreateRandomPointGrid(const char *inFN, const char *outFN, float p, float minVal)
{ return CreateRandomPointGrid(std::string(inFN), std::string(outFN), p, minVal); }
bool CreateRandomPointGrid(const std::string & inFN, const std::string & outFN, float p, float minVal)
{
	try {
		unsigned int nCols, nRows, nCells;
		double cellSize;
		float inVal, outVal, noData;
		int nPoints = 0;
		int nValids = 0;

		//Random number generators
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(minVal, 1.0);
		std::uniform_real_distribution<> dis2(0.0, 1.0);

		//Open grid files and read input header
		std::ifstream inFS(GridFN(inFN), std::ios::in | std::ios::binary);
		std::ofstream outFS(GridFN(inFN), std::ios::out | std::ios::binary);
		if (!inFS.is_open() || !outFS.is_open()) return false;
		if (!ReadGridHeader(HeaderFN(inFN), nCols, nRows, cellSize, noData)) return false;
		nCells = nCols * nRows;

		//p = p * (1 - minVal) + minVal;
		while (inFS.read((char*)&inVal, sizeof(float))) {
			outVal = 0.0f;
			if (inVal >= minVal) nValids++;
			if (inVal == noData) outVal = noData;
			else if (inVal >= dis(gen) && p > dis2(gen)) {
				outVal = 1.0f;
				nPoints++;
			}
			outFS.write((const char*)&outVal, sizeof(float));
		}

		//Display results
		std::cout << "Number of Cells: " << nCells << std::endl;
		std::cout << "Number of Habitat Cells: " << nValids << std::endl;
		std::cout << "Number of Points: " << nPoints << std::endl;
		std::cout << "Ratio: " << float(nPoints) / float(nValids) * 100.0f << std::endl;
		return 0;
	}
	catch (...) { return false; }
}

//Resample a flt grid using GDAL and an output cellsize
bool ResampleGrid(const char *inFN, const char *outFN, double outCellSize, int method) 
{ return ResampleGrid(std::string(inFN), std::string(outFN), outCellSize, method); }
bool ResampleGrid(const std::string & inFN, const std::string & outFN, double outCellSize, int method) 
{
	try {
		unsigned int inCols, inRows, outCols, outRows;
		double inCellSize, inXll, inYll, outXll, outYll, cellFactor;
		float noData;

		//Read grid header info and just copy grid if out cell size is the same
		if (!ReadGridHeader(HeaderFN(inFN), inCols, inRows, inXll, inYll, inCellSize, noData)) return false;
		if (inCellSize == outCellSize) {
			CopyFS(GridFN(inFN), GridFN(outFN), true);
			CopyFS(HeaderFN(inFN), HeaderFN(outFN), true);
			return true;
		}
		//Get output grid dimensions
		cellFactor = inCellSize / outCellSize;
		outCols = int(ceil((double)inCols * cellFactor));
		outRows = int(ceil((double)inRows * cellFactor));

		//Register drivers, get input and create output datasets
		GDALRegister_EHdr();
		GDALDriver *gdalDriver = GetGDALDriverManager()->GetDriverByName("EHdr");
		if (gdalDriver == NULL) return false;
		GDALDataset *inDataset = (GDALDataset *)GDALOpen(GridFN(inFN).c_str(), GA_ReadOnly);
		GDALDataset *outDataset = gdalDriver->Create(GridFN(outFN).c_str(), outCols, outRows, 1, GDT_Float32, NULL);
		if (inDataset == NULL || outDataset == NULL) return false;

		//Set output datasets projection, geo transformation and no data value
		double outGeotransform[6];
		inDataset->GetGeoTransform(outGeotransform);
		outGeotransform[1] /= cellFactor;//Update pixel width
		outGeotransform[5] /= cellFactor;//Update pixel height
		outGeotransform[3] -= ((inCellSize * (double)inRows) - (outCellSize * (double)outRows));//Update top left
		outDataset->SetGeoTransform(outGeotransform);
		outDataset->SetProjection(inDataset->GetProjectionRef());
//#ifndef VS_2013
//		outDataset->GetBands()[0]->SetNoDataValue((double)noData);
//#else
		outDataset->GetRasterBand(1)->SetNoDataValue((double)noData);
		outDataset->GetRasterBand(1)->Fill((double)noData);
//#endif
		//Create warp options
		GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
		psWarpOptions->hSrcDS = (GDALDatasetH)inDataset;
		psWarpOptions->hDstDS = (GDALDatasetH)outDataset;
		psWarpOptions->nBandCount = 1;
		psWarpOptions->panSrcBands = (int *)CPLMalloc(sizeof(int));
		psWarpOptions->panSrcBands[0] = 1;
		psWarpOptions->panDstBands = (int *)CPLMalloc(sizeof(int));
		psWarpOptions->panDstBands[0] = 1;
//#ifndef VS_2013
		psWarpOptions->padfSrcNoDataReal = (double *)CPLMalloc(sizeof(double));
		psWarpOptions->padfSrcNoDataReal[0] = (double)noData;
		psWarpOptions->padfDstNoDataReal = (double *)CPLMalloc(sizeof(double));
		psWarpOptions->padfDstNoDataReal[0] = (double)noData;
#ifdef VS_2013 //Added as GDAL 2.2.4 doesn't like Real nodata with no Imag nodata and warper ignores src/dst band nodata 
		psWarpOptions->padfSrcNoDataImag = (double *)CPLMalloc(sizeof(double));
		psWarpOptions->padfSrcNoDataImag[0] = 0.0;
		psWarpOptions->padfDstNoDataImag = (double *)CPLMalloc(sizeof(double));
		psWarpOptions->padfDstNoDataImag[0] = 0.0;
#endif
		//Added 230919 to resolve 0.0 being assigned to no data cells for large grids 
		psWarpOptions->papszWarpOptions = CSLSetNameValue(psWarpOptions->papszWarpOptions, "INIT_DEST", "NO_DATA");
		psWarpOptions->papszWarpOptions = CSLSetNameValue(psWarpOptions->papszWarpOptions, "NUM_THREADS", "ALL_CPUS");

		psWarpOptions->pfnProgress = msgGDALProgress;
		psWarpOptions->pTransformerArg = GDALCreateGenImgProjTransformer(
			(GDALDatasetH)inDataset,
			inDataset->GetProjectionRef(),
			(GDALDatasetH)outDataset,
			outDataset->GetProjectionRef(),
			FALSE, 0.0, 1);
		psWarpOptions->pfnTransformer = GDALGenImgProjTransform;
		psWarpOptions->eResampleAlg = GDALResampleAlg(method);

		//Resample with warp operation
		GDALWarpOperation operation;
		operation.Initialize(psWarpOptions);
		operation.ChunkAndWarpImage(0, 0, outCols, outRows);

		//Clean up GDAL objects
		GDALDestroyGenImgProjTransformer(psWarpOptions->pTransformerArg);
		GDALDestroyWarpOptions(psWarpOptions);
		GDALClose((GDALDatasetH)inDataset);
		GDALClose((GDALDatasetH)outDataset);

		//Overwrite header with correct ESRI format and remove .aux.xml file
		outXll = outGeotransform[0];
		outYll = outGeotransform[3] + outDataset->GetRasterYSize() * outGeotransform[5];
		if (!WriteGridHeader(HeaderFN(outFN), outCols, outRows, outXll, outYll, outCellSize, noData)) return false;
		remove((GridFN(outFN) + ".aux.xml").c_str());
		return true;
	}
	catch (...) { return false; }
}

//Resample a flt grid using GDAL and template flt grid dimensions
bool ResampleGrid(const char *inFN, const char *outFN, const char *tmpFN, int method)
{
	return ResampleGrid(std::string(inFN), std::string(outFN), std::string(tmpFN), method);
}
bool ResampleGrid(const std::string & inFN, const std::string & outFN, const std::string & tmpFN, int method)
{
	try {
		unsigned int outCols, outRows;
		double outXll, outYll, outCellSize;
		float noData;

		//Read template header info
		if (!ReadGridHeader(HeaderFN(tmpFN), outCols, outRows, outXll, outYll, outCellSize, noData)) return false;

		//Change 220223 Get noData from inFN 
		if (!GetHeader(HeaderFN(inFN), "NODATA_value", noData)) return false;

		//Register drivers, get input, template and create output datasets
		GDALRegister_EHdr();
		GDALDriver *gdalDriver = GetGDALDriverManager()->GetDriverByName("EHdr");
		if (gdalDriver == NULL) return false;
		GDALDataset *inDataset = (GDALDataset *)GDALOpen(GridFN(inFN).c_str(), GA_ReadOnly);
		GDALDataset *tmpDataset = (GDALDataset *)GDALOpen(GridFN(tmpFN).c_str(), GA_ReadOnly);
		//GDALDataset *outDataset = gdalDriver->Create(GridFN(outFN).c_str(), outCols, outRows, 1, GDT_Float32, NULL);

		//Added support for writing to tifs
		//TODO move resample functionality to gstream
		GDALRegister_GTiff();
		GDALDriver *gdalTifDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		if (gdalTifDriver == NULL) return false;
		char ** gdalCreateOptions = nullptr;
		gdalCreateOptions = CSLAddNameValue(gdalCreateOptions, "COMPRESS", "LZW");
		GDALDataset *outDataset;

		std::string fileExt = "";
		GetFileExt(outFN, fileExt);
		if (fileExt.compare("tif") == 0)
			outDataset = gdalTifDriver->Create(outFN.c_str(), outCols, outRows, 1, GDT_Float32, gdalCreateOptions);
		else
			outDataset = gdalDriver->Create(GridFN(outFN).c_str(), outCols, outRows, 1, GDT_Float32, NULL);



		if (inDataset == NULL || tmpDataset == NULL || outDataset == NULL) return false;

		//Set output datasets projection, geo transformation and no data value to match template
		double outGeotransform[6];
		tmpDataset->GetGeoTransform(outGeotransform);
		outDataset->SetGeoTransform(outGeotransform);
		outDataset->SetProjection(tmpDataset->GetProjectionRef());
		//#ifndef VS_2013
		//		outDataset->GetBands()[0]->SetNoDataValue((double)noData);
		//#else
		outDataset->GetRasterBand(1)->SetNoDataValue((double)noData);
		outDataset->GetRasterBand(1)->Fill((double)noData);
		//#endif
		GDALClose((GDALDatasetH)tmpDataset);

		//Create warp options
		GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
		psWarpOptions->hSrcDS = (GDALDatasetH)inDataset;
		psWarpOptions->hDstDS = (GDALDatasetH)outDataset;
		psWarpOptions->nBandCount = 1;
		psWarpOptions->panSrcBands = (int *)CPLMalloc(sizeof(int));
		psWarpOptions->panSrcBands[0] = 1;
		psWarpOptions->panDstBands = (int *)CPLMalloc(sizeof(int));
		psWarpOptions->panDstBands[0] = 1;
		//#ifndef VS_2013
		psWarpOptions->padfSrcNoDataReal = (double *)CPLMalloc(sizeof(double));
		psWarpOptions->padfSrcNoDataReal[0] = (double)noData;
		psWarpOptions->padfDstNoDataReal = (double *)CPLMalloc(sizeof(double));
		psWarpOptions->padfDstNoDataReal[0] = (double)noData;
#ifdef VS_2013
		psWarpOptions->padfSrcNoDataImag = (double *)CPLMalloc(sizeof(double));
		psWarpOptions->padfSrcNoDataImag[0] = 0.0;
		psWarpOptions->padfDstNoDataImag = (double *)CPLMalloc(sizeof(double));
		psWarpOptions->padfDstNoDataImag[0] = 0.0;
#endif
		//Added 230919 to resolve 0.0 being assigned to no data cells for large grids 
		psWarpOptions->papszWarpOptions = CSLSetNameValue(psWarpOptions->papszWarpOptions, "INIT_DEST", "NO_DATA");
		psWarpOptions->papszWarpOptions = CSLSetNameValue(psWarpOptions->papszWarpOptions, "NUM_THREADS", "ALL_CPUS");

		psWarpOptions->pfnProgress = msgGDALProgress;
		psWarpOptions->pTransformerArg = GDALCreateGenImgProjTransformer(
			(GDALDatasetH)inDataset,
			inDataset->GetProjectionRef(),
			(GDALDatasetH)outDataset,
			outDataset->GetProjectionRef(),
			FALSE, 0.0, 1);
		psWarpOptions->pfnTransformer = GDALGenImgProjTransform;
		psWarpOptions->eResampleAlg = GDALResampleAlg(method);

		//Resample with warp operation
		GDALWarpOperation operation;
		operation.Initialize(psWarpOptions);
		operation.ChunkAndWarpImage(0, 0, outCols, outRows);

		//Clean up GDAL objects
		GDALDestroyGenImgProjTransformer(psWarpOptions->pTransformerArg);
		GDALDestroyWarpOptions(psWarpOptions);
		GDALClose((GDALDatasetH)inDataset);
		GDALClose((GDALDatasetH)outDataset);

		//Overwrite header with correct ESRI format and remove .aux.xml file
		//Only for flt grids
		if (fileExt.compare("tif") != 0) {
			if (!WriteGridHeader(HeaderFN(outFN), outCols, outRows, outXll, outYll, outCellSize, noData)) return false;
			remove((GridFN(outFN) + ".aux.xml").c_str());
		}
		return true;
	}
	catch (...) { return false; }
}



//Applies the provided function to each non-nodata value in xFN and writes results to yFN 
bool ApplyFuncToGrid(const char *xFN, const char *yFN, float(*func)(float)) 
{ return ApplyFuncToGrid(std::string(xFN), std::string(yFN), func); }
bool ApplyFuncToGrid(const std::string & xFN, const std::string & yFN, float(*func)(float)) 
{
	try {
		int nCols, nRows;
		float noData;
		std::ifstream xGridFS(GridFN(xFN), std::ios::binary);
		std::ofstream yGridFS(GridFN(yFN), std::ios::binary);
		if (!(xGridFS.is_open() && yGridFS.is_open())) return false;
		if (!GetHeader(HeaderFN(xFN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(xFN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(xFN), "NODATA_value", noData)) return false;
		if (!CopyFS(HeaderFN(xFN), HeaderFN(yFN), true)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			xGridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				dataPtr[c] = func(dataPtr[c]);
			}
			yGridFS.write((const char *)dataPtr.get(), nCols * sizeof(float));
		}
		return true;
	}
	catch (...) { return false; }
}

//Applies the provided lambda to each non-nodata value in xFN and writes results to yFN 
bool ApplyLambdaToGrid(const char *xFN, const char *yFN, std::function<float(float)> lambda) 
{ return ApplyLambdaToGrid(std::string(xFN), std::string(yFN), lambda); }
bool ApplyLambdaToGrid(const std::string & xFN, const std::string & yFN, std::function<float(float)> lambda)
{
	try {
		int nCols, nRows;
		float noData;
		std::ifstream xGridFS(GridFN(xFN), std::ios::binary);
		std::ofstream yGridFS(GridFN(yFN), std::ios::binary);
		if (!(xGridFS.is_open() && yGridFS.is_open())) return false;
		if (!GetHeader(HeaderFN(xFN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(xFN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(xFN), "NODATA_value", noData)) return false;
		if (!CopyFS(HeaderFN(xFN), HeaderFN(yFN), true)) return false;
		std::unique_ptr<float[]> dataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			xGridFS.read((char *)dataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				if (dataPtr[c] == noData) continue;
				dataPtr[c] = lambda(dataPtr[c]);
			}
			yGridFS.write((const char *)dataPtr.get(), nCols * sizeof(float));
		}
		return true;
	}
	catch (...) { return false; }
}

//Applies the provided lambda to each non-nodata value in ioFN then overwrite with the result 
bool ApplyLambdaToGrid(const char *ioFN, std::function<float(float)> lambda) 
{ return ApplyLambdaToGrid(std::string(ioFN), lambda); }
bool ApplyLambdaToGrid(const std::string & ioFN, std::function<float(float)> lambda) 
{
	try {
		std::string tmpFN = RemoveFileExt(ioFN) + "_LambdaTemp";
		if (!ApplyLambdaToGrid(ioFN, tmpFN, lambda)) return false;
		if (!CopyFS(GridFN(tmpFN), GridFN(ioFN), true)) return false;
		DeleteGrid(tmpFN);
		return true;
	}
	catch (...) { return false; }
}

//Applies the provided lambda to each non-nodata value in aFN and bFN and writes results to cFN 
bool ApplyLambdaToTwoGrids(const char *aFN, const char *bFN, const char *cFN, std::function<float(float, float)> lambda) 
{ 
	return ApplyLambdaToTwoGrids(std::string(aFN), std::string(bFN), std::string(cFN), lambda); 
}
bool ApplyLambdaToTwoGrids(const std::string & aFN, const std::string & bFN, const std::string & cFN, std::function<float(float, float)> lambda) 
{
	try {
		int nCols, nRows;
		float noData;
		std::ifstream aGridFS(GridFN(aFN), std::ios::binary);
		std::ifstream bGridFS(GridFN(bFN), std::ios::binary);
		std::ofstream cGridFS(GridFN(cFN), std::ios::binary);
		if (!(aGridFS.is_open() && bGridFS.is_open() && cGridFS.is_open())) return false;
		if (!CompareGridHeaders(HeaderFN(aFN), HeaderFN(bFN))) return false;
		if (!GetHeader(HeaderFN(aFN), "ncols", nCols)) return false;
		if (!GetHeader(HeaderFN(aFN), "nrows", nRows)) return false;
		if (!GetHeader(HeaderFN(aFN), "NODATA_value", noData)) return false;
		if (!CopyFS(HeaderFN(aFN), HeaderFN(cFN), true)) return false;
		std::unique_ptr<float[]> aDataPtr = std::make_unique<float[]>(nCols);
		std::unique_ptr<float[]> bDataPtr = std::make_unique<float[]>(nCols);
		std::unique_ptr<float[]> cDataPtr = std::make_unique<float[]>(nCols);
		for (int r = 0; r < nRows; r++) {
			aGridFS.read((char *)aDataPtr.get(), nCols * sizeof(float));
			bGridFS.read((char *)bDataPtr.get(), nCols * sizeof(float));
			for (int c = 0; c < nCols; c++) {
				cDataPtr[c] = (aDataPtr[c] == noData || bDataPtr[c] == noData) ? noData : lambda(aDataPtr[c], bDataPtr[c]);
			}
			cGridFS.write((const char *)cDataPtr.get(), nCols * sizeof(float));
		}
		return true;
	}
	catch (...) { return false; }
}

//Applies the provided lambda to each non-nodata value in aFN and bFN and writes results to cFN 
//Using gstream objects
//TODO Best to wrap all FltGridFile functionality to gstream class
bool ApplyLambdaToTwoGStreams(const char *aFN, const char *bFN, const char *cFN, std::function<float(float, float)> lambda)
{
	return ApplyLambdaToTwoGStreams(std::string(aFN), std::string(bFN), std::string(cFN), lambda);
}
bool ApplyLambdaToTwoGStreams(const std::string & aFN, const std::string & bFN, const std::string & cFN, std::function<float(float, float)> lambda)
{
	try {
		uint nCols, nRows;
		double cellSize;
		float noData;
		igstream aGridFS(aFN);
		igstream bGridFS(bFN);
		if (!aGridFS.compareHeader(bGridFS)) return false;
		ogstream cGridFS(cFN);
		cGridFS.copyHeader(aGridFS);
		aGridFS.getHeader(nCols, nRows, cellSize, noData);
		//if (!(aGridFS.is_open() && bGridFS.is_open() && cGridFS.is_open())) return false;
		std::unique_ptr<float[]> aDataPtr = std::make_unique<float[]>(nCols);
		std::unique_ptr<float[]> bDataPtr = std::make_unique<float[]>(nCols);
		std::unique_ptr<float[]> cDataPtr = std::make_unique<float[]>(nCols);
		for (uint r = 0; r < nRows; r++) {
			aGridFS.read((char *)aDataPtr.get(), nCols * sizeof(float));
			bGridFS.read((char *)bDataPtr.get(), nCols * sizeof(float));
			for (uint c = 0; c < nCols; c++) {
				cDataPtr[c] = (aDataPtr[c] == noData || bDataPtr[c] == noData) ? noData : lambda(aDataPtr[c], bDataPtr[c]);
			}
			cGridFS.write((const char *)cDataPtr.get(), nCols * sizeof(float));
		}
		return true;
	}
	catch (...) { return false; }
}




//Compare two grids storing the mean absolute cell value difference in meanDif if supplied
bool CompareGrids(const std::string & aFN, const std::string & bFN)
{
	float meanDif;
	return CompareGrids(aFN, bFN, meanDif);
}

bool CompareGrids(const std::string & aFN, const std::string & bFN, float &meanDif)
{
	int nCols, nRows;
	float noData;
	double sumDif = 0.0, cellDif;
	std::ifstream aGridFS(GridFN(aFN), std::ios::binary);
	std::ifstream bGridFS(GridFN(bFN), std::ios::binary);
	if (!aGridFS.is_open() || !bGridFS.is_open()) return false;
	if (!CompareGridHeaders(HeaderFN(aFN), HeaderFN(bFN))) return false;
	if (!GetHeader(HeaderFN(aFN), "ncols", nCols)) return false;
	if (!GetHeader(HeaderFN(aFN), "nrows", nRows)) return false;
	if (!GetHeader(HeaderFN(aFN), "NODATA_value", noData)) return false;
	std::unique_ptr<float[]> aDataPtr = std::make_unique<float[]>(nCols);
	std::unique_ptr<float[]> bDataPtr = std::make_unique<float[]>(nCols);
	for (int r = 0; r < nRows; r++) {
		aGridFS.read((char *)aDataPtr.get(), nCols * sizeof(float));
		bGridFS.read((char *)bDataPtr.get(), nCols * sizeof(float));
		for (int c = 0; c < nCols; c++) {
			if (aDataPtr[c] != bDataPtr[c]) {
				cellDif = abs(double(aDataPtr[c]) - double(bDataPtr[c]));
				if (cellDif > DBL_MAX - sumDif) {
					sumDif = DBL_MAX;
					break;
				}
				sumDif += cellDif;
			}
		}
	}
	if (sumDif != 0.0) {
		meanDif = float(sumDif / (ull(nCols) * ull(nRows)));
		return false;
	}
	meanDif = 0.0f;
	return true;
}
