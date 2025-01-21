/*
FltGridFile.h - Various ESRI format floating point grid (.flt,.hdr) related functions
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
#ifndef FLT_GRID_FILE_H
#define FLT_GRID_FILE_H
#include <functional>

//MACRO for use with ApplyLambdaToGrid function
#define F2FLAMBDA(statement) [&](float x) -> float {return (float) statement }
#define FF2FLAMBDA(statement) [&](float a, float b) -> float {return (float) statement }

typedef unsigned long long ull;
typedef unsigned int uint;

//Grid File Name Functions

//Returns inFN with .flt file extension
std::string GridFN(const std::string & inFN);

//Returns inFN with .hdr file extension
std::string HeaderFN(const std::string & inFN);

//Returns inFN with .prj file extension
std::string PrjFN(const std::string & inFN);

//Returns inFN with .tif file extension
std::string TifFN(const std::string & inFN);

//Get flt and hdr filename 
bool GetFltHdrFileNames(const char*inFN, char *gridFN, char *hdrFN);
bool GetFltHdrFileNames(const std::string & inFN, std::string &gridFN, std::string &hdrFN);

//Delete all grid files associated with inFN
bool DeleteGrid(const std::string &inFN);

//Grid Header Functions
//Get a string from a flt grid header file
bool GetHeader(const char *hdrFN, const char *keyword, char *value);
bool GetHeader(const std::string & hdrFN, const std::string & keyword, std::string & value);

//Get a integer from a flt grid header file
bool GetHeader(const char *hdrFN, const char *keyword, int &value);
bool GetHeader(const std::string & hdrFN, const std::string & keyword, int &value);

//Get an unsigned integer from a flt grid header file
bool GetHeader(const char *hdrFN, const char *keyword, unsigned int &value);
bool GetHeader(const std::string & hdrFN, const std::string & keyword, unsigned int &value);

//Get a float from a flt grid header file
bool GetHeader(const char *hdrFN, const char *keyword, float &value);
bool GetHeader(const std::string & hdrFN, const std::string & keyword, float &value);

//Get a float from a flt grid header file
bool GetHeader(const char *hdrFN, const char *keyword, double &value);
bool GetHeader(const std::string & hdrFN, const std::string & keyword, double &value);

//Read a flt grid header file
bool ReadGridHeader(const char *hdrFN, unsigned int &nCols, unsigned int &nRows, double &xll, double &yll, double &cellSize, float &noData);
bool ReadGridHeader(const std::string & hdrFN, unsigned int &nCols, unsigned int &nRows, double &xll, double &yll, double &cellSize, float &noData);

//Read a flt grid header file ignoring xll and yll
bool ReadGridHeader(const char *hdrFN, unsigned int &nCols, unsigned int &nRows, double &cellSize, float &noData);
bool ReadGridHeader(const std::string & hdrFN, unsigned int &nCols, unsigned int &nRows, double &cellSize, float &noData);

//Write a flt grid header file
bool WriteGridHeader(const char *hdrFN, unsigned int nCols, unsigned int nRows, double xll, double yll, double cellSize, float noData);
bool WriteGridHeader(const std::string & hdrFN, unsigned int nCols, unsigned int nRows, double xll, double yll, double cellSize, float noData);

//Write ascii header to file
bool WriteAsciiHeader(const char *hdrFN, unsigned int nCols, unsigned int nRows, double xll, double yll, double cellSize, float noData);
bool WriteAsciiHeader(const std::string & hdrFN, unsigned int nCols, unsigned int nRows, double xll, double yll, double cellSize, float noData);

//Compare two flt grid header files
bool CompareGridHeaders(const char *hdrFNA, const char *hdrFNB);
bool CompareGridHeaders(const std::string & hdrFNA, const std::string & hdrFNB);
bool CompareGridHeadersNoExt(const std::string & hdrFNA, const std::string & hdrFNB);

//Grid Data Functions

//Get the min value from a flt grid
bool GetGridMin(const char *gridFN, float &min);
bool GetGridMin(const std::string & FN, float & min);

//Get the max value from a flt grid
bool GetGridMax(const char *gridFN, float &max);
bool GetGridMax(const std::string & FN, float & max);

//Get the smallest positive power of 10 greater or equal to the highest grid value
bool GetGridMaxPow10(const std::string gridFN, float &maxPow10);

//Get the min and max value from a flt grid
bool GetGridMinMax(const char *gridFN, float &min, float &max);
bool GetGridMinMax(const std::string & FN, float & min, float & max);

//Get the sum of values from a flt grid
bool GetGridSum(const char *FN, double &sum);
bool GetGridSum(const std::string & FN, double &sum);

//Get the average of values from a flt grid
bool GetGridAvg(const char *FN, double &avg);
bool GetGridAvg(const std::string & FN, double &avg);

//Get the number of data cells in a flt grid
bool CountDataCells(const std::string & FN, ull & nDataCells);

//Rescale flt grid values to between outMin and outMax
bool RescaleGrid(const char *inFN, const char *outFN, float outMin, float outMax);
bool RescaleGrid(const std::string & inFN, const std::string & outFN, float outMin, float outMax);
bool RescaleGrid(const std::string & inFN, const std::string & outFN, float inMin, float inMax, float outMin, float outMax);

//Create a permeability flt grid from a habitat grid
bool CreatePrmGrid(const char *habFN, const char *prmFN, float distMin, float distMax);
bool CreatePrmGrid(const std::string & habFN, const std::string & prmFN, float distMin, float distMax);

//Create a circle flt grid of area in hectares, optional grid dimension dim
bool CreateCircleGrid(const char *outFN, double cellSize, double areaHA, float inVal, float outVal);
bool CreateCircleGrid(const std::string & outFN, double cellSize, double areaHA, float inVal, float outVal);
bool CreateCircleGrid(const char *outFN, double cellSize, double areaHA, int dim, float inVal, float outVal);
bool CreateCircleGrid(const std::string & outFN, double cellSize, double areaHA, int dim, float inVal, float outVal);

//Create a random point grid using a habitat weight grid
bool CreateRandomPointGrid(const char *inFN, const char *outFN, float p, float minVal);
bool CreateRandomPointGrid(const std::string & inFN, const std::string & outFN, float p, float minVal);

//Resample a grid using GDAL (method:) and output cellsize
//Nearest neighbour=0, Bilinear=1, Cubic Convolution=2, Cubic B-Spline=3, Lanczos=4,
//Average=5, Mode=6, GRA_Gauss=7. Max=8, Min=9, Med=10, Q1=11, Q3=12
bool ResampleGrid(const char *inFN, const char *outFN, double outCellSize, int method = 0);

//Nearest neighbour=0, Bilinear=1, Cubic Convolution=2, Cubic B-Spline=3, Lanczos=4,
//Average=5, Mode=6, GRA_Gauss=7. Max=8, Min=9, Med=10, Q1=11, Q3=12
bool ResampleGrid(const std::string & inFN, const std::string & outFN, double outCellSize, int method);

//Resample a grid using GDAL (method:) and template grid dimensions
//Nearest neighbour=0, Bilinear=1, Cubic Convolution=2, Cubic B-Spline=3, Lanczos=4,
//Average=5, Mode=6, GRA_Gauss=7. Max=8, Min=9, Med=10, Q1=11, Q3=12
bool ResampleGrid(const char *inFN, const char *outFN, const char *tmpFN, int method = 0);

//Nearest neighbour=0, Bilinear=1, Cubic Convolution=2, Cubic B-Spline=3, Lanczos=4,
//Average=5, Mode=6, GRA_Gauss=7. Max=8, Min=9, Med=10, Q1=11, Q3=12
bool ResampleGrid(const std::string & inFN, const std::string & outFN, const std::string & tmpFN, int method);

//Applies the provided function to each non-nodata value in xFN and writes results to yFN 
bool ApplyFuncToGrid(const char *xFN, const char *yFN, float(*func)(float));
bool ApplyFuncToGrid(const std::string & xFN, const std::string & yFN, float(*func)(float));

//Applies the provided lambda to each non-nodata value in xFN and writes results to yFN 
//Usage example: ApplyLambdaToGrid(inFN, outFN, F2FLAMBDA(x * 2;))
bool ApplyLambdaToGrid(const char *xFN, const char *yFN, std::function<float(float)> lambda);
bool ApplyLambdaToGrid(const std::string & xFN, const std::string & yFN, std::function<float(float)> lambda);

//Applies the provided lambda to each non-nodata value in xFN and overwrites with results 
//Usage example: ApplyLambdaToGrid(ioFN, F2FLAMBDA(x * 2;))
bool ApplyLambdaToGrid(const char *ioFN, std::function<float(float)> lambda);
bool ApplyLambdaToGrid(const std::string & ioFN, std::function<float(float)> lambda);

//Applies the provided lambda to each non-nodata value in aFN and bFN and writes results to cFN 
//Usage example: ApplyLambdaToGrid(ioFN, FF2FLAMBDA(a * b * 2;))
bool ApplyLambdaToTwoGrids(const char *aFN, const char *bFN, const char *cFN, std::function<float(float, float)> lambda);
bool ApplyLambdaToTwoGrids(const std::string & aFN, const std::string & bFN, const std::string & cFN, std::function<float(float, float)> lambda);
//Using gstream objects
//TODO Best to wrap all FltGridFile functionality to gstream class
bool ApplyLambdaToTwoGStreams(const char *aFN, const char *bFN, const char *cFN, std::function<float(float, float)> lambda);
bool ApplyLambdaToTwoGStreams(const std::string & aFN, const std::string & bFN, const std::string & cFN, std::function<float(float, float)> lambda);

//Compare two grids storing the mean absolute cell value difference in meanDif if supplied
bool CompareGrids(const std::string &aFN, const std::string &bFN);
bool CompareGrids(const std::string &aFN, const std::string &bFN, float &meanDif);
#endif
