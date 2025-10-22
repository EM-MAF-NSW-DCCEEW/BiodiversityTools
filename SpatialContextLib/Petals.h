/*
Petals.h - Updated functions for creating CBA petal files and data structures
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
#ifndef _PETALS_H_
#define _PETALS_H_
//#include <cuda_runtime.h>
//#include <device_launch_parameters.h>
#include "Parameters.h"

int CreatePetals(const std::string & outFN, int srcSize, int dstSize, float zoneRatio, bool createFiles, CBAParams &p);


//Old style functions
bool CreatePetalData(const char *translateFN, const char *lookupFN, const char *baseFN, CBAParams &p);
bool CreatePetalData(const std::string & translateFN, const std::string & lookupFN, const std::string  & baseFN, CBAParams &p);

int GetPetalFiles(const std::string & paramFN, const std::string & paramSection, const std::string workDir, std::string & lookupFN, std::string & translateFN, CBAParams &p);
int GetPetalFiles(const std::string & paramFN, const std::string & paramSection, std::string & lookupFN, std::string & translateFN, CBAParams &p);

int CreateSegments(int radCells, int nRings, int nSlices, double zoneRatio, const std::string & outFN, CBAParams &p);
int CreateSegments(int radCells, int nRings, int nSlices, std::vector<double> radList, const std::string & outFN, CBAParams &p);


//Search window functions
//int CreateSearchWindow(const std::string & paramFN);
//int CreateSearchWindow(const std::string & paramFN, const std::string & paramSection);
int CreateSearchWindow(const std::string & paramFN, const std::string & paramSection, CBAParams & p);


//Petal functions
int CreateLookupAndTranslateFromFiles(std::map<int, std::vector<int2>>& lookup, std::map<int, int2>& translate, uint & srcSize, uint & dstSize, const std::string & lookupFN, const std::string & translateFN);
int CreateLookupAndTranslateFromParameters(std::map<int, std::vector<int2>>& lookup, std::map<int, int2>& translate, const uint & srcSize, const uint & dstSize, const double & zoneRatio);
int CreateGPUDataFromLookupAndTranslate(CBAParams& p, const std::map<int, std::vector<int2>>& lookup, const std::map<int, int2>& translate, const uint& srcSize, const uint& dstSize);

//Segment functions
int CreateSegmentRadListFromParameters(std::vector<double>& radList, const int & radCells, const int & nBands, const double & zoneRatio);
int CreateLookupFromSegmentRadList(std::map<int, std::vector<int2>>& lookup, uint& srcSize, const std::vector<double>& radList);
bool GetSegmentNeighbours(const int segmentID, const int nBands, int * n);
int CreateGPUDataFromSegmentLookup(CBAParams & p, const std::map<int, std::vector<int2>>& lookup, const int & nBands, const uint & srcSize);

//Window grid functions
int CreateWindowGridFromLookup(std::unique_ptr<int[]>& windowGrid, const std::map<int, std::vector<int2>>& lookup, const uint& srcSize);
int WriteWindowGridToFile(const std::unique_ptr<int[]>& windowGrid, const std::string& windowGridFN, const uint& srcSize);

int WritePetalDataToTextFile(const CBAParams& p, const std::string& petalDataFN);


//Deleted functions
//int CreateWindowGridFromGridFile(std::unique_ptr<int[]>& windowGrid, uint& srcSize, uint& dstSize, const std::string& windowGridFN);
//int CreateWindowGridFromSegmentRadList(std::unique_ptr<int[]>& windowGrid, uint & srcSize, const std::vector<double>& radList);
//int CreateLookupFromWindowGrid(std::map<int, std::vector<int2>>& lookup, const std::unique_ptr<int[]>& windowGrid, const uint& srcSize);


#endif
