/*
Petals.cpp - Updated functions for creating CBA petal files and data structures
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

#define _USE_MATH_DEFINES

#include <cstdio>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <ios>
#include <map>
#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <climits>
#include <cfloat>
#include <vector_types.h>

#include "FileUtilsLib.h"
#include "Petals.h"


#ifndef _WIN32
#	define min std::min
#	define max std::max
#endif
//TODO no longer need to create or use petal lookup and translate files.
//TODO replace ints with uints where appropriate

//Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int CreateSearchWindow(const std::string & paramFN, const std::string & paramSection, CBAParams &p)
{
	std::string lookupFN, translateFN, windowGridFN, petalDataFN;
	uint srcSize, dstSize;
	double zoneRatio;
	int radCells, nRings;
	int returnValue = 0;

	//Data structures
	std::map<int, std::vector<int2>> lookup;	//Petal ID, vector(pair(cell.x, cell.y))
	std::map<int, int2> translate; //Petal ID, pair(dst.x, dst.y)
	std::vector<double> radList;
	std::unique_ptr<int[]> windowGrid;


	//Read params to work out what we've got
	//Create petals from lookup and translate files
	if (GetParam(paramFN, paramSection, "LOOKUP_FN", lookupFN) && FileExists(lookupFN) &&
		GetParam(paramFN, paramSection, "TRANSLATE_FN", translateFN) && FileExists(translateFN)) {
	
		CreateLookupAndTranslateFromFiles(lookup, translate, srcSize, dstSize, lookupFN, translateFN);
		CreateGPUDataFromLookupAndTranslate(p, lookup, translate, srcSize, dstSize);
	
	}
	//Create petals from parameters
	else 
	if (GetParam(paramFN, paramSection, "SRCSIZE", srcSize) &&
		GetParam(paramFN, paramSection, "DSTSIZE", dstSize) &&
		GetParam(paramFN, paramSection, "ZONERATIO", zoneRatio)) {

		CreateLookupAndTranslateFromParameters(lookup, translate, srcSize, dstSize, zoneRatio);
		CreateGPUDataFromLookupAndTranslate(p, lookup, translate, srcSize, dstSize);
	}
	//Create segments from parameters
	else if (GetParam(paramFN, paramSection, "RADIUS", radCells) &&
		GetParam(paramFN, paramSection, "N_RINGS", nRings) &&
		GetParam(paramFN, paramSection, "ZONERATIO", zoneRatio)) {

		CreateSegmentRadListFromParameters(radList, radCells, nRings, zoneRatio);
		CreateLookupFromSegmentRadList(lookup, srcSize, radList);
		CreateGPUDataFromSegmentLookup(p, lookup, nRings, srcSize);
	}
	//Create segments from radius list
	else if (GetParam(paramFN, paramSection, "RADIUSLIST", radList)) {
		nRings = (int)radList.size();
		if (nRings < 1) return msgErrorStr(-1, "RADIUSLIST has no values", paramFN);
		CreateLookupFromSegmentRadList(lookup, srcSize, radList);
		CreateGPUDataFromSegmentLookup(p, lookup, nRings, srcSize);
	}
	else return msgErrorStr(-1, "for search window " + paramSection, paramFN);

	if (GetParam(paramFN, paramSection, "WINDOWGRID_FN", windowGridFN)) {
		CreateWindowGridFromLookup(windowGrid, lookup, srcSize);
		WriteWindowGridToFile(windowGrid, windowGridFN, srcSize);
	}

	if (GetParam(paramFN, paramSection, "PETALDATA_FN", petalDataFN)) {
		WritePetalDataToTextFile(p, petalDataFN);
	}

	return returnValue;
}

//Create lookup & translate from files
int CreateLookupAndTranslateFromFiles(std::map<int, std::vector<int2>> & lookup, std::map<int, int2> & translate, uint & srcSize, uint & dstSize, const std::string & lookupFN, const std::string & translateFN)
{
	std::string line, idStr, xStr, yStr;
	int nPetals = 0;
	int xMax = 0, yMax = 0;

	//Read lookup file
	std::ifstream lookupFS(lookupFN, std::ios::in);
	if (!lookupFS.is_open()) return msgErrorStr(-3, lookupFN);
	while (std::getline(lookupFS, line)) {
		std::stringstream lss(line);
		std::getline(lss, idStr, ',');
		std::getline(lss, xStr, ',');
		std::getline(lss, yStr);
		lookup[std::stoi(idStr)].push_back({ std::stoi(xStr), std::stoi(yStr) });
		nPetals = max(nPetals, std::stoi(idStr));
		xMax = max(xMax, (uint) std::stoi(xStr));
		yMax = max(yMax, (uint) std::stoi(yStr));
	}

	//Read translate file
	std::ifstream translateFS(translateFN, std::ios::in);
	if (!translateFS.is_open()) return msgErrorStr(-3, translateFN);
	while (std::getline(translateFS, line)) {
		std::stringstream lss(line);
		std::getline(lss, idStr, ',');
		std::getline(lss, xStr, ',');
		std::getline(lss, yStr);
		translate[std::stoi(idStr)] = { std::stoi(xStr), std::stoi(yStr) };
	}

	srcSize = uint(xMax) + 1;
	double dstSizeDbl = std::sqrt(double(nPetals + 1));
	dstSize = uint(dstSizeDbl);

	//Check petal data
	if (nPetals < 1) return msgErrorStr(-99, "No petals in", lookupFN);
	if (lookup.size() != nPetals) return msgErrorStr(-99, "lookup.size " + toStr(lookup.size()) + " != nPetals " + toStr(nPetals));
	if (xMax < 2 || yMax < 2) return msgErrorStr(-99, "Source window too small " + toStr(xMax) +"x" + toStr(yMax));
	if (xMax != yMax) return msgErrorStr(-99, "Source window not square " + toStr(xMax) + "x" + toStr(yMax));
	if (srcSize % 2 != 1) return msgErrorStr(-99, "Source window size " + toStr(srcSize) + " must be odd");
	if (std::floor(dstSizeDbl) != dstSizeDbl) return msgErrorStr(-99, "Destination window not square " + toStr(nPetals + 1));
	if (dstSize % 2 != 1) return msgErrorStr(-99, "Destination window size " + toStr(dstSize) + " must be odd");
	if (translate.size() != lookup.size() + 1) return msgErrorStr(-99, "Translate and lookup files don't match", toStr(translate.size()), toStr(lookup.size()));

	return 0;
}

//Create lookup & translate from petal parameters
int CreateLookupAndTranslateFromParameters(std::map<int, std::vector<int2>> & lookup, std::map<int, int2> & translate, const uint & srcSize, const uint & dstSize, const double & zoneRatio)
{
	//Check numeric parameters
	if (zoneRatio < 1.0) return msgErrorStr(-99, "Zone ratio must be greater or equal to 1.0");
	if (srcSize % 2 != 1) return msgErrorStr(-99, "Source window size must be odd");
	if (dstSize % 2 != 1) return msgErrorStr(-99, "Destination window size must be odd");
	if (srcSize < dstSize) return msgErrorStr(-99, "Destination window size must be less than or equal to source window size");

	int x, y;
	int srcCentre = int(srcSize / 2);
	int dstCentre = int(dstSize / 2);
	int nPetals = dstSize * dstSize - 1;

	//Create zonefactors
	float rawSum = 0.0f, zoneSize = 0.0f, scaledSum = 0.0f;
	std::map<int, float>zoneFactors;//Zone, factor
	for (int zoneNo = 0; zoneNo < dstCentre; zoneNo++)	rawSum += (float) pow(zoneRatio, zoneNo);
	for (int zoneNo = 1; zoneNo < dstCentre + 1; zoneNo++) {
		zoneSize = std::fmaxf(srcCentre / rawSum * (float) pow(zoneRatio, zoneNo - 1), 1.0f);
		scaledSum += (zoneNo > 1 && zoneSize + scaledSum > srcCentre) ? srcCentre - scaledSum : zoneSize;
		zoneFactors[zoneNo] = std::fmax(scaledSum / zoneNo, 1.0f);
		if (zoneNo > 1 && zoneFactors[zoneNo] * zoneNo <= zoneFactors[zoneNo - 1] * (zoneNo - 1)) return msgErrorStr(-99, "Unable to calculate zone factors");
	}

	//Set a unique petalID for each petal starting from 1
	int petalID = 1;
	int zone, zoneOffset, dstY, dstX;
	float factor, prevOffset, splitOffset, inc, step;
	std::map<int, int2> petalPoints; 
	
	translate[0] = { dstCentre, dstCentre };

	for (auto zoneFactor : zoneFactors) {
		zone = zoneFactor.first;
		factor = zoneFactor.second;
		dstY = dstCentre - zone;
		dstX = dstCentre - zone;

		//Calculate zone offset and increment
		if (zone > 1) {
			prevOffset = ((zone - 1) * zoneFactors[zone - 1]);
			splitOffset = (zone * factor - prevOffset) / 2;
			zoneOffset = int(ceil(prevOffset + splitOffset));
		}
		else zoneOffset = int(ceil(zone * factor));
		inc = float(zoneOffset) / float(zone);

		//Top horizontal x increasing
		y = srcCentre - zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			x = srcCentre - zoneOffset + int(round(step));
			translate[petalID] = { dstX++, dstY };
			petalPoints[petalID] = { x, y };
			petalID++;
		}

		//Right vertical y increasing
		x = srcCentre + zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			y = srcCentre - zoneOffset + int(round(step));
			translate[petalID] = { dstX, dstY++ };
			petalPoints[petalID] = { x, y };
			petalID++;
		}

		//Bottom horizontal x decreasing
		y = srcCentre + zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			x = srcCentre + zoneOffset - int(round(step));
			translate[petalID] = { dstX--, dstY };
			petalPoints[petalID] = { x, y };
			petalID++;
		}

		//Left vertical y decreasing
		x = srcCentre - zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			y = srcCentre + zoneOffset - int(round(step));
			translate[petalID] = { dstX, dstY-- };
			petalPoints[petalID] = { x, y };
			petalID++;
		}
	}

	////Create lookup from petal points
	//double distSqrd = 0.0f, minDistSqrd = 0.0f;
	//int petalID = -1;
	//for (y = 0; y < srcSize; y++) {
	//	for (x = 0; x < srcSize; x++) {
	//		if (x != srcCentre || y != srcCentre) {
	//			minDistSqrd = FLT_MAX;
	//			for (const auto petalPoint : petalPoints) {
	//				if (petalPoint.second.x == x && petalPoint.second.y == y) {
	//					petalID = petalPoint.first;
	//				}
	//				else {
	//					distSqrd = pow(abs(petalPoint.second.x - x), 2) + pow(abs(petalPoint.second.y - y), 2);
	//					if (distSqrd < minDistSqrd) {
	//						minDistSqrd = distSqrd;
	//						petalID = petalPoint.first;
	//					}
	//				}
	//			}
	//			lookup[petalID].push_back({ x, y });
	//		}
	//	}
	//}

	//Create lookup from petal points
	double distSqrd = 0.0f, minDistSqrd = 0.0f;
	petalID = -1;
	for (y = 0; y < srcSize; y++) {
		for (x = 0; x < srcSize; x++) {
			if (x != srcCentre || y != srcCentre) {
				petalID = -1;
				minDistSqrd = FLT_MAX;
				for (const auto petalPoint : petalPoints) {
					distSqrd = pow(abs(petalPoint.second.x - x), 2) + pow(abs(petalPoint.second.y - y), 2);
					if (distSqrd < minDistSqrd) {
						minDistSqrd = distSqrd;
						petalID = petalPoint.first;
					}
				}
				if (petalID < 0) return msgErrorStr(-99, "Petal ID < 1. Error in petal parameters.");
				if (petalID > nPetals) return msgErrorStr(-99, "Petal ID >= nPetals. Error in petal parameters.");
				lookup[petalID].push_back({ x, y });
			}
		}
	}

	return 0;
}


////Create window grid from grid file
//int CreateWindowGridFromGridFile(std::unique_ptr<int[]> & windowGrid, uint & srcSize, uint & dstSize, const std::string & windowGridFN)
//{
//	igstream windowGridGS(windowGridFN);
//	if (!windowGridGS.is_open()) return msgErrorStr(-3, windowGridFN);
//	if (windowGridGS.dataType() != 5) return msgErrorStr(-99, "Data type is not int32", windowGridFN);
//	
//	uint nCols = windowGridGS.nCols();
//	uint nRows = windowGridGS.nRows();
//	if (nCols < 1 || nRows < 1) return msgErrorStr(-99, "0 columns or rows in", windowGridFN);
//	if (nCols != nRows) return msgErrorStr(-99, "Source window not square " + windowGridFN, "\nCols: " + toStr(nCols), " Rows: " + toStr(nRows));
//	if (nCols % 2 != 1) return msgErrorStr(-99, "Source window size is even", toStr(srcSize));
//	uint gridSize = nCols * nRows;
//	srcSize = uint(nCols);
//
//	windowGrid = std::make_unique<int[]>(gridSize);
//	windowGridGS.read((char *)(windowGrid.get()), gridSize * sizeof(int));
//	int nPetals = 0;
//	for (int i = 0; i < gridSize; i++) nPetals = max(nPetals, windowGrid[i]);
//
//	//long double sr = sqrt(long double(nPetals + 1));	 error: expected primary-expression before 'long'
//	double sr = sqrt(double(nPetals + 1));
//	if ((sr - std::floor(sr)) != 0.0L) return msgErrorStr(-99, "Destination window not square", toStr(nPetals + 1));
//	dstSize = uint(sr);
//	if (dstSize % 2 != 1) return msgErrorStr(-99, "Destination window size is even", toStr(dstSize));
//	return 0;
//}


//Create segment radius list from parameters
int CreateSegmentRadListFromParameters(std::vector<double> & radList, const int & radCells, const int & nBands, const double & zoneRatio)
{
	if (zoneRatio < 1.0f) return msgErrorStr(-99, "Zone ratio must be greater or equal to 1.0");
	if (radCells < 1) return msgErrorStr(-99, "Radius in cells must be greater than 1");
	if (nBands < 1) return msgErrorStr(-99, "Number of rings is less than 1");
	if (radCells < nBands) return msgErrorStr(-99, "Radius in cells is less than number of rings");

	//Create ring area list with area within each ring
	double areaFactor = (pow(radCells, 2.0) * M_PI) / ((zoneRatio * (pow(zoneRatio, nBands) - 1.0)) / (zoneRatio - 1.0));
	std::vector<double> areaList;
	areaList.push_back(0.0);
	for (int r = 1; r < nBands + 1; r++) areaList.push_back(pow(zoneRatio, r) * areaFactor);

	//Create ring radius list with outer radius of each ring
	double areaSum = 0.0;
	for (int r = 1; r < nBands + 1; r++) {
		radList.push_back(sqrt((areaList[r] + areaSum) / M_PI));
		areaSum += areaList[r];
	}
	return 0;
}

////Create window grid from segment radius list
//int CreateWindowGridFromSegmentRadList(std::unique_ptr<int[]> & windowGrid, uint & srcSize, const std::vector<double> & radList)
//{
//	//Check if first segments will be too small
//	if (radList[0] < sqrt(2)) return msgErrorStr(-99, "First radius is less than SQRT(2) and would result in missing segments.");
//	if (radList.size() < 1) return msgErrorStr(-99, "No rings in radius list.");
//
//	for (int r = 1; r < radList.size(); r++) {
//		if (radList[r] <= radList[r - 1]) return msgErrorStr(-99, "Radius list not strictly increasing.");
//	}
//
//
//	//TODO check nBands and radCells from radList
//	int nBands = (int) radList.size();
//	int radCells = (int) std::ceil(radList[nBands - 1]);
//	int nSectors = 8;
//	int nSegments = nBands * nSectors;
//	
//	srcSize = radCells * 2 + 1;
//	int halfSrcSize = srcSize / 2;
//
//	//Create window grid
//	windowGrid = std::make_unique<int[]>(srcSize * srcSize);
//	for (int i = 0; i < srcSize * srcSize; i++) windowGrid[i] = -1;
//
//	//Calc cell distance, degrees and segment ID to fill the lookup table, segment cell count and window grid 
//	int segmentID = -1;
//	double cellDist, cellDeg;
//	for (int y = -halfSrcSize; y < 1 + halfSrcSize; y++) {
//		for (int x = -halfSrcSize; x < 1 + halfSrcSize; x++) {
//			segmentID = -1;
//			cellDist = sqrt(pow(double(y), 2) + pow(double(x), 2));
//			cellDeg = 180.0 - atan2(double(y), double(x)) * 180.0 / M_PI;
//
//			if (cellDist < 0.5) segmentID = 0;
//			else {
//				for (int r = 0; r < radList.size(); r++) {
//					if (cellDist <= radList[r]) {
//						segmentID = 1 + (r * nSectors) + int(floor(nSectors * (359 - (int(cellDeg + (180.0 / nSectors)) % 360)) / 360.0));
//						break;
//					}
//				}
//			}
//			if (segmentID > nSegments) return msgErrorStr(-99, "Segment ID > nSegments. Error in segment parameters.");
//			windowGrid[(y + halfSrcSize) * srcSize + (x + halfSrcSize)] = segmentID;
//		}
//	}
//	return 0;
//}

//Create Lookup from segment radius list
int CreateLookupFromSegmentRadList(std::map<int, std::vector<int2>>& lookup, uint& srcSize, const std::vector<double>& radList)
{
	//Check if first segments will be too small
	if (radList.size() < 1) return msgErrorStr(-99, "No bands in radius list.");
	if (radList[0] < sqrt(2)) return msgErrorStr(-99, "First band radius is less than SQRT(2) and would result in missing segments.");
	for (int r = 1; r < radList.size(); r++) {
		if (radList[r] <= radList[r - 1]) return msgErrorStr(-99, "Radius list not strictly increasing.");
	}

	//TODO check nBands and radCells from radList
	int nBands = (int)radList.size();
	int radCells = (int)std::ceil(radList[nBands - 1]);
	int nSectors = 8;
	//int nSegments = nBands * nSectors;

	srcSize = radCells * 2 + 1;
	int halfSrcSize = srcSize / 2;

	lookup.clear();

	//Calc cell distance, degrees and segment ID to fill the lookup table, segment cell count and window grid 
	int segmentID = -1;
	double cellDist, cellDeg;
	for (int y = -halfSrcSize; y < 1 + halfSrcSize; y++) {
		for (int x = -halfSrcSize; x < 1 + halfSrcSize; x++) {
			segmentID = -1;
			cellDist = sqrt(pow(double(y), 2) + pow(double(x), 2));
			cellDeg = 180.0 - atan2(double(y), double(x)) * 180.0 / M_PI;

			if (cellDist < 0.5) 
				segmentID = 0;
			else {
				for (int r = 0; r < radList.size(); r++) {
					if (cellDist <= radList[r]) {
						segmentID = 1 + (r * nSectors) + int(floor(nSectors * (359 - (int(cellDeg + (180.0 / nSectors)) % 360)) / 360.0));
						break;
					}
				}

				if (segmentID != -1) {
					lookup[segmentID].push_back({ x + halfSrcSize, y + halfSrcSize });
				}
			}
			//if (segmentID >= nSegments) return msgErrorStr(-99, "Segment ID >= nSegments. Error in segment parameters.");
		}
	}
	return 0;
}


////Create lookup from window grid
//int CreateLookupFromWindowGrid(std::map<int, std::vector<int2>> & lookup, const std::unique_ptr<int[]> & windowGrid, const uint & srcSize) {
//	int x, y;
//	lookup.clear();
//	for (y = 0; y < srcSize; y++) {
//		for (x = 0; x < srcSize; x++) {
//			if (windowGrid[(y*srcSize) + x] > 0)
//				lookup[windowGrid[(y*srcSize) + x]].push_back({ x, y });
//		}
//	}
//	return 0;
//}

//Create window grid from lookup
int CreateWindowGridFromLookup(std::unique_ptr<int[]> & windowGrid, const std::map<int, std::vector<int2>> & lookup, const uint & srcSize)
{
	//Make window grid from lookup
	uint srcCentre = uint(srcSize / 2);
	windowGrid.reset();
	windowGrid = std::make_unique<int[]>(srcSize * srcSize);
	for (int c = 0; c < srcSize * srcSize; c++) windowGrid[c] = -1;
	windowGrid[srcCentre * srcSize + srcCentre] = 0;
	for (const auto p : lookup) {
		for (const auto c : p.second) {
			windowGrid[c.y * srcSize + c.x] = p.first;
		}
	}
	return 0;
}

//Create grid file from window grid
int WriteWindowGridToFile(const std::unique_ptr<int[]> & windowGrid, const std::string & windowGridFN, const uint & srcSize) {
	
	ogstream windowGridGS(TifFN(windowGridFN));
	windowGridGS.setTifOptions("LZW", 5, 1);
	windowGridGS.setHeader(srcSize, srcSize, 0.0, 0.0, 1.0, -1.0);
	if (!windowGridGS.is_open()) return msgErrorStr(-4, windowGridFN);
	windowGridGS.write((char *)windowGrid.get(), srcSize * srcSize * (sizeof(int)));
	return 0;
}
//Write petal data to text file
int WritePetalDataToTextFile(const CBAParams& p, const std::string& petalDataFN) {

	std::ofstream petalDataFS(petalDataFN, std::ios::out);
	if (!petalDataFS.is_open()) return msgErrorStr(-4, petalDataFN);
	for (uint j = 0; j < p.petalRows; j++)
		for (uint k = 0; k < p.petalCols; k++)
			petalDataFS << "" << p.petalPtr[j * p.petalCols + k].x << "," << p.petalPtr[j * p.petalCols + k].y << (k + 1 < p.petalCols ? "\t" : "\n");
	return 0;
}

////Write petal data to binary file
//bool WritePetalDataToBinFile(const CBAParams& p, const std::string& petalDataFN) {
//	std::ofstream petalDataFS(petalDataFN, std::ios::out | std::ios::binary);
//	if (!petalDataFS.is_open()) return false;
//	petalDataFS.write((char*)&p.fOffset, sizeof(int));
//	petalDataFS.write((char*)&p.nPetals, sizeof(int));
//	petalDataFS.write((char*)&p.petalCols, sizeof(int));
//	petalDataFS.write((char*)&p.petalRows, sizeof(int));
//	petalDataFS.write((char*)&p.firstReduction, sizeof(int));
//	petalDataFS.write((char*)p.petalPtr.get(), p.petalCols * p.petalRows * sizeof(int2));
//	petalDataFS.close();
//	return true;
//}
//
//bool ReadPetalDataFromBinFile(CBAParams& p, const std::string& petalDataFN) {
//	std::ifstream petalDataFS(petalDataFN, std::ios::in | std::ios::binary);
//	if (!petalDataFS.is_open()) return false;
//	petalDataFS.read((char*)&p.fOffset, sizeof(int));
//	petalDataFS.read((char*)&p.nPetals, sizeof(int));
//	petalDataFS.read((char*)&p.petalCols, sizeof(int));
//	petalDataFS.read((char*)&p.petalRows, sizeof(int));
//	petalDataFS.read((char*)&p.firstReduction, sizeof(int));
//	p.petalPtr = std::make_unique<int2[]>(p.petalCols * p.petalRows);
//	petalDataFS.read((char*)p.petalPtr.get(), p.petalCols * p.petalRows * sizeof(int2));
//	petalDataFS.close();
//	return true;
//}

//Create GPU data from lookup and translate
int CreateGPUDataFromLookupAndTranslate(CBAParams &p, const std::map<int, std::vector<int2>> & lookup, const std::map<int, int2> & translate, const uint & srcSize, const uint & dstSize) {
	
	//Set search window parameters and initialise petal data storage
	int maxCells = 1, cellsStart = 5, petalsRow;
	for (auto petal : lookup) maxCells = max(maxCells, int(petal.second.size()));
	p.fOffset = srcSize / 2;
	p.nPetals = dstSize * dstSize - 1;
	p.petalCols = int(ceil(float(p.nPetals) / 32.0f) * 32);
	p.petalRows = maxCells + cellsStart;
	p.firstReduction = int(powf(2, ceil(log2(float(p.petalCols))))) / 2;
	p.petalPtr = std::make_unique<int2[]>(p.petalCols * p.petalRows);
	for (uint i = 0; i < p.petalCols * p.petalRows; i++)	p.petalPtr[i] = { -1, -1 };

	//Fill petals data table
	for (auto pA : translate) {
		if (pA.first == 0) continue;
		for (auto pB : translate) {
			if (pB.first == 0) continue;
			//First four rows contain neighbouring petal indicies
			if (pB.second.x == pA.second.x && pB.second.y == pA.second.y - 1) p.petalPtr[pA.first - 1].x = pB.first - 1;
			else if (pB.second.x == pA.second.x && pB.second.y == pA.second.y + 1) p.petalPtr[pA.first - 1].y = pB.first - 1;
			else if (pB.second.x == pA.second.x - 1 && pB.second.y == pA.second.y) p.petalPtr[p.petalCols + pA.first - 1].x = pB.first - 1;
			else if (pB.second.x == pA.second.x + 1 && pB.second.y == pA.second.y) p.petalPtr[p.petalCols + pA.first - 1].y = pB.first - 1;
			else if (pB.second.x == pA.second.x - 1 && pB.second.y == pA.second.y - 1) p.petalPtr[2 * p.petalCols + pA.first - 1].x = pB.first - 1;
			else if (pB.second.x == pA.second.x + 1 && pB.second.y == pA.second.y - 1) p.petalPtr[2 * p.petalCols + pA.first - 1].y = pB.first - 1;
			else if (pB.second.x == pA.second.x - 1 && pB.second.y == pA.second.y + 1) p.petalPtr[3 * p.petalCols + pA.first - 1].x = pB.first - 1;
			else if (pB.second.x == pA.second.x + 1 && pB.second.y == pA.second.y + 1) p.petalPtr[3 * p.petalCols + pA.first - 1].y = pB.first - 1;
		}

		//From cellStart contains coords of cells in each petal
		petalsRow = cellsStart;
		for (auto lu : lookup.at(pA.first)) {
			p.petalPtr[petalsRow * p.petalCols + pA.first - 1].x = lu.x;
			p.petalPtr[petalsRow * p.petalCols + pA.first - 1].y = lu.y;
			petalsRow++;
		}

		//Fifth row contains number of cells in each petal and petal ID
		p.petalPtr[4 * p.petalCols + pA.first - 1].x = int(lookup.at(pA.first).size());
		p.petalPtr[4 * p.petalCols + pA.first - 1].y = pA.first - 1;
	}
	return 0;
}

bool GetSegmentNeighbours(const int segmentID, const int nBands, int *n) {
	try {
		n[0] = segmentID + 8;//n
		n[1] = segmentID - 8;//s
		n[2] = segmentID % 8 == 0 ? segmentID - 7 : segmentID + 1;//e
		n[3] = (segmentID - 1) % 8 == 0 ? segmentID + 7 : segmentID - 1;//w
		n[4] = n[2] + 8;//ne
		n[5] = n[3] + 8;//nw
		n[6] = n[2] - 8;//se
		n[7] = n[3] - 8;//sw
		if (n[0] > 8 * nBands) {
			n[0] = -1;//n
			n[4] = -1;//ne
			n[5] = -1;//nw
		}
		if (n[1] <= 0) {
			n[1] = -1;//s
			n[6] = -1;//se
			n[7] = -1;//sw
		}
		return true;
	}
	catch (std::exception e) { return false; }
}

//Create GPU data from lookup and neighbours
int CreateGPUDataFromSegmentLookup(CBAParams & p, const std::map<int, std::vector<int2>> & lookup, const int & nBands, const uint & srcSize) {
	
	uint nSegments = (uint) lookup.size(); // nBands * nSectors;
	uint nCols = uint(ceil(float(nSegments) / 32.0f) * 32);

	//Allocate CUDA CBA segment data with maxCells rows and init to (-1, -1)
	uint maxCells = 0;
	for (uint s = 0; s < nSegments; s++) maxCells = max(maxCells, (uint) lookup.at(s + 1).size());//(cellCounts[i]);
	uint nRows = maxCells + 5;

	p.petalPtr = std::make_unique<int2[]>(nRows * nCols);
	for (int i = 0; i < nRows * nCols; i++)	p.petalPtr[i] = { -1, -1 };

	//Create CUDA CBA segment data
	int neighbours[8];
	for (int s = 0; s < nSegments; s++) {
		if (!GetSegmentNeighbours(s + 1, nBands, neighbours)) return msgErrorStr(-99, "Unable to calculate segment neighbours");
		//Neighbours
		p.petalPtr[s] = { neighbours[0], neighbours[1] };
		p.petalPtr[nCols + s] = { neighbours[2], neighbours[3] };
		p.petalPtr[2 * nCols + s] = { neighbours[4], neighbours[5] };
		p.petalPtr[3 * nCols + s] = { neighbours[6], neighbours[7] };
		//Segment ID and cell count
		p.petalPtr[4 * nCols + s] = { int(lookup.at(s + 1).size()), s + 1 };
		//Segment cell's x and y
		for (int c = 0; c < lookup.at(s + 1).size(); c++)
			p.petalPtr[(5 + c) * nCols + s] = { lookup.at(s + 1)[c].x, lookup.at(s + 1)[c].y };
	}

	//Set CUDA CBA search window parameters
	p.petalRows = nRows;
	p.petalCols = nCols;
	p.nPetals = nSegments;
	p.fOffset = srcSize / 2;
	p.firstReduction = uint(powf(2, ceil(log2((float)nCols)))) / 2;
	return 0;
}

