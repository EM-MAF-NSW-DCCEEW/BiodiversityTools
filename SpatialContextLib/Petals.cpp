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

/*
Proposed create search window functions

int CreateSearchWindow(const std::string & paramFN)
int CreateSearchWindow(const std::string & paramFN, const std::string & paramSection)
int CreateSearchWindow(const std::string & paramFN, const std::string & paramSection, CBAParams &p)

int CreatePetals(const std::string & lookupFN, const std::string translateFN)
int CreatePetals(int srcSize, int dstSize, double zoneRatio)

int CreateSegments(int radCells, int nRings, int nSlices, double zoneRatio)
int CreateSegments(std::vector<double> radList)

int CreateSearchWindowFromGrid()
int CreateSearchWindowFromGPUData()

*/


//Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int CreateSearchWindow(const std::string & paramFN)
{
	return CreateSearchWindow(paramFN, "WINDOW");
}

int CreateSearchWindow(const std::string & paramFN, const std::string & paramSection)
{
	CBAParams p;
	return CreateSearchWindow(paramFN, paramSection, p);
}

int CreateSearchWindow(const std::string & paramFN, const std::string & paramSection, CBAParams &p)
{
	std::string lookupFN, translateFN;
	uint srcSize, dstSize;
	double zoneRatio;
	int radCells, nRings;


	std::map<int, std::vector<int2>> lookup;	//Petal ID, vector(pair(cell.x, cell.y))
	std::map<int, int2> translate; //Petal ID, pair(dst.x, dst.y)
	std::vector<double> radList;
	std::unique_ptr<int[]> windowGrid;


	//Read params to work out what we've got
	if (GetParam(paramFN, paramSection, "LOOKUP_FN", lookupFN) && FileExists(lookupFN) &&
		GetParam(paramFN, paramSection, "TRANSLATE_FN", translateFN) && FileExists(translateFN)) {
	
		CreateLookupAndTranslateFromFiles(lookup, translate, srcSize, dstSize, lookupFN, translateFN);
	}
	else if (GetParam(paramFN, paramSection, "SRCSIZE", srcSize) &&
		GetParam(paramFN, paramSection, "DSTSIZE", dstSize) &&
		GetParam(paramFN, paramSection, "ZONERATIO", zoneRatio)) {

		CreateLookupAndTranslateFromParameters(lookup, translate, srcSize, dstSize, zoneRatio);
	}
	else if (GetParam(paramFN, paramSection, "RADIUS", radCells) &&
		GetParam(paramFN, paramSection, "N_RINGS", nRings) &&
		GetParam(paramFN, paramSection, "ZONERATIO", zoneRatio)) {

		CreateSegmentRadListFromParameters(radList, radCells, nRings, zoneRatio);
	}
	else if (GetParam(paramFN, paramSection, "RADIUSLIST", radList)) {
		CreateWindowGridFromSegmentRadList(windowGrid, srcSize, dstSize, radList);
		
	}
	else return msgErrorStr(-1, "for search window " + paramSection, paramFN);
	return 0;
}

//Create lookup & translate from LU & TR files
int CreateLookupAndTranslateFromFiles(std::map<int, std::vector<int2>> & lookup, std::map<int, int2> & translate, uint & srcSize, uint & dstSize, const std::string & lookupFN, const std::string & translateFN)
{
	std::string line, idStr, xStr, yStr;
	int nPetals = 0;

	std::ifstream lookupFS(lookupFN, std::ios::in);
	if (!lookupFS.is_open()) return msgErrorStr(-3, lookupFN);
	while (std::getline(lookupFS, line)) {
		std::stringstream lss(line);
		std::getline(lss, idStr, ',');
		std::getline(lss, xStr, ',');
		std::getline(lss, yStr);
		lookup[std::stoi(idStr)].push_back({ std::stoi(xStr), std::stoi(yStr) });
		nPetals = max(nPetals, std::stoi(idStr));
		srcSize = max(srcSize, (uint) std::stoi(xStr));
		srcSize = max(srcSize, (uint) std::stoi(yStr));
	}

	srcSize += 1;
	if (srcSize % 2 != 1) return msgErrorStr(-99, "Source window size is even", toStr(srcSize));

	//long double sr = std::sqrt(long double(nPetals + 1)); error: expected primary-expression before 'long'
	double sr = std::sqrt(double(nPetals + 1));
	if ((sr - std::floor(sr)) != 0.0L) return msgErrorStr(-99, "Destination window not square", toStr(nPetals + 1));
	dstSize = uint(sr);
	if (dstSize % 2 != 1) return msgErrorStr(-99, "Destination window size is even", toStr(dstSize));

	std::ifstream translateFS(translateFN, std::ios::in);
	if (!translateFS.is_open()) return msgErrorStr(-3, translateFN);
	while (std::getline(translateFS, line)) {
		std::stringstream lss(line);
		std::getline(lss, idStr, ',');
		std::getline(lss, xStr, ',');
		std::getline(lss, yStr);
		translate[std::stoi(idStr)] = { std::stoi(xStr), std::stoi(yStr) };
	}

	return 0;
}

//Create lookup & translate from petal parameters
int CreateLookupAndTranslateFromParameters(std::map<int, std::vector<int2>> & lookup, std::map<int, int2> & translate, const uint & srcSize, const uint & dstSize, const double & zoneRatio)
{
	//Check numeric parameters
	if (zoneRatio < 1.0) return msgErrorStr(-99, "Zone ratio is less than 1.0");
	if (srcSize % 2 != 1 || dstSize % 2 != 1) return msgErrorStr(-99, "Source or destination window size is even");
	if (srcSize < dstSize) return msgErrorStr(-99, "Source window size less than destination window size");

	int x, y;
	int noData = -1;
	int srcCentre = int(srcSize / 2);
	int dstCentre = int(dstSize / 2);

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
		}

		//Right vertical y increasing
		x = srcCentre + zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			y = srcCentre - zoneOffset + int(round(step));
			translate[petalID] = { dstX, dstY++ };
			petalPoints[petalID] = { x, y };
		}

		//Bottom horizontal x decreasing
		y = srcCentre + zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			x = srcCentre + zoneOffset - int(round(step));
			translate[petalID] = { dstX--, dstY };
			petalPoints[petalID] = { x, y };
		}

		//Left vertical y decreasing
		x = srcCentre - zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			y = srcCentre + zoneOffset - int(round(step));
			translate[petalID] = { dstX, dstY-- };
			petalPoints[petalID] = { x, y };
		}
	}

	//Create petal grid searching nearest petal point to allocate petal IDs to cells
	double dist = 0.0f, minDist = 0.0f;
	int pID;
	for (y = 0; y < srcSize; y++) {
		for (x = 0; x < srcSize; x++) {
			if (x != srcCentre || y != srcCentre) {
				minDist = FLT_MAX;
				for (const auto p : petalPoints) {
					if (p.second.x == x && p.second.y == y) {
						pID = p.first;
					}
					else {
						dist = pow(abs(p.second.x - x), 2) + pow(abs(p.second.y - y), 2);
						if (dist < minDist) {
							minDist = dist;
							pID = p.first;
						}
					}
				}
				lookup[pID].push_back({ x, y });
			}
		}
	}
	return 0;
}




//Create window grid from grid file
int CreateWindowGridFromGridFile(std::unique_ptr<int[]> & windowGrid, uint & srcSize, uint & dstSize, const std::string & windowGridFN)
{
	igstream windowGridGS(windowGridFN);
	if (!windowGridGS.is_open()) return msgErrorStr(-3, windowGridFN);
	if (windowGridGS.dataType() != 5) return msgErrorStr(-99, "Data type is not int32", windowGridFN);
	
	uint nCols = windowGridGS.nCols();
	uint nRows = windowGridGS.nRows();
	if (nCols < 1 || nRows < 1) return msgErrorStr(-99, "0 columns or rows in", windowGridFN);
	if (nCols != nRows) return msgErrorStr(-99, "Source window not square " + windowGridFN, "\nCols: " + toStr(nCols), " Rows: " + toStr(nRows));
	if (nCols % 2 != 1) return msgErrorStr(-99, "Source window size is even", toStr(srcSize));
	uint gridSize = nCols * nRows;
	srcSize = uint(nCols);

	windowGrid = std::make_unique<int[]>(gridSize);
	windowGridGS.read((char *)(windowGrid.get()), gridSize * sizeof(int));
	int nPetals = 0;
	for (int i = 0; i < gridSize; i++) nPetals = max(nPetals, windowGrid[i]);

	//long double sr = sqrt(long double(nPetals + 1));	 error: expected primary-expression before 'long'
	double sr = sqrt(double(nPetals + 1));
	if ((sr - std::floor(sr)) != 0.0L) return msgErrorStr(-99, "Destination window not square", toStr(nPetals + 1));
	dstSize = uint(sr);
	if (dstSize % 2 != 1) return msgErrorStr(-99, "Destination window size is even", toStr(dstSize));
	return 0;
}


//Create segment radius list from parameters
int CreateSegmentRadListFromParameters(std::vector<double> & radList, const int & radCells, const int & nRings, const double & zoneRatio)
{
	if (zoneRatio < 1.0f) return msgErrorStr(-99, "Zone ratio is less than 1.0");
	if (radCells < 1) return msgErrorStr(-99, "Radius in cells is less than 1");
	if (nRings < 1) return msgErrorStr(-99, "Number of rings is less than 1");
	if (radCells < nRings) return msgErrorStr(-99, "Radius in cells is less than number of rings");

	//Create ring area list with area within each ring
	double areaFactor = (pow(radCells, 2.0) * M_PI) / ((zoneRatio * (pow(zoneRatio, nRings) - 1.0)) / (zoneRatio - 1.0));
	std::vector<double> areaList;
	areaList.push_back(0.0);
	for (int r = 1; r < nRings + 1; r++) areaList.push_back(pow(zoneRatio, r) * areaFactor);

	//Create ring radius list with outer radius of each ring
	double areaSum = 0.0;
	for (int r = 1; r < nRings + 1; r++) {
		radList.push_back(sqrt((areaList[r] + areaSum) / M_PI));
		areaSum += areaList[r];
	}
	return 0;
}

//Create window grid from segment radius list
int CreateWindowGridFromSegmentRadList(std::unique_ptr<int[]> & windowGrid, uint & srcSize, uint & dstSize, const std::vector<double> & radList)//, const int & radCells, const int & nRings)
{
	//TODO dstSize?

	//Check if first segments will be too small 
	if (radList[0] < sqrt(2)) return msgErrorStr(-99, "First radius is less than SQRT(2) and would result in missing segments.");

	//TODO check nRings and radCells from radList
	int nRings = (int) radList.size();
	int radCells = (int) std::ceil(radList[nRings - 1]);
	int nSlices = 8; //TODO Variable?
	int nSegments = nRings * nSlices;
	
	srcSize = radCells * 2 + 1;
	int halfSrcSize = srcSize / 2;

	//Create window grid
	windowGrid = std::make_unique<int[]>(srcSize * srcSize);
	for (int i = 0; i < srcSize * srcSize; i++) windowGrid[i] = -1;

	//Calc cell distance, degrees and segment ID to fill the lookup table, segment cell count and window grid 
	int segID = -1;
	double cellDist, cellDeg;
	for (int y = -halfSrcSize; y < 1 + halfSrcSize; y++) {
		for (int x = -halfSrcSize; x < 1 + halfSrcSize; x++) {
			segID = -1;
			cellDist = sqrt(pow(double(y), 2) + pow(double(x), 2));
			cellDeg = 180.0 - atan2(double(y), double(x)) * 180.0 / M_PI;

			if (cellDist < 0.5) segID = 0;
			else {
				for (int r = 0; r < radList.size(); r++) {
					if (cellDist <= radList[r]) {
						segID = 1 + (r * nSlices) + int(floor(nSlices * (359 - (int(cellDeg + (180.0 / nSlices)) % 360)) / 360.0));
						break;
					}
				}
			}
			if (segID > nSegments) return msgErrorStr(-99, "Segment ID > nSegments. Error in segment parameters.");
			windowGrid[(y + halfSrcSize) * srcSize + (x + halfSrcSize)] = segID;
		}
	}
	return 0;
}

//Create lookup from window grid
int CreateLookupFromWindowGrid(std::map<int, std::vector<int2>> & lookup, const std::unique_ptr<int[]> & windowGrid, const uint & srcSize) {
	int x, y;
	lookup.clear();
	for (y = 0; y < srcSize; y++) {
		for (x = 0; x < srcSize; x++) {
			if (windowGrid[(y*srcSize) + x] > 0)
				lookup[windowGrid[(y*srcSize) + x]].push_back({ x, y });
		}
	}
	return 0;
}

//Create window grid from lookup
int CreateWindowGridFromLookup(std::unique_ptr<int[]> & windowGrid, const std::map<int, std::vector<int2>> & lookup, const uint & srcSize)
{
	//Make window grid from lookup
	uint srcCentre = uint(srcSize / 2);
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
int CreateGridFileFromWindowGrid(const std::unique_ptr<int[]> & windowGrid, const std::string & windowGridFN, const uint & srcSize) {
	
	ogstream windowGridGS(TifFN(windowGridFN));
	windowGridGS.setTifOptions("LZW", 5, 1);
	windowGridGS.setHeader(srcSize, srcSize, 0.0, 0.0, 1.0, -1.0);
	if (!windowGridGS.is_open()) return msgErrorStr(-4, windowGridFN);
	windowGridGS.write((char *)windowGrid.get(), srcSize * srcSize * (sizeof(int)));
	return 0;
}

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

bool GetSegmentNeighbours(const int segmentID, const int nRings, int *n) {
	try {
		n[0] = segmentID + 8;//n
		n[1] = segmentID - 8;//s
		n[2] = segmentID % 8 == 0 ? segmentID - 7 : segmentID + 1;//e
		n[3] = (segmentID - 1) % 8 == 0 ? segmentID + 7 : segmentID - 1;//w
		n[4] = n[2] + 8;//ne
		n[5] = n[3] + 8;//nw
		n[6] = n[2] - 8;//se
		n[7] = n[3] - 8;//sw
		if (n[0] > 8 * nRings) {
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
int CreateGPUDataFromSegmentsLookup(CBAParams & p, const std::map<int, std::vector<int2>> & lookup, const int & nRings, const uint & srcSize) {
	
	//int nSlices = 8;
	uint nSegments = (uint) lookup.size(); // nRings * nSlices;
	uint nCols = uint(ceil(float(nSegments) / 32.0f) * 32);

	//Allocate CUDA CBA segment data with maxCells rows and init to (-1, -1)
	uint maxCells = 0;
	for (uint s = 0; s < nSegments; s++) maxCells = max(maxCells, (uint) lookup.at(s + 1).size());//(cellCounts[i]);
	uint nRows = maxCells + 5;

	p.petalPtr = std::make_unique<int2[]>(nRows * nCols);
	for (int i = 0; i < nRows * nCols; i++)	p.petalPtr[i] = { -1, -1 };

	//Create CUDA CBA segment data
	int n[8];
	for (int s = 0; s < nSegments; s++) {
		if (!GetSegmentNeighbours(s + 1, nRings, n)) return msgErrorStr(-99, "Unable to calculate segment neighbours");
		//Neighbours
		p.petalPtr[s] = { n[0], n[1] };
		p.petalPtr[nCols + s] = { n[2], n[3] };
		p.petalPtr[2 * nCols + s] = { n[4], n[5] };
		p.petalPtr[3 * nCols + s] = { n[6], n[7] };
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

