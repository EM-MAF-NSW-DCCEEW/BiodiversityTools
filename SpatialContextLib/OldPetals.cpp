/*
OldPetals.cpp - Legacy functions for creating CBA petal files and data structures
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

//Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//TODO no longer need to create or use petal lookup and translate files.
//TODO replace ints with uints where appropriate

int CreatePetals(const std::string & baseFN, int srcSize, int dstSize, float zoneRatio, bool createFiles, CBAParams &p)
{
	//Check numeric parameters
	if (zoneRatio < 1.0f || srcSize % 2 != 1 || dstSize % 2 != 1 || srcSize < dstSize)
		return msgErrorStr(-99, "Source or destination size is even or zone ratio is less than 1.0");

	int i, x, y;
	int noData = -1;
	int srcCentre = int(srcSize / 2);
	int dstCentre = int(dstSize / 2);

	//Create zonefactors
	float rawSum = 0.0f, zoneSize = 0.0f, scaledSum = 0.0f;
	std::map<int, float>zoneFactors;//Zone, factor
	for (int zoneNo = 0; zoneNo < dstCentre; zoneNo++)	rawSum += pow(zoneRatio, zoneNo);
	for (int zoneNo = 1; zoneNo < dstCentre + 1; zoneNo++) {
		zoneSize = std::fmaxf(srcCentre / rawSum * pow(zoneRatio, zoneNo - 1), 1.0f);
		scaledSum += (zoneNo > 1 && zoneSize + scaledSum > srcCentre) ? srcCentre - scaledSum : zoneSize;
		zoneFactors[zoneNo] = std::fmax(scaledSum / zoneNo, 1.0f);
		if (zoneNo > 1 && zoneFactors[zoneNo] * zoneNo <= zoneFactors[zoneNo - 1] * (zoneNo - 1)) return msgErrorStr(-99, "Unable to calculate zone factors");
	}

	//Set a unique petalID for each petal starting from 1
	int petalID = 1;
	int zone, zoneOffset, dstY, dstX;
	float factor, prevOffset, splitOffset, inc, step;
	std::map<int, std::pair<int, int>> petalPoints; //Petal ID, pair(src.y, src.x)
	std::map<int, std::pair<int, int>> translate; //Petal ID, pair(dst.x, dst.y)
	translate[0] = std::make_pair(dstCentre, dstCentre);
	std::unique_ptr<int[]> pointGrid = std::make_unique<int[]>(srcSize * srcSize);
	for (i = 0; i < srcSize * srcSize; i++) pointGrid[i] = noData;

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
			translate[petalID] = std::make_pair(dstX++, dstY);
			petalPoints[petalID] = std::make_pair(x, y);
			pointGrid[(y * srcSize) + x] = petalID++;
		}

		//Right vertical y increasing
		x = srcCentre + zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			y = srcCentre - zoneOffset + int(round(step));
			translate[petalID] = std::make_pair(dstX, dstY++);
			petalPoints[petalID] = std::make_pair(x, y);
			pointGrid[(y * srcSize) + x] = petalID++;
		}

		//Bottom horizontal x decreasing
		y = srcCentre + zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			x = srcCentre + zoneOffset - int(round(step));
			translate[petalID] = std::make_pair(dstX--, dstY);
			petalPoints[petalID] = std::make_pair(x, y);
			pointGrid[(y * srcSize) + x] = petalID++;
		}

		//Left vertical y decreasing
		x = srcCentre - zoneOffset;
		for (step = 0; int(round(step)) < zoneOffset * 2; step += inc) {
			y = srcCentre + zoneOffset - int(round(step));
			translate[petalID] = std::make_pair(dstX, dstY--);
			petalPoints[petalID] = std::make_pair(x, y);
			pointGrid[(y * srcSize) + x] = petalID++;
		}
	}

	//Create petal grid searching nearest petal point to allocate petal IDs to cells
	double dist = 0.0f;
	double minDist = 0.0f;
	std::unique_ptr<int[]> petalGrid = std::make_unique<int[]>(srcSize * srcSize);
	for (y = 0; y < srcSize; y++) {
		for (x = 0; x < srcSize; x++) {
			petalGrid[y * srcSize + x] = noData;
			minDist = FLT_MAX;
			for (const auto p : petalPoints) {
				if (p.second.first == x && p.second.second == y) {
					petalGrid[y * srcSize + x] = p.first;
					break;
				}
				else {
					dist = pow(abs(p.second.first - x), 2) + pow(abs(p.second.second - y), 2);
					if (dist < minDist) {
						minDist = dist;
						petalGrid[y * srcSize + x] = p.first;
					}
				}
			}
		}
	}
	petalGrid[srcCentre * srcSize + srcCentre] = 0;

	//Create petal lookup table
	std::map<int, std::vector<std::pair<int, int>>> lookup;	//Petal ID, vector(pair(cell.x, cell.y))
	for (y = 0; y < srcSize; y++) {
		for (x = 0; x < srcSize; x++) {
			if (petalGrid[(y*srcSize) + x] != 0)
				lookup[petalGrid[(y*srcSize) + x]].push_back(std::make_pair(x, y));
		}
	}

	//Set search window parameters and initialise petal data storage
	int maxCells = 1, cellsStart = 5, petalsRow;
	for (auto p : lookup) maxCells = max(maxCells, int(p.second.size()));
	p.fOffset = srcSize / 2;
	p.nPetals = dstSize * dstSize - 1;
	p.petalCols = int(ceil(float(p.nPetals) / 32.0f) * 32);
	p.petalRows = maxCells + cellsStart;
	p.firstReduction = int(powf(2, ceil(log2(float(p.petalCols))))) / 2;
	p.petalPtr = std::make_unique<int2[]>(p.petalCols * p.petalRows);
	for (i = 0; i < p.petalCols * p.petalRows; i++)	p.petalPtr[i] = { -1, -1 };

	//Fill petals data table
	for (auto pA : translate) {
		if (pA.first == 0) continue;
		for (auto pB : translate) {
			if (pB.first == 0) continue;
			//First four rows contain neighbouring petal indicies
			if (pB.second.first == pA.second.first && pB.second.second == pA.second.second - 1) p.petalPtr[pA.first - 1].x = pB.first - 1;
			if (pB.second.first == pA.second.first && pB.second.second == pA.second.second + 1) p.petalPtr[pA.first - 1].y = pB.first - 1;
			if (pB.second.first == pA.second.first - 1 && pB.second.second == pA.second.second) p.petalPtr[p.petalCols + pA.first - 1].x = pB.first - 1;
			if (pB.second.first == pA.second.first + 1 && pB.second.second == pA.second.second) p.petalPtr[p.petalCols + pA.first - 1].y = pB.first - 1;
			if (pB.second.first == pA.second.first - 1 && pB.second.second == pA.second.second - 1) p.petalPtr[2 * p.petalCols + pA.first - 1].x = pB.first - 1;
			if (pB.second.first == pA.second.first + 1 && pB.second.second == pA.second.second - 1) p.petalPtr[2 * p.petalCols + pA.first - 1].y = pB.first - 1;
			if (pB.second.first == pA.second.first - 1 && pB.second.second == pA.second.second + 1) p.petalPtr[3 * p.petalCols + pA.first - 1].x = pB.first - 1;
			if (pB.second.first == pA.second.first + 1 && pB.second.second == pA.second.second + 1) p.petalPtr[3 * p.petalCols + pA.first - 1].y = pB.first - 1;
		}

		//From cellStart contains coords of cells in each petal
		petalsRow = cellsStart;
		for (auto lu : lookup[pA.first]) {
			p.petalPtr[petalsRow * p.petalCols + pA.first - 1].x = lu.first;
			p.petalPtr[petalsRow * p.petalCols + pA.first - 1].y = lu.second;
			petalsRow++;
		}

		//Fifth row contains number of cells in each petal and petal ID
		p.petalPtr[4 * p.petalCols + pA.first - 1].x = int(lookup[pA.first].size());
		p.petalPtr[4 * p.petalCols + pA.first - 1].y = pA.first - 1;
	}

	if (createFiles) {
		//Write zoneFactor multipliers to file
		std::string zonesFN(baseFN + "_zoneFactors.txt");
		std::ofstream zonesFS(zonesFN, std::ios::out);
		if (!zonesFS.is_open()) return msgErrorStr(-4, zonesFN);
		for (auto z : zoneFactors) zonesFS << z.first << "," << z.second << "\n";

		//Write point and petal grids to ascii files
		std::string pointGridFN(baseFN + "_points.asc");
		std::string petalGridFN(baseFN + "_petals.asc");
		if (!WriteAsciiHeader(pointGridFN, srcSize, srcSize, 0.0, 0.0, 1, float(noData))) return msgErrorStr(-4, pointGridFN);
		if (!WriteAsciiHeader(petalGridFN, srcSize, srcSize, 0.0, 0.0, 1, float(noData))) return msgErrorStr(-4, petalGridFN);
		std::ofstream pointGridFS(pointGridFN, std::ios::app);
		std::ofstream petalGridFS(petalGridFN, std::ios::app);
		if (!pointGridFS.is_open()) return msgErrorStr(-4, pointGridFN);
		if (!petalGridFS.is_open()) return msgErrorStr(-4, petalGridFN);
		for (y = 0; y < srcSize; y++) {
			for (int x = 0; x < srcSize; x++) {
				pointGridFS << pointGrid[(y * srcSize) + x] << (x + 1 < srcSize ? "\t" : "\n");
				petalGridFS << petalGrid[(y * srcSize) + x] << (x + 1 < srcSize ? "\t" : "\n");
			}
		}

		//Write out lookup table as petal id, cell x, cell y and translate table as petal id, petal x, petal y
		std::string lookupFN(baseFN + "_lookup.txt");
		std::string translateFN(baseFN + "_translate.txt");
		std::ofstream luFS(lookupFN, std::ios::out);
		std::ofstream trFS(translateFN, std::ios::out);
		if (!luFS.is_open()) return msgErrorStr(-4, lookupFN);
		if (!trFS.is_open()) return msgErrorStr(-4, translateFN);
		for (auto petal : lookup)
			for (auto cell : petal.second)
				luFS << petal.first << "," << cell.first << "," << cell.second << "\n";
		for (auto petal : translate)
			trFS << petal.first << "," << petal.second.first << "," << petal.second.second << "\n";

		//Write petals data to file
		std::string dataFN(baseFN + "_cuda_petals.txt");
		std::ofstream dataFS(dataFN, std::ios::out);
		if (!dataFS.is_open()) return msgErrorStr(-4, dataFN);
		for (uint j = 0; j < p.petalRows; j++)
			for (uint k = 0; k < p.petalCols; k++)
				dataFS << "" << p.petalPtr[j * p.petalCols + k].x << "," << p.petalPtr[j * p.petalCols + k].y << (k + 1 < p.petalCols ? "\t" : "\n");
	}
	return 0;
}

//Create p.petal data for CUDA CBA from petal lookup and translate files
bool CreatePetalData(const char *translateFN, const char *lookupFN, const char *baseFN, CBAParams &p)
{
	return CreatePetalData(std::string(translateFN), std::string(lookupFN), std::string(baseFN), p);
}
//bool CreatePetalData(const std::string & petalsTRFN, const std::string & petalsLUFN, const std::string  & petalsFN, CBAParams &p)
bool CreatePetalData(const std::string & translateFN, const std::string & lookupFN, const std::string  & baseFN, CBAParams &p)
{
	int cellsStart = 5;
	int	trRows = 0, luRows = 0;
	int petalID = 0;
	int maxCells = 1;
	int srcDim = 0;

	std::string is;
	std::ifstream trfs(translateFN, std::ios::in);
	std::ifstream lufs(lookupFN, std::ios::in);
	if (!trfs.is_open() || !lufs.is_open())	return false;

	//Get tr and lu rows from files
	std::getline(trfs, is);//Skip focal petal
	while (std::getline(trfs, is)) ++trRows;
	while (std::getline(lufs, is)) ++luRows;
	trfs.clear(); trfs.seekg(0);
	lufs.clear(); lufs.seekg(0);

	//Create petal arrays
	std::unique_ptr<int4[]> trA = std::make_unique<int4[]>(trRows);
	std::unique_ptr<int3[]> luA = std::make_unique<int3[]>(luRows);

	//Read tr file
	std::getline(trfs, is);//Skip focal petal
	for (int i = 0; i < trRows; i++) {
		std::getline(trfs, is);
		std::stringstream iss(is);

		//Get petal ID
		std::getline(iss, is, ',');
		trA[i].z = std::stoi(is);

		//Get petal x and y coord, set cell count to 0
		std::getline(iss, is, ',');
		trA[i].x = std::stoi(is);
		std::getline(iss, is);
		trA[i].y = std::stoi(is);
		trA[i].w = 0;
	}

	//Read lu file
	for (int i = 0; i < luRows; i++) {
		std::getline(lufs, is);
		std::stringstream iss(is);

		//Get cell's petal ID, increment cell count, and set maxCells
		std::getline(iss, is, ',');
		luA[i].z = std::stoi(is);
		maxCells = max(++(trA[luA[i].z - 1].w), maxCells);

		//Get cell x and y coords and max srcDim
		std::getline(iss, is, ',');
		luA[i].x = std::stoi(is);
		std::getline(iss, is);
		luA[i].y = std::stoi(is);
		srcDim = max(luA[i].x, srcDim);
	}

	//Close files
	trfs.close();
	lufs.close();

	//Set search window parameters
	p.fOffset = srcDim / 2;
	p.nPetals = trRows;
	p.petalCols = int(ceil(float(p.nPetals) / 32.0f) * 32);
	p.petalRows = maxCells + cellsStart;
	p.firstReduction = int(powf(2, ceil(log2(float(p.petalCols))))) / 2;

	//Init petals storage
	p.petalPtr = std::make_unique<int2[]>(p.petalCols * p.petalRows);
	for (uint i = 0; i < p.petalCols * p.petalRows; i++)	p.petalPtr[i] = { -1, -1 };

	//Fill petals table, written late at night...
	int trRowA = 0, trRowB = 0, luRow = 0, petalsRow;
	for (trRowA = 0; trRowA < trRows; trRowA++) {

		//First four rows contain neighbouring petal indicies
		for (trRowB = 0; trRowB < trRows; trRowB++) {
			if (trA[trRowB].x == trA[trRowA].x && trA[trRowB].y == trA[trRowA].y - 1) p.petalPtr[trRowA].x = trRowB;
			if (trA[trRowB].x == trA[trRowA].x && trA[trRowB].y == trA[trRowA].y + 1) p.petalPtr[trRowA].y = trRowB;
			if (trA[trRowB].x == trA[trRowA].x - 1 && trA[trRowB].y == trA[trRowA].y) p.petalPtr[p.petalCols + trRowA].x = trRowB;
			if (trA[trRowB].x == trA[trRowA].x + 1 && trA[trRowB].y == trA[trRowA].y) p.petalPtr[p.petalCols + trRowA].y = trRowB;
			if (trA[trRowB].x == trA[trRowA].x - 1 && trA[trRowB].y == trA[trRowA].y - 1) p.petalPtr[2 * p.petalCols + trRowA].x = trRowB;
			if (trA[trRowB].x == trA[trRowA].x + 1 && trA[trRowB].y == trA[trRowA].y - 1) p.petalPtr[2 * p.petalCols + trRowA].y = trRowB;
			if (trA[trRowB].x == trA[trRowA].x - 1 && trA[trRowB].y == trA[trRowA].y + 1) p.petalPtr[3 * p.petalCols + trRowA].x = trRowB;
			if (trA[trRowB].x == trA[trRowA].x + 1 && trA[trRowB].y == trA[trRowA].y + 1) p.petalPtr[3 * p.petalCols + trRowA].y = trRowB;
		}

		//Coords of cells in each petal
		petalsRow = cellsStart;
		while (luA[luRow].z == trA[trRowA].z) {
			p.petalPtr[petalsRow * p.petalCols + trRowA].x = luA[luRow].x;
			p.petalPtr[petalsRow * p.petalCols + trRowA].y = luA[luRow].y;
			luRow++;
			petalsRow++;
		}

		//Number of cells in each petal
		p.petalPtr[4 * p.petalCols + trRowA].x = petalsRow - cellsStart;
		p.petalPtr[4 * p.petalCols + trRowA].y = trRowA;
	}

	//Write petals to a file
	std::string dataFN(baseFN + "_data.txt");
	std::ofstream dataFS(dataFN, std::ios::out);
	if (!dataFS.is_open()) return false;
	for (uint j = 0; j < p.petalRows; j++)
		for (uint k = 0; k < p.petalCols; k++)
			dataFS << "" << p.petalPtr[j * p.petalCols + k].x << "," << p.petalPtr[j * p.petalCols + k].y << (k + 1 < p.petalCols ? "\t" : "\n");

	return true;
}

//Gets the segments neighbour segment IDs
//(n, s, e, w, ne, nw, se, sw), n = outwards
bool GetSegmentNeighbours(const int segmentID, const int nSlices, const int nRings, int *n) {
	try {
		n[0] = segmentID + 8;//n
		n[1] = segmentID - 8;//s
		n[2] = segmentID % nSlices == 0 ? segmentID - 7 : segmentID + 1;//e
		n[3] = (segmentID - 1) % nSlices == 0 ? segmentID + 7 : segmentID - 1;//w
		n[4] = n[2] + 8;//ne
		n[5] = n[3] + 8;//nw
		n[6] = n[2] - 8;//se
		n[7] = n[3] - 8;//sw
		if (n[0] > nSlices * nRings) {
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

//Create p.petalData for CUDA CBA from radial segments
int CreateSegments(int radCells, int nRings, int nSlices, double zoneRatio, const std::string & outFN, CBAParams &params) {

	//Create ring area list with area within each ring
	double areaFactor = (pow(radCells, 2.0) * M_PI) / ((zoneRatio * (pow(zoneRatio, nRings) - 1.0)) / (zoneRatio - 1.0));
	std::vector<double> areaList;
	areaList.push_back(0.0);
	for (int r = 1; r < nRings + 1; r++) areaList.push_back(pow(zoneRatio, r) * areaFactor);

	//Create ring radius list with outer radius of each ring
	std::vector<double> radList;
	double areaSum = 0.0;
	for (int r = 1; r < nRings + 1; r++) {
		radList.push_back(sqrt((areaList[r] + areaSum) / M_PI));
		areaSum += areaList[r];
	}

	return CreateSegments(radCells, nRings, nSlices, radList, outFN, params);
}

int CreateSegments(int radCells, int nRings, int nSlices, std::vector<double> radList, const std::string & outFN, CBAParams &p) {

	//Check if first segments will be too small 
	if (radList[0] < sqrt(2)) return msgErrorStr(-99, "First radius is less than SQRT(2) and would result in missing segments.");
	//Search window parameters
	std::string gridFN(outFN + "_segments.asc");
	std::string tableFN(outFN + "_segments.txt");

	int nSegments = nRings * nSlices;
	int nCols = int(ceil(float(nSegments) / 32.0f) * 32);
	int windowDim = radCells * 2 + 1;
	int halfWinDim = windowDim / 2;
	int cellsStart = 5;

	std::vector<int> depthList;
	depthList.push_back(int(radList[0]));
	for (int r = 1; r < nRings; r++) {
		depthList.push_back(int(round(radList[r] - radList[r - 1])));
	}

	//Create the lookup table, segment cell count and window grid
	std::map<int, std::vector<std::pair<int, int>>> lookup;
	std::map<int, int> cellCounts;
	std::unique_ptr<int[]> windowGrid = std::make_unique<int[]>(windowDim * windowDim);
	for (int i = 1; i < nSegments + 1; i++)	cellCounts[i] = 0;
	for (int i = 0; i < windowDim * windowDim; i++) windowGrid[i] = -1;

	//Calc cell distance, degrees and segment ID to fill the lookup table, segment cell count and window grid 
	int segID = -1;
	double cellDist, cellDeg;
	for (int y = -halfWinDim; y < 1 + halfWinDim; y++) {
		for (int x = -halfWinDim; x < 1 + halfWinDim; x++) {
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
				if (segID != -1) {
					lookup[segID].push_back(std::make_pair(x + halfWinDim, y + halfWinDim));
					cellCounts[segID] += 1;
				}
			}
			if (segID > nSegments) return msgErrorStr(-99, "Segment ID > nSegments. Error in segment parameters.");
			windowGrid[(y + halfWinDim) * windowDim + (x + halfWinDim)] = segID;
		}
	}
	//Check for empty segments
	for (int k = 1; k <= nSegments; k++)
		if (lookup.count(k) == 0) return msgErrorStr(-99, "Segment lookup contains no cells for segment", toStr(k));

	//Write windowGrid to file
	if (!WriteAsciiHeader(gridFN, windowDim, windowDim, 0.0, 0.0, 1, -1)) return msgErrorStr(-4, gridFN);
	std::ofstream windowGridFS(gridFN, std::ios::app);
	if (!windowGridFS.is_open()) return msgErrorStr(-4, gridFN);
	for (int y = 0; y < windowDim; y++)
		for (int x = 0; x < windowDim; x++)
			windowGridFS << windowGrid[(y * windowDim) + x] << (x + 1 < windowDim ? "\t" : "\n");
	windowGrid.reset();

	//Allocate CUDA CBA segment data with maxCells rows and init to (-1, -1)
	int maxCells = 0;
	for (int i = 1; i < nSegments + 1; i++)	maxCells = max(maxCells, cellCounts[i]);
	int nRows = maxCells + cellsStart;

	p.petalPtr = std::make_unique<int2[]>(nRows * nCols);
	for (int i = 0; i < nRows * nCols; i++)	p.petalPtr[i] = { -1, -1 };

	//Create CUDA CBA segment data
	int n[8];
	for (int i = 0; i < nSegments; i++) {
		if (!GetSegmentNeighbours(i + 1, nSlices, nRings, n)) return msgErrorStr(-99, "Unable to calculate segment neighbours");
		//Neighbours
		p.petalPtr[i] = { n[0], n[1] };
		p.petalPtr[nCols + i] = { n[2], n[3] };
		p.petalPtr[2 * nCols + i] = { n[4], n[5] };
		p.petalPtr[3 * nCols + i] = { n[6], n[7] };
		//Segment ID and cell count
		p.petalPtr[4 * nCols + i] = { cellCounts[i + 1], depthList[i / 8] };
		//Segment cell's x and y
		for (int c = 0; c < cellCounts[i + 1]; c++)
			p.petalPtr[(5 + c) * nCols + i] = { lookup[i + 1][c].first, lookup[i + 1][c].second };
	}

	//Write CUDA CBA segment data to file
	std::ofstream segmentFS(tableFN, std::ios::out);
	if (!segmentFS.is_open()) return msgErrorStr(-4, tableFN);
	for (int i = 0; i < nRows; i++) {
		for (int k = 0; k < nCols; k++)
			segmentFS << p.petalPtr[i * nCols + k].x << "," << p.petalPtr[i * nCols + k].y << "\t";
		segmentFS << std::endl;
	}

	//Set CUDA CBA search window parameters
	p.petalRows = nRows;
	p.petalCols = nCols;
	p.nPetals = nSegments;
	p.fOffset = halfWinDim;
	p.firstReduction = int(powf(2, ceil(log2((float)nCols)))) / 2;
	return 0;
}

//Get petal files from parameter file and section
int GetPetalFiles(const std::string & paramFN, const std::string & paramSection, const std::string workDir, std::string & lookupFN, std::string & translateFN, CBAParams & params) {
	if (GetParam(paramFN, paramSection, "LOOKUP_FN", lookupFN) && FileExists(lookupFN) &&
		GetParam(paramFN, paramSection, "TRANSLATE_FN", translateFN) && FileExists(translateFN)) {
		return 0;
	}
	else {
		int winSrcSize, winDstSize;
		float winZoneRatio;
		bool saveTmpFiles;

		//Get window parameters
		if (!GetParam(paramFN, paramSection, "SRCSIZE", winSrcSize)) return msgErrorStr(-1, "SRCSIZE", paramFN);
		if (!GetParam(paramFN, paramSection, "DSTSIZE", winDstSize)) return msgErrorStr(-1, "DSTSIZE", paramFN);
		if (!GetParam(paramFN, paramSection, "ZONERATIO", winZoneRatio)) return msgErrorStr(-1, "ZONERATIO", paramFN);
		if (!GetParam(paramFN, paramSection, "SAVETMPFILES", saveTmpFiles)) saveTmpFiles = true;
		std::string petalsFN = workDir + pathSep + "petals_" + toStr(winSrcSize) + "_" + toStr(winDstSize) + "_" + toStr(int(round(winZoneRatio * 10)));

		//Run CreatePetals to create lookup and translate
		int result = CreatePetals(petalsFN, winSrcSize, winDstSize, winZoneRatio, saveTmpFiles, params);
		if (result != 0) return msgErrorStr(-15, "CreatePetals()", toStr(result));

		//Set petal lookup and translate filenames in source parameter file
		lookupFN = petalsFN + "_lookup.txt";
		translateFN = petalsFN + "_translate.txt";
		if (!SetParam(paramFN, paramSection, "LOOKUP_FN", lookupFN))return msgErrorStr(-14, "LOOKUP_FN", paramFN);
		if (!SetParam(paramFN, paramSection, "TRANSLATE_FN", translateFN)) return msgErrorStr(-14, "TRANSLATE_FN", paramFN);
		msgText(("Created petal files for " + toStr(winSrcSize) + " cells to " + toStr(winDstSize) + " petals with a zone ratio of " + toStr(winZoneRatio)).c_str());
		return 0;
	}
	return -1;
}

int GetPetalFiles(const std::string & paramFN, const std::string & paramSection, std::string & lookupFN, std::string & translateFN, CBAParams & params) {
	std::string workDir;
	if (!GetFilesFolder(paramFN, workDir)) return msgErrorStr(-99, "Unable to get directory from ", paramFN);;
	return GetPetalFiles(paramFN, paramSection, workDir, lookupFN, translateFN, params);
}

int ReadLUTRFiles(const std::string & lookupFN, const std::string & translateFN, std::unique_ptr<int[]> & windowGrid) {

	std::map<int, std::vector<std::pair<int, int>>> lookup;	//Petal ID, vector(pair(cell.x, cell.y))
	std::map<int, std::pair<int, int>> translate; //Petal ID, pair(dst.x, dst.y)

	std::ifstream lookupFS, translateFS;
	uint luRows = 0, trRows = 0;

	std::string line, idStr, xStr, yStr;

	uint nPetals = 0, maxCellXY = 0;

	lookupFS.open(lookupFN, std::ios::in);
	if (!lookupFS.is_open()) return msgErrorStr(-3, lookupFN);

	while (std::getline(lookupFS, line)) {
		std::stringstream lss(line);
		std::getline(lss, idStr, ',');
		std::getline(lss, xStr, ',');
		std::getline(lss, yStr);
		lookup[std::stoi(idStr)].push_back(std::make_pair(std::stoi(xStr), std::stoi(yStr)));
		luRows++;
		maxCellXY = max(maxCellXY, uint(std::stoi(xStr)));
		maxCellXY = max(maxCellXY, uint(std::stoi(yStr)));
	}

	translateFS.open(translateFN, std::ios::in);
	if (!translateFS.is_open()) return msgErrorStr(-3, translateFN);

	int id = 0, maxCells = 0;

	while (std::getline(translateFS, line)) {
		std::stringstream lss(line);
		std::getline(lss, idStr, ',');
		std::getline(lss, xStr, ',');
		std::getline(lss, yStr);
		translate[std::stoi(idStr)] = std::make_pair(std::stoi(xStr), std::stoi(yStr));
		trRows++;
		nPetals = max(nPetals, uint(std::stoi(idStr)));//no stou
	}

	//Make window grid from lookup
	uint srcSize = maxCellXY + 1;
	uint dstSize = uint(round(sqrt(double(nPetals + 1))));
	uint srcCentre = uint(srcSize / 2);
	uint dstCentre = uint(dstSize / 2);

	windowGrid = std::make_unique<int[]>(srcSize * srcSize);
	for (int c = 0; c < srcSize * srcSize; c++) windowGrid[c] = -1;
	windowGrid[srcCentre * srcSize + srcCentre] = 0;
	for (const auto p : lookup) {
		for (const auto c : p.second) {
			windowGrid[c.second * srcSize + c.first] = p.first;
		}
	}
	return 0;
}

