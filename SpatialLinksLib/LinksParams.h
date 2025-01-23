/*
LinksParams.h - parameter struct for complete and random sampling Spatial Links implementaion using new custom heap
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

//Notes
//TODO 

#pragma once
#ifndef LINKS_PARAMS_H
#define LINKS_PARAMS_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>

#include "LinksCellNodeGridHeap.h"
#include "FileUtilsLib.h"


using uint = unsigned int;
using ull = unsigned long long;
struct int2 { int x; int y; };
struct PointPair { int2 p1, p2; };

const double sqrt2{ std::sqrt(2.0) };
const double halfSqrt2 = sqrt2 * 0.5;

struct LinksParams {

	//For similarity calcs
	int nTxGrids;
	std::vector<igstream> txSrcGSs;
	std::vector<igstream> txDstGSs;
	
	
	//Grid streams
	igstream habGS, cstGS;
	ogstream lnkGS, cxtGS, pntGS;

	//Complete Sampling Alt Dst Hab
	igstream srcHabGS, dstHabGS;

	//Cost grid header parameters
	uint nCols, nRows;
	ull nCells;
	double xll, yll;
	double cellSize;
	float noData;

	bool useHabGrid;
	bool useMaxED;
	bool doLinks, doContext;
	bool useSrcDstHabWeight;
	bool useFCellHabWeight;

	//Complete sampling search parameters
	double searchRadius;
	float habMin;
	double maxED;
	
	//Random sampling search parameters
	std::string habGFN, cstGFN, lnkGFN, pntGFN;
	std::string  pointPairsTabFN;
	bool havePointPairsTabFN;
	bool pointPairsTabExists;
	bool writePointPairsTab;
	bool writePointGrid;
	double maxPairDist;
	int nPointPairs;
	int nPaths;
	int writeInterval;
	float P;
	std::vector<PointPair> pointPairs;

	double minPairDistSqrd, maxPairDistSqrd;
		
	double d_i, d_a;

	double lnkExponent;

	//Multithreading
	bool useMTProcessing;
	int nThreads;


};

//Distance factor for orthagonal and diagonal neighbours
const double neighbourDists[8]{ 
	0.5, 0.5, 0.5, 0.5,//Orthogonal
	halfSqrt2, halfSqrt2, halfSqrt2, halfSqrt2//Diagonal
};

//x and y offsets of neighbouring pixels
const int2 neighbourDirs[8]{
	//{ x, y }
	{ 0,-1 },//N
	{ 1,0 },//E
	{ 0,1 },//S
	{ -1,0 },//W
	{ 1,-1 },//NE
	{ 1,1 },//SE
	{ -1,1 },//SW
	{ -1,-1 }//NW
};


#endif