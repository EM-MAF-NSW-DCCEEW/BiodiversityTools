/*
LinksRandomSampling.cpp - random sampling Spatial Links implementaion using new custom heap
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
//Standardise parameter file format
//CreatePointPairs writes point grid
//PointTab is either input or output parameter
//Is it time to save


#include "LinksRandomSampling.h"
#include <sstream>
#include <ctime>
#include <random>

int LinksRandomSampling(std::string paramFN, LinksParams &p) {


	//Cost grid FN is required
	if (!GetParam(paramFN, "INPUTS", "CostGrid", p.cstGFN)) msgErrorStr(-1, "CostGrid", paramFN);
	//Habitat grid FN is optional
	p.useHabGrid = GetParam(paramFN, "INPUTS", "HabGrid", p.habGFN);
	//Point pairs table FN is optional
	p.havePointPairsTabFN = GetParam(paramFN, "INPUTS", "PointTab", p.pointPairsTabFN);
	//Links grid FN is required
	if (!GetParam(paramFN, "OUTPUTS", "LinksGrid", p.lnkGFN)) msgErrorStr(-1, "LinksGrid", paramFN);
	//Point grid FN is optional
	p.writePointGrid = GetParam(paramFN, "OUTPUTS", "PointGrid", p.pntGFN);
	//Maximum pair distance is required
	if (!GetParam(paramFN, "SEARCH", "MaxPairDist", p.maxPairDist)) msgErrorStr(-1, "MaxPairDist", paramFN);
	//Minimum habitat threshold is requred
	if (!GetParam(paramFN, "SEARCH", "MinHabThreshold", p.habMin)) msgErrorStr(-1, "MinHabThreshold", paramFN);
	//Maximum effective distance is optional
	p.useMaxED = GetParam(paramFN, "SEARCH", "MaxED", p.maxED);
	//Number of point pairs is required
	if (!GetParam(paramFN, "SEARCH", "nPointPairs", p.nPointPairs)) msgErrorStr(-1, "nPointPairs", paramFN);
	//Probability factor is optional, default is 1.0
	if (!GetParam(paramFN, "SEARCH", "ProbabilityFactor", p.P)) p.P = 1.0;

	p.maxPairDistSqrd = p.maxPairDist * p.maxPairDist;
	p.minPairDistSqrd = 0;

	//Decay I is requried 
	if (!GetParam(paramFN, "DECAY", "I", p.d_i)) msgErrorStr(-1, "I", paramFN);
	//Decay A is required
	if (!GetParam(paramFN, "DECAY", "A", p.d_a)) msgErrorStr(-1, "A", paramFN);
	//Links exponent is optional, default is 0.25
	if (!GetParam(paramFN, "OPTIONS", "LinksExponent", p.lnkExponent)) p.lnkExponent = 0.25;
	//Write interval is optional, default is 0 
	if (!GetParam(paramFN, "OPTIONS", "WriteInterval", p.writeInterval)) p.writeInterval = 0;
	//Number of threads is optional, default is single threaded
	p.useMTProcessing = GetParam(paramFN, "OPTIONS", "nthreads", p.nThreads);

	//Open cost grid
	p.cstGS.open(p.cstGFN);
	if (!p.cstGS.is_open()) return msgErrorStr(-3, p.cstGFN);
	p.cstGS.getHeader(p.nCols, p.nRows, p.xll, p.yll, p.cellSize, p.noData);
	p.nCells = p.nCols * p.nRows;

	//Open hab grid if using
	if (p.useHabGrid) {
		p.habGS.open(p.habGFN);
		if (!p.habGS.is_open()) return msgErrorStr(-3, p.habGFN);
		if (!p.habGS.compareHeader(p.cstGS)) return msgErrorStr(-7, p.habGFN, p.cstGFN);
	}

	//Read or create point pairs
	//if pointsPairTab exists try loading point pairs
	int retval = 0;
	if (p.havePointPairsTabFN && FileExists(p.pointPairsTabFN)) {
		retval = LoadPointPairs(p);
	}
	else if (p.useHabGrid) {
		retval = CreatePointPairs(p);
	}
	else {
		return msgErrorStr(-99, "Need either existing PointTab to load, or HabGrid to create point pairs");
	}
	//TODO test retval


	//p.pntGS.close();
	p.habGS.open(p.habGFN);

	int result = 0;
	if (p.useMTProcessing) {
		//TODO multithreaded implementation;
		msgErrorStr(-99, "TODO multithreaded implementation");
	}
	else
	{
		result = LinksRandomSampling_ST(p);
	}
	return result;
}

int LinksRandomSampling_ST(LinksParams &p) {
	msgText("Performing single-threaded Spatial Links analysis with random grid sampling");

	//Path Loop Variables
	int saveEvery = 1;//TODO
	int nextSave = 0;
	int pathsCompleted = 0;// num completed paths
	int pathsFailed = 0;// num failed paths
	
	double nAccumCost = 0.0;
	double srcDstHab = 1.0;
	double d_e;
	float pathValue;
	int n, nx, ny;


	//Initialise CellNodes
	msgText("Allocating memory");
	CellNodeGrid nodes((size_t)p.nRows, (size_t)p.nCols);
	CellNodeHeap minHeap(p.nCells);
	CellNode *nCell, *cCell, *pCell, *srcCell, *dstCell;

	//Init data
	std::vector<std::vector<float>> cstData(p.nRows, std::vector<float>(p.nCols, p.noData));
	std::vector<std::vector<float>> habData(p.nRows, std::vector<float>(p.nCols, p.noData));
	std::vector<std::vector<float>> lnkData(p.nRows, std::vector<float>(p.nCols, 0.0f));
	
	msgText("Reading input data to memory");
	for (uint r = 0; r < p.nRows; r++) {
		p.cstGS.read((char *)(cstData[r].data()), sizeof(float) * p.nCols);
	}

	if (p.useHabGrid) {
		for (uint r = 0; r < p.nRows; r++) {
			p.habGS.read((char *)(habData[r].data()), sizeof(float) * p.nCols);
		}
	}

	msgText("Processing least cost paths");
	for (int i = 0; i < p.pointPairs.size(); i++) {
		msgProgress("Percent complete: ", i * 100 / p.pointPairs.size());
		//msgString("Path:\t" + toStr(i) + "\tSrc:\t " + toStr(srcCell->x) + ", " + toStr(srcCell->y) + "\tDst:\t " + toStr(dstCell->x) + ", " + toStr(dstCell->y) + "\n");

		//Reset nodes, heapSize and cell pointers
		nAccumCost = 0.0;
		nodes.reset();
		minHeap.reset();

		cCell = nullptr;//Current search cell
		nCell = nullptr;//Neighbour cell
		srcCell = &(nodes[p.pointPairs[i].p1.y][p.pointPairs[i].p1.x]);
		dstCell = &(nodes[p.pointPairs[i].p2.y][p.pointPairs[i].p2.x]);

		//Get source and destination data
		srcCell->cost = double(cstData[srcCell->y][srcCell->x]);
		srcDstHab = p.useHabGrid ? double(habData[srcCell->y][srcCell->x]) * double(habData[dstCell->y][dstCell->x]) : 1.0;

		//Insert source node into heap
		minHeap.insert(srcCell);

		//Perform Least cost path algorithm breaking when cCell = dstCell or cCell is nullptr
		while (true)
		{
			//Get lowest accum cost cell from tree
			cCell = minHeap.popMin();
			//Return nullptr if minHeap has run out of cells. 
			//Should only happen if src and dst are completely seperated by noData
			if (cCell == nullptr) break;

			//return cCell if dst cell is found
			if (cCell->x == dstCell->x && cCell->y == dstCell->y) break;

			//Check we're in some search radius from dst
			//I don't think this is needed, artificial restraint to reduce searching too far in the wrong direction
			//This could be an implementation of the AStar optimisation
			//if (((currentSCell->x - dstPCell->x) * (currentSCell->x - dstPCell->x)) + ((currentSCell->y - dstPCell->y) * (currentSCell->y - dstPCell->y)) > max_path_cells_sqrd) continue;

			//For each neighbour cell
			for (n = 0; n < 8; n++) {
				//Get neighbour's spatial location
				nx = cCell->x + neighbourDirs[n].x;
				ny = cCell->y + neighbourDirs[n].y;

				//Check were in grid bounds and not nodata
				if ((nx | ny) < 0 ||
					nx >= p.nCols ||
					ny >= p.nRows)
					continue;
				
				if (cstData[ny][nx] == p.noData)
					continue;

				//Get neighbour cell and skip if already searched or outside radius
				nCell = &(nodes[ny][nx]);
				if (nCell->popped) continue;

				//Calculate neighbour's accumulated cost and skip if greater than max ED
				nAccumCost = cCell->cost + (double(cstData[ny][nx]) + double(cstData[cCell->y][cCell->x])) * neighbourDists[n];
				//if (p.useMaxED && nAccumCost > p.maxED) continue;

				//Insert nCell into minHeap if it isn't already
				if (nCell->heapPos == 0) {
					nCell->previous = cCell;
					nCell->cost = nAccumCost;
					minHeap.insert(nCell);
					continue;
				}
				//Otherwise update nCell position in heap
				if (nAccumCost < nCell->cost) {
					nCell->previous = cCell;
					nCell->cost = nAccumCost;
					minHeap.decrease(nCell);
					continue;
				}
			}//End for(n)
		}//End while(true)


		 //Save Path to grid
		if (cCell != nullptr)
		{
			pathsCompleted++;
			//apply decay function and habitat multiplier
			d_e = exp(-p.d_i * (cCell->cost - p.d_a));
			pathValue = float(srcDstHab * (d_e) / (1 + d_e));

			//walk path from dst to src through cell->previous pointers
			pCell = cCell;
			while (pCell != nullptr)
			{
				lnkData[pCell->y][pCell->x] += pathValue;
				pCell = pCell->previous;
			}

			//is it time to save?
			//if (i % saveEvery == 0)
			//{
			//	p.lnkGS.open(p.lnkFN);
			//	for (int r = 0; r < p.nRows; r++) {
			//		p.lnkGS.write((const char *)(lnkData[r].data()), sizeof(float) * p.nCols);
			//	}
			//	p.lnkGS.close();
			//}
		}
		else {
			pathsFailed++;
		}


	}//End for loop

	 //Save output and cleanup
	msgString("\rPercent complete: 100\n");
	msgString("Paths completed: " + toStr(pathsCompleted) + "\n");
	msgString("Paths failed: " + toStr(pathsFailed) + "\n");

	p.lnkGS.open(p.lnkGFN);
	p.lnkGS.copyHeader(p.cstGS);
	for (int r = 0; r < p.nRows; r++) {
		p.lnkGS.write((const char *)(lnkData[r].data()), sizeof(float) * p.nCols);
	}
	p.lnkGS.close();

	msgString("Complete!\n");
	return 0;
}


int LoadPointPairs(LinksParams &p) {
	msgString("Loading Point Pairs from " + p.pointPairsTabFN);
	int2 p1, p2;
	std::string line;
	std::ifstream pointPairsIFS(p.pointPairsTabFN, std::ios::in);
	if (!pointPairsIFS.is_open()) return msgErrorStr(-4, p.pointPairsTabFN);

	p.nPointPairs = 0;
	while (std::getline(pointPairsIFS, line)) {
		std::stringstream lineSS(line);
		lineSS >> p1.x >> p1.y >> p2.x >> p2.y;
		p.pointPairs.push_back({ p1, p2 });
		p.nPointPairs++;
	}
	return 0;
}


int CreatePointPairs(LinksParams &p) {
	msgString("Creating " + toStr(p.nPointPairs) + " Point Pairs\n");

	//if (p.writePointGrid)
	//{
	//	p.pntGS.open(p.pntGFN);
	//	p.pntGS.copyHeader(p.cstGS);
	//}
	std::random_device rd;
	std::mt19937 gen(rd());

	p.habGS.open(p.habGFN);
	if (!p.habGS.is_open()) return msgErrorStr(-3, p.habGFN);

	std::ofstream pointPairsOFS;
	if (p.havePointPairsTabFN) {
		pointPairsOFS.open(p.pointPairsTabFN, std::ios::out);
		if (!pointPairsOFS.is_open()) return msgErrorStr(-4, p.pointPairsTabFN);
	}

	//Sample points > habMin to make vector of points
	//TODO Create point grid?
	msgString("Sampling points\n");
	std::vector<float> habRow(p.nCols, 0);
	//std::vector<float> pntRow(p.nCols, 0);
	std::vector<int2> points;

	std::uniform_real_distribution<> uniformReals(0.0, 1.0);

	for (int y = 0; y < p.nRows; y++) {
		msgProgress("Percent complete: ", y * 100 / p.nRows);
		p.habGS.read((char *)habRow.data(), sizeof(float) * p.nCols);
		//std::fill(pntRow.begin(), pntRow.end(), 0);

		for (int x = 0; x < p.nCols; x++) {
			
			if (habRow[x] != float(p.habGS.noData()) &&
				habRow[x] >= p.habMin && 
				habRow[x] * p.P > uniformReals(gen))
			{
				points.push_back({ x, y });
				//pntRow[x] = 1;
			}
		}

		//if (p.writePointGrid)
		//	p.pntGS.write((const char*) pntRow.data(), sizeof(float) * p.nCols);
	}

	msgString("\rPercent complete: 100\n");
	msgString("Number of points: " + toStr(points.size()) + "\n");

	p.habGS.close();

	//Pair up points that are within distance
	msgString("Pairing points\n");
	int p1, p2;
	double pairDistSqrd;

	std::uniform_int_distribution<> uniformInts(0, points.size() - 1);

	for (int n = 0; n < p.nPointPairs; n++) {
		msgProgress("Percent complete: ", n * 100 / p.nPointPairs);
		do {
			p1 = uniformInts(gen);//% points.size();
			p2 = uniformInts(gen);//% points.size();

			pairDistSqrd =
				std::pow(points[p1].x - points[p2].x, 2.0) +
				std::pow(points[p1].y - points[p2].y, 2.0);

		} while (p1 == p2 || pairDistSqrd < p.minPairDistSqrd || pairDistSqrd > p.maxPairDistSqrd);
		p.pointPairs.push_back({ points[p1], points[p2] });
	}

	msgString("\rPercent complete: 100\n");
	msgString("Number of point pairs: " + toStr(p.pointPairs.size()) + "\n");

	//Write tab seperated values without header:
	if (p.havePointPairsTabFN) {
		for (int n = 0; n < p.pointPairs.size(); n++) {
			pointPairsOFS << p.pointPairs[n].p1.x << "\t" << p.pointPairs[n].p1.y << "\t" << p.pointPairs[n].p2.x << "\t" << p.pointPairs[n].p2.y << std::endl;
		}
		pointPairsOFS.close();
	}
	return 0;
}
