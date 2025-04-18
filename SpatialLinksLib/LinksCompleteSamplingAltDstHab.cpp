/*
LinksCompleteSamplingAltDstHab.cpp - complete sampling Spatial Links implementaion with alt dst habitat using new custom heap
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

//Notes
//TODO 
//Standardise parameter file format

#include "LinksCompleteSamplingAltDstHab.h"

int LinksCompleteSamplingAltDstHab(std::string paramFN, LinksParams &p) {

	std::string srcHabGFN, dstHabGFN, cstGFN, pntGFN, lnkGFN, cxtGFN, outParamFN;

	if (!GetParam(paramFN, "INPUTS", "CostGrid", cstGFN)) msgErrorStr(-1, "CostGrid", paramFN);
	if (!GetParam(paramFN, "INPUTS", "SrcHabGrid", srcHabGFN)) msgErrorStr(-1, "SrcHabGrid", paramFN);
	if (!GetParam(paramFN, "INPUTS", "DstHabGrid", dstHabGFN)) msgErrorStr(-1, "DstHabGrid", paramFN);
	p.useHabGrid = true;

	p.doLinks = GetParam(paramFN, "OUTPUTS", "LinksGrid", lnkGFN);
	p.doContext = GetParam(paramFN, "OUTPUTS", "CxtGrid", cxtGFN);
	if (!p.doLinks && !p.doContext)	return msgErrorStr(-99, "No LinksGrid or CxtGrid parameters Nothing to do");

	if (!GetParam(paramFN, "SEARCH", "SearchRadius", p.searchRadius)) msgErrorStr(-1, "SearchRadius", paramFN);
	if (!GetParam(paramFN, "SEARCH", "MinHabThreshold", p.habMin)) msgErrorStr(-1, "MinHabThreshold", paramFN);
	p.useMaxED = GetParam(paramFN, "SEARCH", "MaxED", p.maxED);

	if (!GetParam(paramFN, "DECAY", "I", p.d_i)) msgErrorStr(-1, "I", paramFN);
	if (!GetParam(paramFN, "DECAY", "A", p.d_a)) msgErrorStr(-1, "A", paramFN);

	if (!GetParam(paramFN, "OPTIONS", "UseHabWeightLinks", p.useSrcDstHabWeight)) p.useSrcDstHabWeight = p.useHabGrid;
	if (!GetParam(paramFN, "OPTIONS", "UseHabWeightCxt", p.useFCellHabWeight)) p.useFCellHabWeight = p.useHabGrid;
	if (!GetParam(paramFN, "OPTIONS", "LinksExponent", p.lnkExponent)) p.lnkExponent = 0.25;
	if (!GetParam(paramFN, "OPTIONS", "nthreads", p.nThreads)) p.nThreads = 1;
	p.useMTProcessing = p.nThreads > 1;

	//Open cost grid
	p.cstGS.open(cstGFN);
	if (!p.cstGS.is_open()) return msgErrorStr(-3, cstGFN);
	p.cstGS.getHeader(p.nCols, p.nRows, p.xll, p.yll, p.cellSize, p.noData);
	p.nCells = p.nCols * p.nRows;

	//Open src and dst hab grids
	p.srcHabGS.open(srcHabGFN);
	if (!p.srcHabGS.is_open()) return msgErrorStr(-3, srcHabGFN);
	if (!p.srcHabGS.compareHeader(p.cstGS)) return msgErrorStr(-7, srcHabGFN, cstGFN);

	p.dstHabGS.open(dstHabGFN);
	if (!p.dstHabGS.is_open()) return msgErrorStr(-3, dstHabGFN);
	if (!p.dstHabGS.compareHeader(p.cstGS)) return msgErrorStr(-7, dstHabGFN, cstGFN);

	//Open links grid for writing
	if (p.doLinks) {
		p.lnkGS.open(lnkGFN);
		p.lnkGS.copyHeader(p.cstGS);
		//TODO Test is_open()
		//if (!lnkGS.is_open()) return msgErrorStr(-4, lnkGFN);
	}

	//Open context grid for writing
	if (p.doContext) {
		p.cxtGS.open(cxtGFN);
		p.cxtGS.copyHeader(p.cstGS);
		//TODO Test is_open()
		//if (!cxtGS.is_open()) return msgErrorStr(-4, cxtGFN);
	}

	int result = 0;
	if (p.useMTProcessing) {
		//TODO multithreaded implementation;
		msgErrorStr(-99, "TODO multithreaded implementation");
	}
	else
	{
		result = LinksCompleteSamplingAltDstHab_ST(p);
	}
	return result;
}

int LinksCompleteSamplingAltDstHab_ST(LinksParams &p) {
	msgText("Performing single-threaded Spatial Links analysis with complete grid sampling using alt dst hab");

	//Calculate search parameters
	double searchRadiusCellsSqrd = (p.searchRadius / p.cellSize) * (p.searchRadius / p.cellSize);
	int windowSize = (int)(p.searchRadius / p.cellSize * 2);
	windowSize += windowSize % 2 == 0 ? 3 : 4;
	int windowCenter = (windowSize + 1) / 2;
	int lastRow = windowSize - 1;

	int windowOffset = 0;
	double nAccumCost = 0.0;
	double d_e;
	float srcHab, dstHab;
	double srcDstHab = 1.0;
	float pathValue;
	int n, nx, ny;

	//Initialise CellNodes
	msgText("Allocating memory");
	CellNodeGrid nodes(windowSize, searchRadiusCellsSqrd);
	CellNodeHeap minHeap(windowSize * windowSize);
	CellNode *nCell, *cCell, *pCell;

	//Init data
	std::vector<std::vector<float>> cstData(windowSize, std::vector<float>(p.nCols, p.noData));
	std::vector<std::vector<float>> srcHabData(windowSize, std::vector<float>(p.nCols, p.noData));
	std::vector<std::vector<float>> dstHabData(windowSize, std::vector<float>(p.nCols, p.noData));
	std::vector<std::vector<float>> lnkData(windowSize, std::vector<float>(p.nCols, 0.0f));
	std::vector<std::vector<float>> cxtData(windowSize, std::vector<float>(p.nCols, 0.0f));

	std::vector<float> zeroRow(p.nCols, 0.0f);
	std::vector<float> oneRow(p.nCols, 1.0f);
	std::vector<float> nullRow(p.nCols, p.noData);

	//Start reading input cost and habitat data
	msgText("Reading input data to memory");
	for (int r = windowCenter; r < windowSize; r++) {
		p.cstGS.read((char*)(cstData[r].data()), sizeof(float) * p.nCols);
		p.srcHabGS.read((char*)(srcHabData[r].data()), sizeof(float) * p.nCols);
		p.dstHabGS.read((char*)(dstHabData[r].data()), sizeof(float) * p.nCols);
	}

	//Process each row in input data
	msgText("Processing least cost paths");
	for (int focalRow = 0; focalRow < p.nRows; focalRow++) {
		msgProgress("Processing row: ", focalRow);

		//Process each focal cell in row
		for (int focalCol = 0; focalCol < p.nCols; focalCol++) {

			//Focal cell is no data
			//Focal cell src or dst hab?
			srcHab = srcHabData[windowCenter][focalCol];
			if (cstData[windowCenter][focalCol] == p.noData || srcHab == p.noData) {
				if (p.doLinks) lnkData[windowCenter][focalCol] = p.noData;
				if (p.doContext) cxtData[windowCenter][focalCol] = p.noData;
				continue;
			}

			//Don't process focal cells below habitat threshold
			if (srcHab < p.habMin) continue;

			//Reset nodes then insert focal cell to root of tree
			nAccumCost = 0.0;
			nodes.reset();
			nCell = nullptr;
			cCell = &(nodes[windowCenter][windowCenter]);
			cCell->cost = double(cstData[windowCenter][focalCol]);
			minHeap.insert(cCell);
			windowOffset = focalCol - windowCenter;

			//Process Window while cells to search
			while (true) {
				//Pop lowest cost cell or break when empty (when every cell searched)
				cCell = minHeap.popMin();
				if (cCell == nullptr) break;
				dstHab = dstHabData[cCell->y][cCell->x + windowOffset];
				
				//Output Context by applying decay kernel then adding to focal cell
				if (p.doContext)
					cxtData[windowCenter][focalCol] += float(exp(-(1.0 / p.d_a) * cCell->cost) * double(dstHab));

				//Output Links by applying decay kernel then adding to each (previous) cell in LCP 
				if (p.doLinks && dstHab >= p.habMin) {
					srcDstHab = double(srcHab) * double(dstHab);
					d_e = exp(-p.d_i * (double(cCell->cost) - p.d_a));
					pathValue = float(d_e / (1.0 + d_e) * srcDstHab);
					pCell = cCell;

					while (pCell != nullptr) {
						lnkData[pCell->y][pCell->x + windowOffset] += pathValue;
						pCell = pCell->previous;
					}
				}

				//For each neighbour cell
				for (n = 0; n < 8; n++) {
					//Get neighbour's spatial location
					nx = cCell->x + neighbourDirs[n].x;
					ny = cCell->y + neighbourDirs[n].y;

					//Check we're in grid bounds
					if (nx + windowOffset >= (signed)p.nCols ||
						nx + windowOffset < 0 ||
						ny >= windowSize ||
						ny < 0)
						continue;

					//Check for no data
					if (cstData[ny][nx + windowOffset] == p.noData) continue;

					//Get neighbour cell and skip if already solved or outside radius
					nCell = &(nodes[ny][nx]);
					if (nCell->popped || !nCell->inRadius) continue;

					//Calculate neighbour's accumulated cost and skip if greater than max ED
					nAccumCost = cCell->cost + (double(cstData[ny][nx + windowOffset]) + double(cstData[cCell->y][cCell->x + windowOffset])) * neighbourDists[n];
					if (p.useMaxED && nAccumCost > p.maxED) continue;

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

				}//End For each Neighbour
			}//End while(true)

			 //Context Focal Cell Mod
			if (p.doContext) cxtData[windowCenter][focalCol] *= p.useFCellHabWeight ? srcHabData[windowCenter][focalCol] : 1.0F;

		}//End for each focalCol

		 //Shuffle up input data
		for (int r = 0; r < lastRow; r++) {
			std::memcpy(cstData[r].data(), cstData[r + 1].data(), p.nCols * sizeof(float));
			std::memcpy(srcHabData[r].data(), srcHabData[r + 1].data(), p.nCols * sizeof(float));
			std::memcpy(dstHabData[r].data(), dstHabData[r + 1].data(), p.nCols * sizeof(float));
		}

		//Read next row from files into input data or pad input data beyond end of files
		if (focalRow + windowCenter <= p.nRows) {
			p.cstGS.read((char*)(cstData[lastRow].data()), sizeof(float) * p.nCols);
			p.srcHabGS.read((char*)(srcHabData[lastRow].data()), sizeof(float) * p.nCols);
			p.dstHabGS.read((char*)(dstHabData[lastRow].data()), sizeof(float) * p.nCols);
		}
		else {
			std::memcpy(cstData[lastRow].data(), nullRow.data(), p.nCols * sizeof(float));
			std::memcpy(srcHabData[lastRow].data(), nullRow.data(), p.nCols * sizeof(float));
			std::memcpy(dstHabData[lastRow].data(), nullRow.data(), p.nCols * sizeof(float));
		}

		//Write top row to files then shuffle up
		if (p.doLinks) {
			if (focalRow >= windowCenter) {

				//Apply lnkExponent before writing to file
				for (uint c = 0; c < p.nCols; c++) {
					if (lnkData[0][c] != p.noData)
						lnkData[0][c] = (float)std::pow(lnkData[0][c], p.lnkExponent);
				}
				p.lnkGS.write((const char*)(lnkData[0].data()), sizeof(float) * p.nCols);
			}

			for (int r = 0; r < lastRow; r++)
				std::memcpy(lnkData[r].data(), lnkData[r + 1].data(), p.nCols * sizeof(float));
			std::memcpy(lnkData[lastRow].data(), zeroRow.data(), p.nCols * sizeof(float));
		}

		if (p.doContext) {
			if (focalRow >= windowCenter)
				p.cxtGS.write((const char*)(cxtData[0].data()), sizeof(float) * p.nCols);

			for (int r = 0; r < lastRow; r++)
				std::memcpy(cxtData[r].data(), cxtData[r + 1].data(), p.nCols * sizeof(float));
			std::memcpy(cxtData[lastRow].data(), zeroRow.data(), p.nCols * sizeof(float));
		}

	}//End for each focalRow

	msgText("\rFinished processing rows.");
	//Write last rows to files
	if (p.doLinks) {
		for (int r = 0; r < windowCenter - 1; r++) {

			//Apply lnkExponent before writing to file
			for (uint c = 0; c < p.nCols; c++) {
				if (lnkData[r][c] != p.noData)
					lnkData[r][c] = (float)std::pow(lnkData[r][c], p.lnkExponent);
			}
			p.lnkGS.write((const char*)(lnkData[r].data()), sizeof(float) * p.nCols);
		}
	}

	if (p.doContext) {
		for (int r = 0; r < windowCenter - 1; r++)
			p.cxtGS.write((const char*)(cxtData[r].data()), sizeof(float) * p.nCols);
	}

	p.lnkGS.close();
	p.cxtGS.close();

	msgText("Finished writing outputs.");
	msgText("CompleteLinksAltDstHabST() Complete!");
	return 0;
}