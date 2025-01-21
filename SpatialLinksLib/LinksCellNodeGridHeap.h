/*
LinksCellNodeGridHeap.h - defines CellNode, CellNodeGrid and CellNodeHeap classes for implementing Spatial Links custom heap
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
//Implements Dijkstra's least cost path algorithm using a min heap.
//Complete Links needs a node for each pixel in the search window, 
//Random Links needs a node for each pixel in the raster domain
//inRadius maintains binary mask for those in/out of search radius
//nodes maintain their heap index (heapPos) so they can be addressed in O(1)

//TODO
//Move definitions to cpp file?

#pragma once
#ifndef LINKS_CELL_NODE_GRID_HEAP_H
#define LINKS_CELL_NODE_GRID_HEAP_H

#include <vector>

class CellNode {
	//No privates
private:
	//Expose all members for easy access
public:
	//x, y = window location. k = heap position. >1 in heap =0 not in heap;
	CellNode * previous;//8
	double cost;//8
	size_t heapPos;//8
	int x, y;//8
	bool popped;//1
	bool inRadius;//1

	CellNode() : CellNode(0, 0) {}
	CellNode(int winX, int winY) {
		previous = nullptr;
		x = winX;
		y = winY;
		heapPos = 0;
		cost = 0.0;
		popped = false;
		inRadius = false;
	}

	void reset() {
		heapPos = 0;
		cost = 0.0;
		popped = false;
		previous = nullptr;
	}

	~CellNode() { previous = nullptr; }
};

//CellNodeGrid owns the CellNodes, CellNodeHeap just holds pointers
class CellNodeGrid {

public:
	size_t nRows, nCols, centre, r, c;
	std::vector<std::vector<CellNode>> cellNodes;

	//Constructor for square grid(with centre) resizes cellNodes and assigns node locations
	CellNodeGrid() : CellNodeGrid(0) {}
	CellNodeGrid(size_t windowSize) {
		nRows = windowSize;
		nCols = windowSize;
		centre = (windowSize + 1) / 2;
		cellNodes.resize(nRows, std::vector<CellNode>(nCols));
		for (r = 0; r < nRows; ++r) {
			for (c = 0; c < nCols; ++c) {
				cellNodes[r][c].x = int(c);
				cellNodes[r][c].y = int(r);
			}
		}
	}
	CellNodeGrid(size_t windowSize, double searchRadiusCellsSqrd) : CellNodeGrid(windowSize) {
		setSearchRadiusSqrd(searchRadiusCellsSqrd);
	}

	//Constructor for non-square grid (no centre) resizes cellNodes and assigns node locations
	CellNodeGrid(size_t rows, size_t cols) {
		nRows = rows;
		nCols = cols;
		centre = 0;//Don't use centre or inRadius
		cellNodes.resize(nRows, std::vector<CellNode>(nCols));
		for (r = 0; r < nRows; ++r) {
			for (c = 0; c < nCols; ++c) {
				cellNodes[r][c].x = int(c);
				cellNodes[r][c].y = int(r);
			}
		}
	}

	//Reset all cellNodes heap/search info
	void reset() {
		for (r = 0; r < nRows; r++) {
			for (c = 0; c < nCols; c++)	{
				cellNodes[r][c].reset();
			}
		}
	}

	//Flag all cellNodes that are within the search radius
	void setSearchRadiusSqrd(double searchRadiusCellsSqrd) {
		if (centre > 0) {
			for (r = 0; r < nRows; r++) {
				for (c = 0; c < nCols; c++) {
					cellNodes[r][c].inRadius = ((r - centre) * (r - centre)) + ((c - centre) * (c - centre)) < searchRadiusCellsSqrd;
				}
			}
		}
	}
	//Provide [row][col] vector access
	std::vector<CellNode> & operator[](size_t row) { return cellNodes[row]; }
	const std::vector<CellNode> & operator[](size_t row) const { return cellNodes[row]; }

	~CellNodeGrid() {}

};

//CellNodeHeap implements a binary heap using a vector of CellNode pointers 
//This class makes assumptions that favor speed over safety
//That exceed heapMaxSize or in == nullptr will never happen on insert()
//That heapSize < 1 will only happen after construction, reset or last node has been popped.
class CellNodeHeap {
private:
	size_t heapSize;
	size_t parent, min, left, right;
	CellNode *tmp;
	CellNode * minNode;
	std::vector<CellNode*> minHeap;

	//Swap node at k with parent node at k / 2 while its cost is lower
	void upHeap(size_t k) {
		parent = k / 2;
		while (k > 1 && minHeap[parent]->cost > minHeap[k]->cost) {
			
			//parent cost is higher so swap
			std::swap(minHeap[k], minHeap[parent]);

			minHeap[k]->heapPos = k;
			minHeap[parent]->heapPos = parent;

			k = parent;
			parent = k / 2;
		}
	}

	//Swap node at k with child node at 2k or 2k+1 while its cost is higher
	void downHeap(size_t k)
	{
		while (true) {
			min = k;
			left = 2 * k;
			right = 2 * k + 1;
			if (left <= heapSize && minHeap[left]->cost < minHeap[min]->cost) min = left;
			if (right <= heapSize && minHeap[right]->cost < minHeap[min]->cost) min = right;
			if (min == k) break;

			//child with lower cost found so swap
			std::swap(minHeap[k], minHeap[min]);

			minHeap[min]->heapPos = min;
			minHeap[k]->heapPos = k;

			k = min;
		}
	}

public:

	CellNodeHeap() : CellNodeHeap(0) {}
	//Constructor resizes minHeap to size + 1 since heap range is 1..n
	CellNodeHeap(size_t size) {
		heapSize = 0;
		parent = 0;
		min = 0;
		left = 0;
		right = 0;
		tmp = nullptr;
		minNode = nullptr;
		minHeap.resize(size + 1, nullptr);
	}
	
	void reset() {
		//This works since heapSize is only incremented on insert limiting access to only 'new' pointers
		heapSize = 0;
	}

	//Increment heapSize by inserting node at end of heap then upHeap it as needed
	void insert(CellNode *in) {
		minHeap[++heapSize] = in;
		minHeap[heapSize]->heapPos = heapSize;
		upHeap(heapSize);
	}
		
	//Move node up the heap if needed
	//void decrease(int k) {
	void decrease(CellNode *node) {
			upHeap(node->heapPos);
	}

	//Deincrement heapSize by returning root node replacing with last node in heap then downheap it as needed 
	CellNode* popMin() {
		minNode = nullptr;
		//heapSize < 1 should only happen after construction, reset or last node has been popped.
		if (heapSize > 0) {
			minNode = minHeap[1];
			minNode->heapPos = 0;
			minNode->popped = true;

			minHeap[1] = minHeap[heapSize];
			minHeap[1]->heapPos = 1;
			minHeap[heapSize] = nullptr;

			if (--heapSize > 1) downHeap(1);
		}
		//std::cout << "heap pop: " << minNode->cost << std::endl;
		return minNode;
	}
};

#endif
