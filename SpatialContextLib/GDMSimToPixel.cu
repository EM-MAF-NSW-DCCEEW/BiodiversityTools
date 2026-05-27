/*
GDMSimToPixel.cu - Experimental CUDA functions for calculating GDM similarity to a single pixel
Copyright(C) 2023 State of New South Wales and Department of Planning and Environment
Author: Jamie Love, SEI Metrics and Forecasting

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


#ifndef _RESILIENCE_DENOM_KERNEL_CU_
#define _RESILIENCE_DENOM_KERNEL_CU_
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "Common.cuh"
#include "Parameters.h"
#include "FileUtilsLib.h"


//Compares sample of txDst grids to each pixel in a row of txSrc data and calculates the GDM similarity (Sij) for each pixel in the row

//TODO process multiple rows
//Sij = sum(exp(sum(-|srcTxi - dstTxj|)))
//This kernel function gets called on the device operating for every pixel in a row of data in parallel
__global__ void GDMSimToPixel_kernel(float *txSrcData, float *txDstSamples, float *outData)
{
	size_t cell = size_t(blockIdx.x * blockDim.x + threadIdx.x);
	size_t nColsAligned = size_t(d_nColsAligned);
	float sumDist = 0.0f, sumSij = 0.0f;
	
	//TODO address with txAddress as per Uniqueness
	//size_t txAddress = cell;

	//Check this thread is in the row
	if (cell < d_nCols) {

		float txSrcPixel0 = txSrcData[cell];
		if (txSrcPixel0 != -9999.0f) {

			sumDist = 0.0f;
			sumDist -= fabsf(txSrcPixel0 - txDstSamples[0ULL]);
			sumDist -= fabsf(txSrcData[1ULL * nColsAligned + cell] - txDstSamples[1ULL]);
			sumDist -= fabsf(txSrcData[2ULL * nColsAligned + cell] - txDstSamples[2ULL]);
			sumDist -= fabsf(txSrcData[3ULL * nColsAligned + cell] - txDstSamples[3ULL]);
			sumDist -= fabsf(txSrcData[4ULL * nColsAligned + cell] - txDstSamples[4ULL]);
			sumDist -= fabsf(txSrcData[5ULL * nColsAligned + cell] - txDstSamples[5ULL]);
			sumDist -= fabsf(txSrcData[6ULL * nColsAligned + cell] - txDstSamples[6ULL]);
			sumDist -= fabsf(txSrcData[7ULL * nColsAligned + cell] - txDstSamples[7ULL]);
			sumDist -= fabsf(txSrcData[8ULL * nColsAligned + cell] - txDstSamples[8ULL]);
			sumDist -= fabsf(txSrcData[9ULL * nColsAligned + cell] - txDstSamples[9ULL]);
			sumDist -= fabsf(txSrcData[10ULL * nColsAligned + cell] - txDstSamples[10ULL]);
			sumDist -= fabsf(txSrcData[11ULL * nColsAligned + cell] - txDstSamples[11ULL]);
			sumDist -= fabsf(txSrcData[12ULL * nColsAligned + cell] - txDstSamples[12ULL]);
			sumDist -= fabsf(txSrcData[13ULL * nColsAligned + cell] - txDstSamples[13ULL]);
			sumDist -= fabsf(txSrcData[14ULL * nColsAligned + cell] - txDstSamples[14ULL]);
			sumDist -= fabsf(txSrcData[15ULL * nColsAligned + cell] - txDstSamples[15ULL]);
			sumDist -= fabsf(txSrcData[16ULL * nColsAligned + cell] - txDstSamples[16ULL]);
			sumDist -= fabsf(txSrcData[17ULL * nColsAligned + cell] - txDstSamples[17ULL]);
			sumDist -= fabsf(txSrcData[18ULL * nColsAligned + cell] - txDstSamples[18ULL]);
			sumDist -= fabsf(txSrcData[19ULL * nColsAligned + cell] - txDstSamples[19ULL]);
			sumDist -= fabsf(txSrcData[20ULL * nColsAligned + cell] - txDstSamples[20ULL]);
			sumDist -= fabsf(txSrcData[21ULL * nColsAligned + cell] - txDstSamples[21ULL]);
			sumDist -= fabsf(txSrcData[22ULL * nColsAligned + cell] - txDstSamples[22ULL]);
			sumDist -= fabsf(txSrcData[23ULL * nColsAligned + cell] - txDstSamples[23ULL]);
			sumDist -= fabsf(txSrcData[24ULL * nColsAligned + cell] - txDstSamples[24ULL]);
			sumDist -= fabsf(txSrcData[25ULL * nColsAligned + cell] - txDstSamples[25ULL]);
			sumDist -= fabsf(txSrcData[26ULL * nColsAligned + cell] - txDstSamples[26ULL]);
			sumDist -= fabsf(txSrcData[27ULL * nColsAligned + cell] - txDstSamples[27ULL]);
			sumDist -= fabsf(txSrcData[28ULL * nColsAligned + cell] - txDstSamples[28ULL]);
			sumDist -= fabsf(txSrcData[29ULL * nColsAligned + cell] - txDstSamples[29ULL]);
			sumDist -= fabsf(txSrcData[30ULL * nColsAligned + cell] - txDstSamples[30ULL]);
			sumDist -= fabsf(txSrcData[31ULL * nColsAligned + cell] - txDstSamples[31ULL]);
			//Calculate Sij for this pixel
			sumSij += __expf(sumDist);

			//Write Sij to the pixel in the output row
			outData[cell] = sumSij;
		}
		else outData[cell] = -9999.0f;
	}
}

int CUDAGDMSimToPixel(CBAParams &p)
{
	CUDA(cudaDeviceSetCacheConfig(cudaFuncCache::cudaFuncCachePreferL1));

	float *d_txDstSample;
	float *d_txSrcData;
	float *d_outData;

	uint nThreads = 256;
	uint nBlocks = p.nCols / nThreads + 1;
	dim3 blockDim{ nThreads, 1U, 1U };
	dim3 gridDim{ nBlocks, 1U, 1U };
	uint nColsAligned = int(ceil(float(p.nCols) / 32.0f) * 32);

	std::unique_ptr<float[]> h_txDstSample = std::make_unique<float[]>(p.maxTxGrids);
	std::unique_ptr<float[]> h_txSrcData = std::make_unique<float[]>(nColsAligned * p.maxTxGrids);
	std::unique_ptr<float[]> h_outData = std::make_unique<float[]>(nColsAligned);

	for (uint i = 0; i < p.maxTxGrids; i++) h_txDstSample[i] = 0.0f;
	for (uint i = 0; i < nColsAligned * p.maxTxGrids; i++) h_txSrcData[i] = 0.0f;

	std::unique_ptr<float[]> pixelGrid = std::make_unique<float[]>(p.nCells);
	p.pixelFS.read((char*)pixelGrid.get(), p.nCells * sizeof(float));

	//Sample all txDst grids
	std::unique_ptr<float[]> txData = std::make_unique<float[]>(p.nCells);
	std::cout << "Sampling " << p.nTxGrids << " grids\n";

	uint pixel = 0;
	for (; pixel < p.nCells; pixel++) if (pixelGrid[pixel] > FLT_EPSILON) break;
	pixelGrid.reset();

	uint txSampled = 0;

	for (uint tx = 0; tx < p.nTxGrids; tx++) {
		msgProgress("Percent complete: ", tx * 100 / p.nTxGrids);
		p.txDstFSs[tx].read((char*)txData.get(), p.nCells * sizeof(float));
		if (txData[pixel] == p.noData) continue;
		h_txDstSample[tx] = txData[pixel];
		txSampled++;
	}
	msgText("\rPercent complete: 100");

	txData.reset();
	if (txSampled < p.nTxGrids) return -1;

	//Allocate device memory and copy data from host
	CUDA(cudaMalloc(&d_txDstSample, p.maxTxGrids * sizeof(float)));
	CUDA(cudaMalloc(&d_txSrcData, nColsAligned * p.maxTxGrids * sizeof(float)));
	CUDA(cudaMalloc(&d_outData, nColsAligned * sizeof(float)));

	CUDA(cudaMemcpy(d_txDstSample, h_txDstSample.get(), p.maxTxGrids * sizeof(float), cudaMemcpyHostToDevice));
	CUDA(cudaMemcpyToSymbol(d_nCols, &(p.nCols), sizeof(uint)));
	CUDA(cudaMemcpyToSymbol(d_nColsAligned, &(nColsAligned), sizeof(uint)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(p.noData), sizeof(uint)));
	

	std::cout << "Processing " << p.nCells << " grid cells on the GPU\n";
	for (uint r = 0; r < p.nRows; r++) {
		msgProgress("Percent complete: ", r * 100 / p.nRows);

		//Set the output row to 0
		CUDA(cudaMemset(d_outData, 0, nColsAligned * sizeof(float)));

		//Read in  txSrc rows and copy to device
		for (uint tx = 0; tx < p.nTxGrids; tx++)
			p.txSrcFSs[tx].read((char*)&(h_txSrcData.get()[tx * nColsAligned]), p.nCols * sizeof(float));
		CUDA(cudaMemcpy(d_txSrcData, h_txSrcData.get(), nColsAligned * p.maxTxGrids * sizeof(float), cudaMemcpyHostToDevice));

		//Call kernel to calculate Sij for all pixels in the row
		GDMSimToPixel_kernel << < gridDim, blockDim >> > (d_txSrcData, d_txDstSample, d_outData);

		//Copy row from device to host then write to disk
		CUDA(cudaMemcpy(h_outData.get(), d_outData, nColsAligned * sizeof(float), cudaMemcpyDeviceToHost));
		p.cxtOutFS.write((const char*)h_outData.get(), p.nCols * sizeof(float));
	}
	msgText("\rPercent complete: 100");
	p.cxtOutFS.close();

	//Clean up device memory
	CUDA(cudaFree(d_txDstSample));
	CUDA(cudaFree(d_txSrcData));
	CUDA(cudaFree(d_outData));

	cudaError_t cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDAResilienceDenom() Complete!");
	return int(cudaStatus);
}

int CUDAGDMSimToPixel(CBAParams& p, uint x, uint y)
{
	CUDA(cudaDeviceSetCacheConfig(cudaFuncCache::cudaFuncCachePreferL1));

	float* d_txDstSample;
	float* d_txSrcData;
	float* d_outData;

	uint nThreads = 256;
	uint nBlocks = p.nCols / nThreads + 1;
	dim3 blockDim{ nThreads, 1U, 1U };
	dim3 gridDim{ nBlocks, 1U, 1U };
	uint nColsAligned = int(ceil(float(p.nCols) / 32.0f) * 32);

	//Check x and y are within bounds and return error if not
	if (x >= p.nCols || y >= p.nRows) return msgErrorStr(-99, "Pixel coordinates are out of bounds");

	std::vector<float> h_txDstSample(p.maxTxGrids, 0.0f);
	std::vector<float> h_txSrcData(nColsAligned * p.maxTxGrids, 0.0f);
	std::vector<float> h_outData(nColsAligned, 0.0f);

	uint pixel = y * p.nCols + x;

	for (uint tx = 0; tx < p.nTxGrids; tx++) {
		p.txDstFSs[tx].seekg(pixel * sizeof(float), std::ios::beg);
		p.txDstFSs[tx].read((char*)&(h_txDstSample[tx]), sizeof(float));
		if (h_txDstSample[tx] == p.noData) return msgErrorStr(-99, "noData value found at sample pixel in TX" + toStr(tx) + "DST");
	}

	//Allocate device memory and copy data from host
	CUDA(cudaMalloc(&d_txDstSample, p.maxTxGrids * sizeof(float)));
	CUDA(cudaMalloc(&d_txSrcData, nColsAligned * p.maxTxGrids * sizeof(float)));
	CUDA(cudaMalloc(&d_outData, nColsAligned * sizeof(float)));

	CUDA(cudaMemcpy(d_txDstSample, h_txDstSample.data(), p.maxTxGrids * sizeof(float), cudaMemcpyHostToDevice));
	CUDA(cudaMemcpyToSymbol(d_nCols, &(p.nCols), sizeof(uint)));
	CUDA(cudaMemcpyToSymbol(d_nColsAligned, &(nColsAligned), sizeof(uint)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(p.noData), sizeof(uint)));

	//Read in  first txSrc rows
	for (uint tx = 0; tx < p.nTxGrids; tx++)
		p.txSrcFSs[tx].read((char*)&(h_txSrcData[tx * nColsAligned]), p.nCols * sizeof(float));

	for (uint r = 0; r < p.nRows; r++) {
		msgProgress("Percent complete: ", r * 100 / p.nRows);

		//Set the output row to 0
		CUDA(cudaMemset(d_outData, 0, nColsAligned * sizeof(float)));

		//Copy txSrc rows to device
		CUDA(cudaMemcpy(d_txSrcData, h_txSrcData.data(), nColsAligned * p.maxTxGrids * sizeof(float), cudaMemcpyHostToDevice));

		//Call kernel to calculate Sij for all pixels in the row
		GDMSimToPixel_kernel << < gridDim, blockDim >> > (d_txSrcData, d_txDstSample, d_outData);

		//Read in next txSrc rows
		for (uint tx = 0; tx < p.nTxGrids; tx++)
			p.txSrcFSs[tx].read((char*)&(h_txSrcData[tx * nColsAligned]), p.nCols * sizeof(float));

		//Copy row from device to host then write to disk
		CUDA(cudaMemcpy(h_outData.data(), d_outData, nColsAligned * sizeof(float), cudaMemcpyDeviceToHost));
		p.cxtOutFS.write((const char*)h_outData.data(), p.nCols * sizeof(float));
	}
	msgText("\rPercent complete: 100");
	p.cxtOutFS.close();

	//Clean up device memory
	CUDA(cudaFree(d_txDstSample));
	CUDA(cudaFree(d_txSrcData));
	CUDA(cudaFree(d_outData));

	cudaError_t cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDAResilienceDenom() Complete!");
	return int(cudaStatus);
}
#endif