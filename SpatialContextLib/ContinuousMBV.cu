/*
ContinuousMBV.cpp - Experimental CUDA functions for calculating marginal biodiversity values from continuous GDM data
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

//Sij = sum(exp(sum(-|srcTxi - dstTxj|)))
//This kernel function gets called on the device operating for every pixel in a row of data in parallel
//For 25 tx grids and 500,000 samples this performs 888,427 billion subtractions!

//txSrcData[0 + cell] and txDstSamples[s + 0] store habitat condition for the focal pixel and samples

__global__ void ContinuousMBV_kernel(size_t nSampled, float *txSrcData, float *txDstSamples, float2 *outData)
{
	size_t cell = size_t(blockIdx.x * blockDim.x + threadIdx.x);
	size_t nColsAligned = size_t(d_nColsAligned);
	float sumDist = 0.0f, sumSij = 0.0f, sumHjSij = 0.0f, Hj;

	//Check this thread is in the row
	if (cell < d_nCols) {
		/*
		Get pixel txSrc data, unrolled to use registers rather than cached local memory. First txSrc value is checked for no data.
		txSrcData are stored as end to end nColsAligned wide (pitched) grid data rows so that all 32 threads in a warp read consecutive data from a single coalesced 128byte memory transaction
		First nColsAligned width in bytes is calculated then each line is compiled into two instructions, a 64 bit add incrementing the address by nColsAligned width then a load

		Each line is compiled into two PTX instructions:
		add.s64 	txAddress, txPrevAddress, nColsAlignedBytes;
		ld.global.f32 	txSrcVal, [txAddress];
		*/
		float Hi = txSrcData[cell];
		if (Hi != d_noData) {
			float txSrcPixel1 = txSrcData[1ULL * nColsAligned + cell];
			float txSrcPixel2 = txSrcData[2ULL * nColsAligned + cell];
			float txSrcPixel3 = txSrcData[3ULL * nColsAligned + cell];
			float txSrcPixel4 = txSrcData[4ULL * nColsAligned + cell];
			float txSrcPixel5 = txSrcData[5ULL * nColsAligned + cell];
			float txSrcPixel6 = txSrcData[6ULL * nColsAligned + cell];
			float txSrcPixel7 = txSrcData[7ULL * nColsAligned + cell];
			float txSrcPixel8 = txSrcData[8ULL * nColsAligned + cell];
			float txSrcPixel9 = txSrcData[9ULL * nColsAligned + cell];
			float txSrcPixel10 = txSrcData[10ULL * nColsAligned + cell];
			float txSrcPixel11 = txSrcData[11ULL * nColsAligned + cell];
			float txSrcPixel12 = txSrcData[12ULL * nColsAligned + cell];
			float txSrcPixel13 = txSrcData[13ULL * nColsAligned + cell];
			float txSrcPixel14 = txSrcData[14ULL * nColsAligned + cell];
			float txSrcPixel15 = txSrcData[15ULL * nColsAligned + cell];
			float txSrcPixel16 = txSrcData[16ULL * nColsAligned + cell];
			float txSrcPixel17 = txSrcData[17ULL * nColsAligned + cell];
			float txSrcPixel18 = txSrcData[18ULL * nColsAligned + cell];
			float txSrcPixel19 = txSrcData[19ULL * nColsAligned + cell];
			float txSrcPixel20 = txSrcData[20ULL * nColsAligned + cell];
			float txSrcPixel21 = txSrcData[21ULL * nColsAligned + cell];
			float txSrcPixel22 = txSrcData[22ULL * nColsAligned + cell];
			float txSrcPixel23 = txSrcData[23ULL * nColsAligned + cell];
			float txSrcPixel24 = txSrcData[24ULL * nColsAligned + cell];
			float txSrcPixel25 = txSrcData[25ULL * nColsAligned + cell];
			//Currently only processing 25 tx grids after habitat so cut this short
			//float txSrcPixel26 = txSrcData[26ULL * nColsAligned + cell];
			//float txSrcPixel27 = txSrcData[27ULL * nColsAligned + cell];
			//float txSrcPixel28 = txSrcData[28ULL * nColsAligned + cell];
			//float txSrcPixel29 = txSrcData[29ULL * nColsAligned + cell];
			//float txSrcPixel30 = txSrcData[30ULL * nColsAligned + cell];
			//float txSrcPixel31 = txSrcData[31ULL * nColsAligned + cell];

			/*
			Get txSrc - txDst distances, subtracting from 0, also manually unrolled for tx
			All 32 threads in a warp read txDst from the same address for each tx which should be then broadcast only once to all threads
			As the addresses for each txDst are consecutive they're cached in a single coalesced 128byte memory transaction and should be hit in L1 (or L2) for 1..31

			Each line is compiled into four PTX instructions:
			ld.global.f32 	txDstVal, [txDstSamples + txByteOffset];
			sub.f32 	txDif, txSrcVal, txDstVal;
			abs.f32 	txAbsDif, txDif;
			sub.f32 	sumDist, sumDist, txAbsDif;
			*/
			for (size_t s = 0ULL; s < nSampled; s++) {
				sumDist = 0.0f;
				Hj = txDstSamples[s * 32ULL + 0ULL];
				sumDist -= fabsf(txSrcPixel1 - txDstSamples[s * 32ULL + 1ULL]);
				sumDist -= fabsf(txSrcPixel2 - txDstSamples[s * 32ULL + 2ULL]);
				sumDist -= fabsf(txSrcPixel3 - txDstSamples[s * 32ULL + 3ULL]);
				sumDist -= fabsf(txSrcPixel4 - txDstSamples[s * 32ULL + 4ULL]);
				sumDist -= fabsf(txSrcPixel5 - txDstSamples[s * 32ULL + 5ULL]);
				sumDist -= fabsf(txSrcPixel6 - txDstSamples[s * 32ULL + 6ULL]);
				sumDist -= fabsf(txSrcPixel7 - txDstSamples[s * 32ULL + 7ULL]);
				sumDist -= fabsf(txSrcPixel8 - txDstSamples[s * 32ULL + 8ULL]);
				sumDist -= fabsf(txSrcPixel9 - txDstSamples[s * 32ULL + 9ULL]);
				sumDist -= fabsf(txSrcPixel10 - txDstSamples[s * 32ULL + 10ULL]);
				sumDist -= fabsf(txSrcPixel11 - txDstSamples[s * 32ULL + 11ULL]);
				sumDist -= fabsf(txSrcPixel12 - txDstSamples[s * 32ULL + 12ULL]);
				sumDist -= fabsf(txSrcPixel13 - txDstSamples[s * 32ULL + 13ULL]);
				sumDist -= fabsf(txSrcPixel14 - txDstSamples[s * 32ULL + 14ULL]);
				sumDist -= fabsf(txSrcPixel15 - txDstSamples[s * 32ULL + 15ULL]);
				sumDist -= fabsf(txSrcPixel16 - txDstSamples[s * 32ULL + 16ULL]);
				sumDist -= fabsf(txSrcPixel17 - txDstSamples[s * 32ULL + 17ULL]);
				sumDist -= fabsf(txSrcPixel18 - txDstSamples[s * 32ULL + 18ULL]);
				sumDist -= fabsf(txSrcPixel19 - txDstSamples[s * 32ULL + 19ULL]);
				sumDist -= fabsf(txSrcPixel20 - txDstSamples[s * 32ULL + 20ULL]);
				sumDist -= fabsf(txSrcPixel21 - txDstSamples[s * 32ULL + 21ULL]);
				sumDist -= fabsf(txSrcPixel22 - txDstSamples[s * 32ULL + 22ULL]);
				sumDist -= fabsf(txSrcPixel23 - txDstSamples[s * 32ULL + 23ULL]);
				sumDist -= fabsf(txSrcPixel24 - txDstSamples[s * 32ULL + 24ULL]);
				sumDist -= fabsf(txSrcPixel25 - txDstSamples[s * 32ULL + 25ULL]);
				//Currently only processing 25 tx grids after habitat so cut this short
				//sumDist -= fabsf(txSrcPixel26 - txDstSamples[s * 32ULL + 26ULL]);
				//sumDist -= fabsf(txSrcPixel27 - txDstSamples[s * 32ULL + 27ULL]);
				//sumDist -= fabsf(txSrcPixel28 - txDstSamples[s * 32ULL + 28ULL]);
				//sumDist -= fabsf(txSrcPixel29 - txDstSamples[s * 32ULL + 29ULL]);
				//sumDist -= fabsf(txSrcPixel30 - txDstSamples[s * 32ULL + 30ULL]);
				//sumDist -= fabsf(txSrcPixel31 - txDstSamples[s * 32ULL + 31ULL]);
				//Calculate Sij for this pixel
				sumSij += __expf(sumDist);
				sumHjSij += Hj * __expf(sumDist);
			}
			//Write Sij to the pixel in the output row
			outData[cell].x = sumSij;
			outData[cell].y = sumHjSij;
		}
		else {
			outData[cell].x = d_noData;
			outData[cell].y = d_noData;
		}
	}
}

int CUDAContinuousMBV(CBAParams &p)
{
	CUDA(cudaDeviceSetCacheConfig(cudaFuncCache::cudaFuncCachePreferL1));

	float *d_txDstSamples;
	float *d_txSrcData;
	float2 *d_outData;
	uint nThreads = 512;
	uint nBlocks = p.nCols / nThreads + 1;
	dim3 blockDim{ nThreads, 1U, 1U };
	dim3 gridDim{ nBlocks, 1U, 1U };
	uint nColsAligned = int(ceil(float(p.nCols) / 32.0f) * 32);

	std::unique_ptr<float[]> habData = std::make_unique<float[]>(p.nCells);
	std::unique_ptr<float[]> h_txDstSamples = std::make_unique<float[]>(p.maxSamples * p.maxTxGrids);
	std::unique_ptr<float[]> h_txSrcData = std::make_unique<float[]>(nColsAligned * p.maxTxGrids);
	std::unique_ptr<float2[]> h_outData = std::make_unique<float2[]>(nColsAligned);
	std::unique_ptr<float[]> sampleGrid = std::make_unique<float[]>(p.nCells);
	std::unique_ptr<float[]> txData = std::make_unique<float[]>(p.nCells);

	for (uint i = 0; i < p.maxSamples * p.maxTxGrids; i++) h_txDstSamples[i] = 0.0f;
	for (uint i = 0; i < nColsAligned * p.maxTxGrids; i++) h_txSrcData[i] = 0.0f;
	for (uint i = 0; i < p.nCells; i++) sampleGrid[i] = 0.0f;

	//Sample all txDst grids
	std::cout << "Sampling " << p.nTxGrids << " grids\n";

	float txVal = p.noData;
	float habVal = p.noData;
	float habNoData = p.noData;
	float txNoData = p.noData;
	uint nSampled = 0, prevSampled = 0;

	p.habInFS.read((char*)(habData.get()), p.nCells * sizeof(float));

	for (uint tx = 0; tx < p.nTxGrids; tx++) {
		msgProgress("Percent complete: ", tx * 100 / p.nTxGrids);
		p.txDstGSs[tx].read((char*)(txData.get()), p.nCells * sizeof(float));
		txNoData = float(p.txDstGSs[tx].noData());
		nSampled = 0;
		for (uint y = p.sampleStep; y < p.nRows - p.sampleStep; y += p.sampleStep) {
			for (uint x = p.sampleStep; x < p.nCols - p.sampleStep; x += p.sampleStep) {
				txVal = txData[y * p.nCols + x];
				
				if (txVal == txNoData) continue;
				h_txDstSamples[nSampled * p.maxTxGrids + tx + 1] = txVal;// + 1 for habitat
				sampleGrid[y * p.nCols + x] = 1.0f;
				nSampled++;

				if (tx == 0) {
					habVal = habData[y * p.nCols + x];
					if (habVal == habNoData) continue;
					h_txDstSamples[nSampled * p.maxTxGrids + tx] = habVal;
				}

				if (nSampled == p.maxSamples) break;
			}
			if (nSampled == p.maxSamples) break;
		}
		if (tx > 0 && nSampled != prevSampled) return -1;
		prevSampled = nSampled;
	}
	msgText("\rPercent complete: 100");

	habData.reset();
	txData.reset();

	if (nSampled == 0) return -1;
	std::cout << "\npixels sampled: \t\t" << nSampled << "\n";

	//Write samples to output table and grid
	for (uint s = 0; s < nSampled; s++)
		for (uint tx = 0; tx < p.nTxGrids; tx++)
			p.txDstSamplesFS << h_txDstSamples[s * p.maxTxGrids + tx] << (tx + 1 < p.nTxGrids ? ", " : "\n");

	p.sampleGridFS.write((const char *)(sampleGrid.get()), p.nCols * p.nRows * sizeof(float));
	p.sampleGridFS.close();
	sampleGrid.reset();

	//Allocate device memory and copy data from host
	CUDA(cudaMalloc(&d_txDstSamples, nSampled * p.maxTxGrids * sizeof(float)));
	CUDA(cudaMalloc(&d_txSrcData, nColsAligned * p.maxTxGrids * sizeof(float)));
	CUDA(cudaMalloc(&d_outData, nColsAligned * sizeof(float2)));
	CUDA(cudaMemcpyToSymbol(d_nCols, &(p.nCols), sizeof(uint)));
	CUDA(cudaMemcpyToSymbol(d_nColsAligned, &(nColsAligned), sizeof(uint)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(habNoData), sizeof(uint)));
	CUDA(cudaMemcpy(d_txDstSamples, h_txDstSamples.get(), nSampled * p.maxTxGrids * sizeof(float), cudaMemcpyHostToDevice));

	std::cout << "Processing " << p.nCells << " grid cells on the GPU\n";
	for (uint r = 0; r < p.nRows; r++) {
		msgProgress("Percent complete: ", r * 100 / p.nRows);

		//Set the output row to 0
		CUDA(cudaMemset(d_outData, 0, nColsAligned * sizeof(float2)));

		//Read in  txSrc rows and copy to device
		for (uint tx = 0; tx < p.nTxGrids; tx++)
			p.txSrcGSs[tx].read((char*)&(h_txSrcData[(tx + 1) * nColsAligned]), p.nCols * sizeof(float));// + 1 for habitat
		CUDA(cudaMemcpy(d_txSrcData, h_txSrcData.get(), nColsAligned * p.maxTxGrids * sizeof(float), cudaMemcpyHostToDevice));

		//Call kernel to calculate Sij for all pixels in the row
		ContinuousMBV_kernel <<< gridDim, blockDim >>> (size_t(nSampled), d_txSrcData, d_txDstSamples, d_outData);

		//Copy row from device to host then write to disk
		CUDA(cudaMemcpy(h_outData.get(), d_outData, nColsAligned * sizeof(float2), cudaMemcpyDeviceToHost));
		for (uint c = 0; c < p.nCols; c++) {
			p.denOutFS.write((const char*)&(h_outData[c].x), sizeof(float));
			p.numOutFS.write((const char*)&(h_outData[c].y), sizeof(float));
		}
	}
	msgText("\rPercent complete: 100");

	p.cxtOutFS.close();
	//Clean up device memory
	CUDA(cudaFree(d_txDstSamples));
	CUDA(cudaFree(d_txSrcData));
	CUDA(cudaFree(d_outData));

	cudaError_t cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDAResilienceDenom() Complete!");
	return int(cudaStatus);
}
#endif