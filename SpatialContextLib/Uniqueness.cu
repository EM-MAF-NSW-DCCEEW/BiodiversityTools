/*
Uniqueness.cu - CUDA functions for calculating uniqueness from continuous GDM data
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

#ifndef _UNIQUENESS_KERNEL_CU_
#define _UNIQUENESS_KERNEL_CU_
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "Common.cuh"
#include "Parameters.h"
#include "FileUtilsLib.h"

//Sij = exp(b + sum(-|srcTxi - dstTxj|))
//Dij = sum(1 - exp(b + sum(-|srcTxi - dstTxj|)))
//This kernel function gets called on the device operating for every pixel in a row of data in parallel

__global__ void Uniqueness_kernel(size_t nSampled, float* txSrcData, float* txDstSamples, float* outData)
{
	size_t cell = size_t(blockIdx.x * blockDim.x + threadIdx.x);
	size_t nColsAligned = size_t(d_nColsAligned);
	float sumDist = 0.0f, sumDij = 0.0f;
	size_t txAddress = cell;

	//Check this thread is in the row
	if (cell < d_nCols) {
		/*
		Get pixel txSrc data, unrolled to use registers rather than cached local memory. First txSrc value is checked for no data.
		txSrcData are stored as end to end nColsAligned wide (pitched) grid data rows so that all 32 threads in a warp should read consecutive data from a single coalesced 128byte memory transaction
		First nColsAligned width in bytes is calculated then each line is compiled into two instructions, a 64 bit add incrementing the address by nColsAligned width then a load

		(Current PTX not reviewed) Each line was compiled into two PTX instructions:
		add.s64 	txAddress, txPrevAddress, nColsAlignedBytes;
		ld.global.f32 	txSrcVal, [txAddress];
		*/
		float txSrcPixel0 = txSrcData[txAddress];
		if (txSrcPixel0 != d_noData) {
			float txSrcPixel1 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel2 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel3 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel4 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel5 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel6 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel7 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel8 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel9 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel10 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel11 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel12 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel13 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel14 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel15 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel16 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel17 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel18 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel19 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel20 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel21 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel22 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel23 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel24 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel25 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel26 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel27 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel28 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel29 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel30 = txSrcData[txAddress += nColsAligned];
			float txSrcPixel31 = txSrcData[txAddress += nColsAligned];
			
			/*
			Get txSrc - txDst distances, subtracting from 0, also manually unrolled for tx
			All 32 threads in a warp read txDst from the same address for each tx which should be  broadcast only once to all threads
			As the addresses for each txDst are consecutive they should be cached in a single coalesced 128byte memory transaction and hit in L1 (or L2) for 1..31

			(Current PTX not reviewed) Each line was compiled into four PTX instructions:
			ld.global.f32 	txDstVal, [txDstSamples + txByteOffset];
			sub.f32 	txDif, txSrcVal, txDstVal;
			abs.f32 	txAbsDif, txDif;
			sub.f32 	sumDist, sumDist, txAbsDif;
			*/
			for (size_t s = 0ULL; s < nSampled; s++) {
				txAddress = s * 32ULL;
				sumDist = 0.0f;
				sumDist -= d_xIntercept;
				sumDist -= fabsf(txSrcPixel0 - txDstSamples[txAddress + 0ULL]);
				sumDist -= fabsf(txSrcPixel1 - txDstSamples[txAddress + 1ULL]);
				sumDist -= fabsf(txSrcPixel2 - txDstSamples[txAddress + 2ULL]);
				sumDist -= fabsf(txSrcPixel3 - txDstSamples[txAddress + 3ULL]);
				sumDist -= fabsf(txSrcPixel4 - txDstSamples[txAddress + 4ULL]);
				sumDist -= fabsf(txSrcPixel5 - txDstSamples[txAddress + 5ULL]);
				sumDist -= fabsf(txSrcPixel6 - txDstSamples[txAddress + 6ULL]);
				sumDist -= fabsf(txSrcPixel7 - txDstSamples[txAddress + 7ULL]);
				sumDist -= fabsf(txSrcPixel8 - txDstSamples[txAddress + 8ULL]);
				sumDist -= fabsf(txSrcPixel9 - txDstSamples[txAddress + 9ULL]);
				sumDist -= fabsf(txSrcPixel10 - txDstSamples[txAddress + 10ULL]);
				sumDist -= fabsf(txSrcPixel11 - txDstSamples[txAddress + 11ULL]);
				sumDist -= fabsf(txSrcPixel12 - txDstSamples[txAddress + 12ULL]);
				sumDist -= fabsf(txSrcPixel13 - txDstSamples[txAddress + 13ULL]);
				sumDist -= fabsf(txSrcPixel14 - txDstSamples[txAddress + 14ULL]);
				sumDist -= fabsf(txSrcPixel15 - txDstSamples[txAddress + 15ULL]);
				sumDist -= fabsf(txSrcPixel16 - txDstSamples[txAddress + 16ULL]);
				sumDist -= fabsf(txSrcPixel17 - txDstSamples[txAddress + 17ULL]);
				sumDist -= fabsf(txSrcPixel18 - txDstSamples[txAddress + 18ULL]);
				sumDist -= fabsf(txSrcPixel19 - txDstSamples[txAddress + 19ULL]);
				sumDist -= fabsf(txSrcPixel20 - txDstSamples[txAddress + 20ULL]);
				sumDist -= fabsf(txSrcPixel21 - txDstSamples[txAddress + 21ULL]);
				sumDist -= fabsf(txSrcPixel22 - txDstSamples[txAddress + 22ULL]);
				sumDist -= fabsf(txSrcPixel23 - txDstSamples[txAddress + 23ULL]);
				sumDist -= fabsf(txSrcPixel24 - txDstSamples[txAddress + 24ULL]);
				sumDist -= fabsf(txSrcPixel25 - txDstSamples[txAddress + 25ULL]);
				sumDist -= fabsf(txSrcPixel26 - txDstSamples[txAddress + 26ULL]);
				sumDist -= fabsf(txSrcPixel27 - txDstSamples[txAddress + 27ULL]);
				sumDist -= fabsf(txSrcPixel28 - txDstSamples[txAddress + 28ULL]);
				sumDist -= fabsf(txSrcPixel29 - txDstSamples[txAddress + 29ULL]);
				sumDist -= fabsf(txSrcPixel30 - txDstSamples[txAddress + 30ULL]);
				sumDist -= fabsf(txSrcPixel31 - txDstSamples[txAddress + 31ULL]);

				//Calculate Dij for this pixel
				sumDij += 1.0 - __expf(sumDist);
			}
			
			//Write Dij to the pixel in the output row
			outData[cell] = sumDij;
		}
		else {
			outData[cell] = d_noData;
		}
	}
}

int CUDAUniqueness(CBAParams& p)
{
	CUDA(cudaDeviceSetCacheConfig(cudaFuncCache::cudaFuncCachePreferL1));

	float* d_txDstSamples;
	float* d_txSrcData;
	float* d_outData;
	uint nThreads = 256;
	uint nBlocks = p.nCols / nThreads + 1;
	dim3 blockDim{ nThreads, 1U, 1U };
	dim3 gridDim{ nBlocks, 1U, 1U };
	uint nColsAligned = int(ceil(float(p.nCols) / 32.0f) * 32);

	//Allocate host memory
	std::unique_ptr<float[]> h_txDstSamples = std::make_unique<float[]>(p.maxSamples * p.maxTxGrids);
	std::unique_ptr<float[]> h_txSrcData = std::make_unique<float[]>(nColsAligned * p.maxTxGrids);
	std::unique_ptr<float[]> h_outData = std::make_unique<float[]>(nColsAligned);
	std::unique_ptr<float[]> sampleGrid = std::make_unique<float[]>(p.nCells);
	std::unique_ptr<float[]> txData = std::make_unique<float[]>(p.nCells);

	//Initialize host memory
	for (uint i = 0; i < p.maxSamples * p.maxTxGrids; i++) h_txDstSamples[i] = 0.0f;
	for (uint i = 0; i < nColsAligned * p.maxTxGrids; i++) h_txSrcData[i] = 0.0f;
	for (uint i = 0; i < p.nCells; i++) sampleGrid[i] = 0.0f;

	//Sample all txDst grids
	std::cout << "Sampling " << p.nTxGrids << " grids\n";

	float txVal = float(p.txDstGSs[0].noData());
	float txNoData = float(p.txDstGSs[0].noData());
	uint nSampled = 0, prevSampled = 0;

	for (uint tx = 0; tx < p.nTxGrids; tx++) {
		msgProgress("Percent complete: ", tx * 100 / p.nTxGrids);
		
		p.txDstGSs[tx].read((char*)(txData.get()), p.nCells * sizeof(float));

		nSampled = 0;
		
		for (uint y = p.sampleStep; y < p.nRows - p.sampleStep; y += p.sampleStep) {
			for (uint x = p.sampleStep; x < p.nCols - p.sampleStep; x += p.sampleStep) {
				
				txVal = txData.get()[y * p.nCols + x];
				if (txVal == txNoData) continue;
				h_txDstSamples[nSampled * p.maxTxGrids + tx] = txVal;
				sampleGrid[y * p.nCols + x] = float(nSampled);
				nSampled++;
				if (nSampled == p.maxSamples) break;
			}

			if (nSampled == p.maxSamples) break;
		}

		if (tx > 0 && nSampled != prevSampled) return msgErrorStr(-99, "Number of samples from " + p.txDstGSs[tx].fileName() + " don't match");
		prevSampled = nSampled;
	}

	msgText("\rPercent complete: 100");

	txData.reset();

	if (nSampled == 0) return -1;
	std::cout << "\npixels sampled: \t\t" << nSampled << "\n";

	//Write sample data to file
	p.txDstSamplesFS << "Samples,";
	for (uint tx = 0; tx < p.nTxGrids; tx++)
		p.txDstSamplesFS << p.txDstGSs[tx].fileName() << (tx + 1 < p.nTxGrids ? "," : "\n");
	
	for (uint s = 0; s < nSampled; s++) {
		p.txDstSamplesFS << s << ",";
		for (uint tx = 0; tx < p.nTxGrids; tx++)
			p.txDstSamplesFS << h_txDstSamples[s * p.maxTxGrids + tx] << (tx + 1 < p.nTxGrids ? "," : "\n");
	}
	p.txDstSamplesFS.close();

	//Write sample grid to file
	p.sampleGridGS.write((const char*)(sampleGrid.get()), p.nCells * sizeof(float));
	p.sampleGridGS.close();
	sampleGrid.reset();

	//Allocate device memory and copy data from host
	CUDA(cudaMalloc(&d_txDstSamples, nSampled * p.maxTxGrids * sizeof(float)));
	CUDA(cudaMalloc(&d_txSrcData, nColsAligned * p.maxTxGrids * sizeof(float)));
	CUDA(cudaMalloc(&d_outData, nColsAligned * sizeof(float)));
	
	CUDA(cudaMemcpyToSymbol(d_nCols, &(p.nCols), sizeof(uint)));
	CUDA(cudaMemcpyToSymbol(d_nColsAligned, &(nColsAligned), sizeof(uint)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(txNoData), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_xIntercept, &(p.xIntercept), sizeof(float)));

	CUDA(cudaMemcpy(d_txDstSamples, h_txDstSamples.get(), nSampled * p.maxTxGrids * sizeof(float), cudaMemcpyHostToDevice));

	//Process each row on the GPU
	std::cout << "Processing " << p.nCells << " grid cells on the GPU\n";

	//Overlap reading txSrc rows and processing
	//Read in txSrc rows and copy to device
	for (uint tx = 0; tx < p.nTxGrids; tx++)
		p.txSrcGSs[tx].read((char*)&(h_txSrcData[tx * nColsAligned]), p.nCols * sizeof(float));

	for (uint r = 0; r < p.nRows; r++) {
		msgProgress("Percent complete: ", r * 100 / p.nRows);

		//Copy txSrc rows to device
		CUDA(cudaMemcpy(d_txSrcData, h_txSrcData.get(), nColsAligned * p.maxTxGrids * sizeof(float), cudaMemcpyHostToDevice));

		//Call kernel to calculate Sij for all pixels in the row
		Uniqueness_kernel << < gridDim, blockDim >> > (size_t(nSampled), d_txSrcData, d_txDstSamples, d_outData);

		//Read in txSrc rows and copy to device
		if (r + 1 < p.nRows) {
			for (uint tx = 0; tx < p.nTxGrids; tx++)
				p.txSrcGSs[tx].read((char*)&(h_txSrcData[tx * nColsAligned]), p.nCols * sizeof(float));
		}
		//Copy output row from device to host then write to disk
		CUDA(cudaMemcpy(h_outData.get(), d_outData, nColsAligned * sizeof(float), cudaMemcpyDeviceToHost));
		p.outUniquenessGS.write((const char*)h_outData.get(), p.nCols * sizeof(float));

	}

	msgText("\rPercent complete: 100");

	p.outUniquenessGS.close();
	//Clean up device memory
	CUDA(cudaFree(d_txDstSamples));
	CUDA(cudaFree(d_txSrcData));
	CUDA(cudaFree(d_outData));

	cudaError_t cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDAUniqueness() Complete!");
	return int(cudaStatus);
}
#endif 
/* _UNIQUENESS_KERNEL_CU_ */
