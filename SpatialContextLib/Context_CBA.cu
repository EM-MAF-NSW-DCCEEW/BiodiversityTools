/*
Context_CBA.cu - CUDA CBA functions for performing BFT context analysis
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

#ifndef _CONTEXT_CBA_KERNEL_CU_
#define _CONTEXT_CBA_KERNEL_CU_
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "Common.cuh"
#include "Parameters.h"
#include "FileUtilsLib.h"

#ifndef USE_PRECALC_EDGE_KERNELS
__global__ void ContextCBA_kernel(int inOffset, int outOffset, int2 *d_petalData, float *d_outData, cudaTextureObject_t tex2Obj)
{
	//Focal cell
	//float2 focal = tex2D(tex2Ref, blockIdx.x, blockIdx.y + inOffset);
	float2 focal = tex2D<float2>(tex2Obj, blockIdx.x, blockIdx.y + inOffset);
	if ((focal.x == 0.0f && d_multFocal) || focal.y == 0.0f)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = 0.0f;
	else if (focal.x == d_noData || focal.y == d_noData)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = d_noData;
	else {
		//Declare thread local storage
		int i;
		float Hj = 0.0f;
		int2 intPairs[4];
		float2 t{ -1.0f, -1.0f };

		//Address shared memory
		extern __shared__ float pCost[];		//Permeability of this petal
		float *pDist = &pCost[blockDim.x];		//Highest permeability back to the focal cell
		float *pCount = &pCost[2 * blockDim.x];	//Count of the number of petals used
		pCost[threadIdx.x] = 0.0f;
		pDist[threadIdx.x] = 0.0f;
		pCount[threadIdx.x] = 0.0f;

		//Aggregate cell values into petals
		intPairs[0] = d_petalData[4 * blockDim.x + threadIdx.x];      //petal cell count, valid cell count = 0
		intPairs[0].y = 0;
		for (i = 0; i < intPairs[0].x; i++) {
			intPairs[1] = d_petalData[(i + 5) * blockDim.x + threadIdx.x];//plus 5 to jump neighbours storage
			if (intPairs[1].x != -1 && intPairs[1].y != -1) {
				//t = tex2D(tex2Ref, float(int(blockIdx.x) + intPairs[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + intPairs[1].y - d_focalOffset));
				t = tex2D<float2>(tex2Obj, float(int(blockIdx.x) + intPairs[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + intPairs[1].y - d_focalOffset));
				if (t.x >= 0.0f && t.y > 0.0f) {
					Hj += t.x;
					pCost[threadIdx.x] += t.y;
					++(intPairs[0].y);
				}
			}
		}

		//Normalise petal cost and set petal counter
		if (intPairs[0].y > 0) {
			pCost[threadIdx.x] = powf(pCost[threadIdx.x] / float(intPairs[0].y), 2.0f * sqrtf(float(intPairs[0].y) * d_oneOnPi));
			pCount[threadIdx.x] = 1.0f;
		}

		//Get petal neighbour indices N, S, W, E, NW, NE, SW, SE
		intPairs[0] = d_petalData[0 * blockDim.x + threadIdx.x];
		intPairs[1] = d_petalData[1 * blockDim.x + threadIdx.x];
		intPairs[2] = d_petalData[2 * blockDim.x + threadIdx.x];
		intPairs[3] = d_petalData[3 * blockDim.x + threadIdx.x];

		//Calculate permeability for eight focal cell neighbours first
		if (threadIdx.x < 8)
			pDist[threadIdx.x] = threadIdx.x % 2 == 0 ? powf(pCost[threadIdx.x] * focal.y, d_diag) : sqrtf(pCost[threadIdx.x] * focal.y);

		//Loop until all permeabilities are maximised (i==0)
		i = 1;
		while (__syncthreads_or(i)) {
			i = 0;
			if (intPairs[0].x != -1 && pDist[intPairs[0].x] * sqrtf(pCost[threadIdx.x] * pCost[intPairs[0].x]) > pDist[threadIdx.x]) {//N                
				pDist[threadIdx.x] = pDist[intPairs[0].x] * sqrtf(pCost[threadIdx.x] * pCost[intPairs[0].x]);
				i++;
			}
			if (intPairs[0].y != -1 && pDist[intPairs[0].y] * sqrtf(pCost[threadIdx.x] * pCost[intPairs[0].y]) > pDist[threadIdx.x]) {//S                
				pDist[threadIdx.x] = pDist[intPairs[0].y] * sqrtf(pCost[threadIdx.x] * pCost[intPairs[0].y]);
				i++;
			}
			if (intPairs[1].x != -1 && pDist[intPairs[1].x] * sqrtf(pCost[threadIdx.x] * pCost[intPairs[1].x]) > pDist[threadIdx.x]) {//E                
				pDist[threadIdx.x] = pDist[intPairs[1].x] * sqrtf(pCost[threadIdx.x] * pCost[intPairs[1].x]);
				i++;
			}
			if (intPairs[1].y != -1 && pDist[intPairs[1].y] * sqrtf(pCost[threadIdx.x] * pCost[intPairs[1].y]) > pDist[threadIdx.x]) {//W                
				pDist[threadIdx.x] = pDist[intPairs[1].y] * sqrtf(pCost[threadIdx.x] * pCost[intPairs[1].y]);
				i++;
			}
			if (intPairs[2].x != -1 && pDist[intPairs[2].x] * powf(pCost[threadIdx.x] * pCost[intPairs[2].x], d_diag) > pDist[threadIdx.x]) {//NW                
				pDist[threadIdx.x] = pDist[intPairs[2].x] * powf(pCost[threadIdx.x] * pCost[intPairs[2].x], d_diag);
				i++;
			}
			if (intPairs[2].y != -1 && pDist[intPairs[2].y] * powf(pCost[threadIdx.x] * pCost[intPairs[2].y], d_diag) > pDist[threadIdx.x]) {//NE               
				pDist[threadIdx.x] = pDist[intPairs[2].y] * powf(pCost[threadIdx.x] * pCost[intPairs[2].y], d_diag);
				i++;
			}
			if (intPairs[3].x != -1 && pDist[intPairs[3].x] * powf(pCost[threadIdx.x] * pCost[intPairs[3].x], d_diag) > pDist[threadIdx.x]) {//SW               
				pDist[threadIdx.x] = pDist[intPairs[3].x] * powf(pCost[threadIdx.x] * pCost[intPairs[3].x], d_diag);
				i++;
			}
			if (intPairs[3].y != -1 && pDist[intPairs[3].y] * powf(pCost[threadIdx.x] * pCost[intPairs[3].y], d_diag) > pDist[threadIdx.x]) {//SE               
				pDist[threadIdx.x] = pDist[intPairs[3].y] * powf(pCost[threadIdx.x] * pCost[intPairs[3].y], d_diag);
				i++;
			}
		}

		//Calculate HjWij
		pDist[threadIdx.x] *= Hj;

		//Sum all HjWij and petal count
		__syncthreads();
		for (i = d_firstReduction; i > 0; i /= 2) {
			if (threadIdx.x < i && threadIdx.x + i < d_nPetals) {
				pDist[threadIdx.x] += pDist[threadIdx.x + i];
				pCount[threadIdx.x] += pCount[threadIdx.x + i];
			}
			__syncthreads();
		}

		//Final focal cell calculations then write NHA to outData
		if (threadIdx.x == 0) {
			pDist[0] += focal.x;
			pDist[0] *= (pCount[0] < float(d_nPetals) ? float(d_nPetals) / (pCount[0] + 1.0f) : 1.0f);
			pDist[0] *= d_multFocal ? powf(focal.x, d_focalPower) : 1.0f;
			d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = powf(pDist[0], d_sumPower);
		}
	}
}

#else
//NOTE Do NOT use causes non - reproducable errors
__global__ void ContextCBA_kernel(int inOffset, int outOffset, int2 *d_petalData, float *d_outData)
{
	//Get an check the focal cell
	float2 focal = tex2D(tex2Ref, blockIdx.x, blockIdx.y + inOffset);
	if (focal.y == 0.0f)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = 0.0f;
	else if (focal.x == d_noData || focal.y == d_noData)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = d_noData;
	else {

		//Declare thread local storage
		int i;
		float Hj = 0.0f;
		int2 intPairs[4];
		float2 costNS{ -1.0f, -1.0f }; //also used for fetching cell values
		float2 costWE{ -1.0f, -1.0f };
		float2 costND{ -1.0f, -1.0f };
		float2 costSD{ -1.0f, -1.0f };

		//Address shared memory
		extern __shared__ float pDist[];
		int *pCount = (int*)&pDist[blockDim.x];
		pDist[threadIdx.x] = 0.0f; //Cost, Distance
		pCount[threadIdx.x] = 0; //Active, Count


								 //Calculate petal values from cells
		intPairs[0] = d_petalData[4 * blockDim.x + threadIdx.x]; //petal cell count
		intPairs[0].y = 0; //valid cell count = 0
		for (i = 0; i < intPairs[0].x; i++) {
			intPairs[1] = d_petalData[(i + 5) * blockDim.x + threadIdx.x]; //plus 5 to jump neighbours storage
			if (intPairs[1].x != -1 && intPairs[1].y != -1) {
				costNS = tex2D(tex2Ref, float(int(blockIdx.x) + intPairs[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + intPairs[1].y - d_focalOffset));
				if (costNS.y > d_fltEps) {
					Hj += costNS.x;
					pDist[threadIdx.x] += costNS.y;
					++(intPairs[0].y);
				}
			}
		}

		//Normalise petal cost and set petal counter
		if (intPairs[0].y > 0) {
			pDist[threadIdx.x] = powf(pDist[threadIdx.x] / (float)intPairs[0].y, 2.0f * sqrtf((float)intPairs[0].y * d_oneOnPi));
			pCount[threadIdx.x] = 1;
		}

		//Get petal neighbour indices
		intPairs[0] = d_petalData[0 * blockDim.x + threadIdx.x];//N, S
		intPairs[1] = d_petalData[1 * blockDim.x + threadIdx.x];//W, E
		intPairs[2] = d_petalData[2 * blockDim.x + threadIdx.x];//NW, NE
		intPairs[3] = d_petalData[3 * blockDim.x + threadIdx.x];//SW, SE

															   //Calculate petal edge costs
		__syncthreads();
		costNS.x = intPairs[0].x != -1 ? sqrtf(pDist[threadIdx.x] * pDist[intPairs[0].x]) : -1.0f;
		costNS.y = intPairs[0].y != -1 ? sqrtf(pDist[threadIdx.x] * pDist[intPairs[0].y]) : -1.0f;
		costWE.x = intPairs[1].x != -1 ? sqrtf(pDist[threadIdx.x] * pDist[intPairs[1].x]) : -1.0f;
		costWE.y = intPairs[1].y != -1 ? sqrtf(pDist[threadIdx.x] * pDist[intPairs[1].y]) : -1.0f;
		costND.x = intPairs[2].x != -1 ? powf(pDist[threadIdx.x] * pDist[intPairs[2].x], d_diag) : -1.0f;
		costND.y = intPairs[2].y != -1 ? powf(pDist[threadIdx.x] * pDist[intPairs[2].y], d_diag) : -1.0f;
		costSD.x = intPairs[3].x != -1 ? powf(pDist[threadIdx.x] * pDist[intPairs[3].x], d_diag) : -1.0f;
		costSD.y = intPairs[3].y != -1 ? powf(pDist[threadIdx.x] * pDist[intPairs[3].y], d_diag) : -1.0f;

		//Calculate permeability for focal cell neighbours first
		if (threadIdx.x < 8)
			pDist[threadIdx.x] = threadIdx.x % 2 == 0 ? powf(pDist[threadIdx.x] * focal.y, d_diag) : sqrtf(pDist[threadIdx.x] * focal.y);
		else
			pDist[threadIdx.x] = 0.0f;

		//Loop until all permeabilities are maximised (i==0) then multiply by Hj 
		i = 1;
		while (__syncthreads_or(i)) {
			i = 0;
			if (intPairs[0].x != -1 && pDist[intPairs[0].x] * costNS.x > pDist[threadIdx.x]) {//N edge
				pDist[threadIdx.x] = pDist[intPairs[0].x] * costNS.x;
				i = 1;
			}
			if (intPairs[0].y != -1 && pDist[intPairs[0].y] * costNS.y > pDist[threadIdx.x]) {//S edge
				pDist[threadIdx.x] = pDist[intPairs[0].y] * costNS.y;
				i = 1;
			}
			if (intPairs[1].x != -1 && pDist[intPairs[1].x] * costWE.x > pDist[threadIdx.x]) {//W edge
				pDist[threadIdx.x] = pDist[intPairs[1].x] * costWE.x;
				i = 1;
			}
			if (intPairs[1].y != -1 && pDist[intPairs[1].y] * costWE.y > pDist[threadIdx.x]) {//E edge
				pDist[threadIdx.x] = pDist[intPairs[1].y] * costWE.y;
				i = 1;
			}
			if (intPairs[2].x != -1 && pDist[intPairs[2].x] * costND.x > pDist[threadIdx.x]) {//NW edge
				pDist[threadIdx.x] = pDist[intPairs[2].x] * costND.x;
				i = 1;
			}
			if (intPairs[2].y != -1 && pDist[intPairs[2].y] * costND.y > pDist[threadIdx.x]) {//NE edge
				pDist[threadIdx.x] = pDist[intPairs[2].y] * costND.y;
				i = 1;
			}
			if (intPairs[3].x != -1 && pDist[intPairs[3].x] * costSD.x > pDist[threadIdx.x]) {//SW edge
				pDist[threadIdx.x] = pDist[intPairs[3].x] * costSD.x;
				i = 1;
			}
			if (intPairs[3].y != -1 && pDist[intPairs[3].y] * costSD.y > pDist[threadIdx.x]) {//SE edge
				pDist[threadIdx.x] = pDist[intPairs[3].y] * costSD.y;
				i = 1;
			}
		}
		pDist[threadIdx.x] *= Hj;

		//Sum all HjWij and petal count
		__syncthreads();
		for (i = d_firstReduction; i > 0; i /= 2) {
			if (threadIdx.x < i && threadIdx.x + i < d_nPetals) {
				pDist[threadIdx.x] += pDist[threadIdx.x + i];
				pCount[threadIdx.x] += pCount[threadIdx.x + i];
			}
			__syncthreads();
		}

		//Final focal cell calcs then write result to outData
		if (threadIdx.x == 0) {
			pDist[0] += focal.x;
			if (pCount[0] < d_nPetals) pDist[0] *= float(d_nPetals) / (float(pCount[0]) + 1.0f);
			if (d_multFocal) pDist[0] *= powf(focal.x, d_focalPower);
			d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = powf(pDist[0], d_sumPower);
		}
	}
}
#endif


////////////////////////////////////////////////////////////
//Context CBA small
////////////////////////////////////////////////////////////
int CUDAContextCBA_S(CBAParams &p)
{
	cudaError_t cudaStatus;
	msgText("Performing Context CBA");

	//Host and Device data
	float2 *h_inData, *d_inData;
	float *h_outData, *d_outData;
	int2 *d_petalData;

	//Local variables
	size_t d_pitch;
	dim3 gridSize{ p.nCols, 1U, 1U };
	uint maxKernelTime = h_maxKernelTime; // 100000;
	uint kernelTime = maxKernelTime;
	uint nSharedSlots = 4;
	uint i;
	
	//Host 3 x 4bytes x nCells, Device 3 x 4bytes x nCells
	//Malloc pinned host and device memory 
	msgText("Allocating host and device memory");
	CUDA(cudaMallocHost(&h_inData, p.nCells * sizeof(float2)));
	CUDA(cudaMallocHost(&h_outData, p.nCells * sizeof(float)));
	CUDA(cudaMalloc(&d_petalData, sizeof(int2) *p. petalRows * p.petalCols));
	CUDA(cudaMallocPitch(&d_inData, &d_pitch, p.nCols * sizeof(float2), p.nRows));
	CUDA(cudaMalloc(&d_outData, p.nCells * sizeof(float)));

	//Read h_inData from disk
	msgText("Reading input data to host memory");
	std::unique_ptr<float[]> inData = std::make_unique<float[]>(p.nCells);
	p.habInFS.read((char*)inData.get(), p.nCells * sizeof(float));
	for (i = 0; i < p.nCells; i++) { h_inData[i].x = inData[i]; }
	p.prmInFS.read((char*)inData.get(), p.nCells * sizeof(float));
	for (i = 0; i < p.nCells; i++) { h_inData[i].y = inData[i]; }
	inData.reset();

	//Move input data to device
	msgText("Moving input data to device memory");
	CUDA(cudaMemcpy(d_petalData, p.petalPtr.get(), sizeof(int2) * p.petalRows * p.petalCols, cudaMemcpyHostToDevice));
	CUDA(cudaMemcpy2D(d_inData, d_pitch, h_inData, p.nCols * sizeof(float2), p.nCols * sizeof(float2), p.nRows, cudaMemcpyHostToDevice));
	CUDA(cudaFreeHost(h_inData));

	//Copy constants to device
	msgText("Setting device parameters");
	CUDA(cudaMemcpyToSymbol(d_nPetals, &(p.nPetals), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_focalOffset, &(p.fOffset), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_multFocal, &(p.multFocal), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_firstReduction, &(p.firstReduction), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_sumPower, &(p.sumPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_focalPower, &(p.focalPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(p.noData), sizeof(float)));

	//Testing dispersal probability
	CUDA(cudaMemcpyToSymbol(d_doAveragePetalHabitat, &(p.doAveragePetalHabitat), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_doDispersalProbability, &(p.doDispersalProbability), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_dispersalYears, &(p.dispersalYears), sizeof(float)));
	
	//Texture Reference
	//msgText("Setting texture reference");
	//tex2Ref.filterMode     = cudaFilterModePoint;
	//tex2Ref.addressMode[0] = cudaAddressModeBorder;
	//tex2Ref.addressMode[1] = cudaAddressModeBorder;
	//cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float2>();
	//CUDA(cudaBindTexture2D (NULL, &tex2Ref, d_inData, &channelDesc, p.nCols, p.nRows, d_pitch));

	//Replace Texture Reference with Texture Object
	struct cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypePitch2D;
	resDesc.res.pitch2D.desc = cudaCreateChannelDesc<float2>();
	resDesc.res.pitch2D.devPtr = d_inData;
	resDesc.res.pitch2D.width = p.nCols;
	resDesc.res.pitch2D.height = p.nRows;
	resDesc.res.pitch2D.pitchInBytes = d_pitch;
	
	struct cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.filterMode = cudaFilterModePoint;
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.readMode = cudaReadModeElementType;

	cudaTextureObject_t tex2Obj = 0;
	CUDA(cudaCreateTextureObject(&tex2Obj, &resDesc, &texDesc, NULL));
	

	//Profiled Kernel Call Loop
	msgText("Processing data on device");
	Profiler profiler(1000000), kernelTimer(1000000);
	profiler.Start();
	for (i = 0; i < p.nRows; i += gridSize.y) {
		//Set gridSize.y
		gridSize.y = max(gridSize.y * maxKernelTime / kernelTime, 1);
		gridSize.y = min(gridSize.y, p.nRows - i);
		kernelTimer.Start();
		ContextCBA_kernel <<<gridSize, p.petalCols, p.petalCols * nSharedSlots * sizeof(float) >>>(i, i, d_petalData, d_outData, tex2Obj);
		msgProgress("Percent complete: ", i * 100 / p.nRows);
		CUDA(cudaDeviceSynchronize());
		kernelTime = int(kernelTimer.Stop());
	}
	profiler.Stop();
	msgText("\rPercent complete: 100");
	msgText(("Processing time: " + toStr(profiler.Total())).c_str());

	//Copy device output data to host and write to disk
	msgText("Writing output data to disk");
	CUDA(cudaMemcpy(h_outData, d_outData, p.nCells * sizeof(float), cudaMemcpyDeviceToHost));
	p.cxtOutFS.write((const char*)h_outData, p.nCells * sizeof(float));

	//Free pinned host and device memory
	msgText("Freeing host and device memory");
	CUDA(cudaDestroyTextureObject(tex2Obj));
	CUDA(cudaFreeHost(h_outData));
	CUDA(cudaFree(d_petalData));
	CUDA(cudaFree(d_inData));
	CUDA(cudaFree(d_outData));
	
	//Check CUDA errors
	cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDAContextCBA_S() Complete!");
	return int(cudaStatus);
}

////////////////////////////////////////////////////////////
//Context CBA large
////////////////////////////////////////////////////////////
int CUDAContextCBA_L(CBAParams &p)
{
	cudaError_t cudaStatus;
	msgText("Performing Context CBA");

	//Host and Device data and buffers
	float2 *h_inData, *d_inData;
	float2 *h_inBuf, *d_inBuf;
	float *h_outBuf, *d_outBuf;
	int2 *d_petalData;

	size_t d_pitch;
	uint nBufferRows = p.fOffset * 2;
	uint nDataRows = p.fOffset * 4;
	uint nFirstRows = nDataRows - p.fOffset;
	ull nReads = 0;
	dim3 gridSize( (unsigned int) p.nCols, 1U, 1U);
	uint maxKernelTime = h_maxKernelTime;
	uint kernelTime = maxKernelTime;
	uint nSharedSlots = 4;
	uint i, j;

	//Malloc pinned host memory
	msgText("Allocating host memory");
	CUDA(cudaMallocHost(&h_inData, p.nCols * nDataRows * sizeof(float2)));
	CUDA(cudaMallocHost(&h_inBuf, p.nCols * nBufferRows * sizeof(float2)));
	CUDA(cudaMallocHost(&h_outBuf, p.nCols * nBufferRows * sizeof(float)));
	
	//Malloc device data
	msgText("Allocating device memory");
	CUDA(cudaMalloc(&d_petalData, sizeof(int2) * p.petalRows * p.petalCols));
	CUDA(cudaMallocPitch(&d_inData, &d_pitch, p.nCols * sizeof(float2), nDataRows));
	CUDA(cudaMallocPitch(&d_inBuf, &d_pitch, p.nCols * sizeof(float2), nBufferRows));
	CUDA(cudaMalloc(&d_outBuf, p.nCols * nBufferRows * sizeof(float)));
	
	//Read h_inData from disk
	msgText("Reading input data to host memory");
	for (i = 0; i < nDataRows * p.nCols; i++) {
		h_inData[i].x = 0;
		h_inData[i].y = 0;
		if (i >= p.fOffset * p.nCols) {
			p.habInFS.read((char*)&(h_inData[i].x), sizeof(float));
			p.prmInFS.read((char*)&(h_inData[i].y), sizeof(float));
		}
	}

	//Copy data to device
	msgText("Moving input data to device memory");
	CUDA(cudaMemcpy(d_petalData, p.petalPtr.get(), sizeof(int2) * p.petalRows * p.petalCols, cudaMemcpyHostToDevice));
	CUDA(cudaMemcpy2D(d_inData, d_pitch, h_inData, p.nCols * sizeof(float2), p.nCols * sizeof(float2), nDataRows, cudaMemcpyHostToDevice));
	
	//Copy constants to device
	msgText("Setting device parameters");
	CUDA(cudaMemcpyToSymbol(d_nPetals, &(p.nPetals), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_focalOffset, &(p.fOffset), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_multFocal, &(p.multFocal), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_firstReduction, &(p.firstReduction), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_sumPower, &(p.sumPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_focalPower, &(p.focalPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(p.noData), sizeof(float)));
	
	//Testing dispersal probability
	CUDA(cudaMemcpyToSymbol(d_doAveragePetalHabitat, &(p.doAveragePetalHabitat), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_doDispersalProbability, &(p.doDispersalProbability), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_dispersalYears, &(p.dispersalYears), sizeof(float)));

	//Init Texture Reference
	//msgText("Setting texture reference");
	//tex2Ref.filterMode     = cudaFilterModePoint;
	//tex2Ref.addressMode[0] = cudaAddressModeBorder;
	//tex2Ref.addressMode[1] = cudaAddressModeBorder;
	//cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float2>();
	//CUDA(cudaBindTexture2D (NULL, &tex2Ref, d_inData, &channelDesc, p.nCols, nDataRows, d_pitch));
	
	//Replace Texture Reference with Texture Object
	struct cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypePitch2D;
	resDesc.res.pitch2D.devPtr = d_inData;
	resDesc.res.pitch2D.width = p.nCols;
	resDesc.res.pitch2D.height = nDataRows;
	resDesc.res.pitch2D.desc = cudaCreateChannelDesc<float2>();
	resDesc.res.pitch2D.pitchInBytes = d_pitch;

	struct cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.filterMode = cudaFilterModePoint;
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;

	cudaTextureObject_t tex2Obj;
	cudaCreateTextureObject(&tex2Obj, &resDesc, &texDesc, NULL);
	
	//Profiled Kernel Call Loop
	msgText("Processing data on device");
	Profiler profiler(1000000), kernelTimer(1000000);
	profiler.Start();
	for (i = 0; i < p.nRows; i += gridSize.y) {
		//Set gridSize.y
		gridSize.y = max(min(gridSize.y * maxKernelTime / kernelTime, nBufferRows), 1);//No Linux profiler
		if (i + gridSize.y >= p.nRows) gridSize.y = p.nRows - i;

		//Parallel CBA kernel call
		kernelTimer.Start();
		ContextCBA_kernel <<<gridSize, p.petalCols, p.petalCols * nSharedSlots * sizeof(float) >>>(p.fOffset, 0, d_petalData, d_outBuf, tex2Obj);
		CUDA(cudaDeviceSynchronize());
		kernelTime = (uint)kernelTimer.Stop();
		msgProgress("Percent complete: ", i * 100 / p.nRows);

		//Clear gridSize.y rows in the host buffer then read in any new data
		nReads = (nFirstRows + i) * p.nCols;
		for (j = 0; j < p.nCols * gridSize.y; j++) {
			h_inBuf[j].x = 0.0f;
			h_inBuf[j].y = 0.0f;
			if (nReads + j < p.nCells) {
				p.habInFS.read((char*)&(h_inBuf[j].x), sizeof(float));
				p.prmInFS.read((char*)&(h_inBuf[j].y), sizeof(float));
			}
		}
		//Copy host buffer to device buffer then shuffle into inData
		CUDA(cudaMemcpy2D(d_inBuf, d_pitch, h_inBuf, p.nCols * sizeof(float2), p.nCols * sizeof(float2), gridSize.y, cudaMemcpyHostToDevice));
		ShuffleUp(d_inData, d_inBuf, nDataRows, gridSize.y, int(d_pitch / sizeof(float2)), true);

		//Move device output buffer to host then write to disk
		CUDA(cudaMemcpy(h_outBuf, d_outBuf, p.nCols * gridSize.y * sizeof(float), cudaMemcpyDeviceToHost));
		CUDA(cudaMemset(d_outBuf, 0, p.nCols * gridSize.y * sizeof(float)));
		p.cxtOutFS.write((const char*)h_outBuf, p.nCols * gridSize.y * sizeof(float));

		CUDA(cudaDeviceSynchronize());
	}
	profiler.Stop();
	msgText("\rPercent complete: 100");
	msgText(("Processing time: " + toStr(profiler.Total())).c_str());

	//Free pinned Host Memory and device memory
	msgText("\nFreeing host and device memory");
	CUDA(cudaDestroyTextureObject(tex2Obj));
	CUDA(cudaFreeHost(h_inData));
	CUDA(cudaFreeHost(h_inBuf));
	CUDA(cudaFreeHost(h_outBuf));
	CUDA(cudaFree(d_petalData));
	CUDA(cudaFree(d_inData));
	CUDA(cudaFree(d_inBuf));
	CUDA(cudaFree(d_outBuf));

	//Check CUDA errors
	cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDAContextCBA_L() Complete!");
	return int(cudaStatus);
}


/**/
#endif