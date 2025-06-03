/*
Benefits_CBA.cu - CUDA CBA functions for performing BFT benefits analysis
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

#ifndef _BENEFITS_CBA_KERNEL_CU_
#define _BENEFITS_CBA_KERNEL_CU_
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "Common.cuh"
#include "Parameters.h"
#include "FileUtilsLib.h"

//TODO Explore using grid dim z for processing each petal's alt values
// Calculate Ci first and need memory to store result
//Call nRows x nCols x nPetals threads 


#ifndef USE_PRECALC_EDGE_KERNELS
__global__ void BenefitCBA_kernel(int inOffset, int outOffset, float *d_mbvData, int2 *d_petalData, float *d_outData, cudaTextureObject_t tex4Obj)
{
	//Focal cell
	float4 focal = tex2D<float4>(tex4Obj, blockIdx.x, inOffset + blockIdx.y);
	if ((focal.x == 0.0f && d_multFocal) || focal.y == 0.0f)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = 0.0f;
	else if (
		focal.x == d_noData || focal.x < 0.0 ||
		focal.y == d_noData || focal.y < 0.0 ||
		focal.z == d_noData || focal.z < 0.0 ||
		focal.w == d_noData || focal.w < 0.0 ||
		d_mbvData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] == d_noData
		)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = d_noData;
	//MBV can be less than zero under future scenarios if class EHA increases over its OHA 
	//In this case, we set benefits to 0
	//TODO Is ths the best option....
	else if (d_mbvData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] <= 0.0f)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = 0.0f;
	else {
		//Declare thread local storage
		int i, petal, activeCells = 0;
		int2 tmpPair[4];
		float orig = 0.0f, delta = 0.0f, tmpCst = 0.0f, tmpBen = 0.0f;
		float4 t;

		//Address shared memory
		extern __shared__ float	pBenefit[];
		float *pCost     =		&pBenefit[blockDim.x];
		float *pDist     =		&pBenefit[2 * blockDim.x];
		float *pCount    =	    &pBenefit[3 * blockDim.x];
		pBenefit[threadIdx.x] = 0.0f;
		pCost[threadIdx.x] = 0.0f;
		pDist[threadIdx.x] = 0.0f;
		pCount[threadIdx.x] = 0.0f;

		//Aggregate cell values into petals
		tmpPair[0] = d_petalData[4 * blockDim.x + threadIdx.x];      //petal cell count, valid cell count = 0
		tmpPair[0].y = 0;
		for (i = 0; i < tmpPair[0].x; i++)     {
			tmpPair[1] = d_petalData[(i + 5) * blockDim.x + threadIdx.x];//plus 5 to jump neighbours storage
			if (tmpPair[1].x != -1 && tmpPair[1].y != -1) {
				t = tex2D<float4>(tex4Obj, float(int(blockIdx.x) + tmpPair[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + tmpPair[1].y - d_focalOffset));
				if (t.x >= 0.0f && t.y > 0.0f) {
					pBenefit[threadIdx.x] += t.x;
					pCost[threadIdx.x] += t.y;
					tmpBen += t.z;
					tmpCst += t.w;
					++(tmpPair[0].y);
				}
			}
		}

		//Normalise petal cost and set petal counter
		if (tmpPair[0].y > 0){
			pCost[threadIdx.x] = powf(pCost[threadIdx.x] / (float)tmpPair[0].y, 2.0f * sqrtf((float)tmpPair[0].y * d_oneOnPi));
			tmpCst = powf(tmpCst / (float)tmpPair[0].y, 2.0f * sqrtf((float)tmpPair[0].y * d_oneOnPi));
			pCount[threadIdx.x] = 1.0f;
			activeCells = tmpPair[0].y;
		}

		//Get petal neighbour indices N, S, W, E, NW, NE, SW, SE
		tmpPair[0] = d_petalData[0 * blockDim.x + threadIdx.x];//N, S
		tmpPair[1] = d_petalData[1 * blockDim.x + threadIdx.x];//W, E
		tmpPair[2] = d_petalData[2 * blockDim.x + threadIdx.x];//NW, NE
		tmpPair[3] = d_petalData[3 * blockDim.x + threadIdx.x];//SW, SE

		//Swap values for each petal, petal = -1 is original, petal = d_nPetals is focal cell
		for (petal = -1; petal <= d_nPetals; petal++) {

			//If this thread is the current petal, swap cost and benefit values with alt
			if (petal == threadIdx.x){
				t.x = pBenefit[threadIdx.x];
				t.y = pCost[threadIdx.x];
				pBenefit[threadIdx.x] = tmpBen;
				pCost[threadIdx.x] = tmpCst;
				tmpBen = t.x;
				tmpCst = t.y;
			}
			
			//if petal == nPetals swap focal cell values
			if (petal == d_nPetals) {
				focal.x = focal.z;
				focal.y = powf(focal.w, 2.0f * sqrtf(d_oneOnPi));
			}

			//Reset distance then calculate permeability for eight focal cell neighbours first
			pDist[threadIdx.x] = 0.0f;
			if (threadIdx.x < 8)
				pDist[threadIdx.x] = threadIdx.x % 2 == 0 ? powf(pCost[threadIdx.x] * focal.y, d_diag) : sqrtf(pCost[threadIdx.x] * focal.y);

			//Loop until all permeabilities are maximised (i==0)
			i = 1;
			while (__syncthreads_or(i)) {
				i = 0;
				if (tmpPair[0].x != -1 && pDist[tmpPair[0].x] * sqrtf(pCost[threadIdx.x] * pCost[tmpPair[0].x]) > pDist[threadIdx.x]) {//N                
					pDist[threadIdx.x] = pDist[tmpPair[0].x] * sqrtf(pCost[threadIdx.x] * pCost[tmpPair[0].x]);
					i++;
				}
				if (tmpPair[0].y != -1 && pDist[tmpPair[0].y] * sqrtf(pCost[threadIdx.x] * pCost[tmpPair[0].y]) > pDist[threadIdx.x]) {//S                
					pDist[threadIdx.x] = pDist[tmpPair[0].y] * sqrtf(pCost[threadIdx.x] * pCost[tmpPair[0].y]);
					i++;
				}
				if (tmpPair[1].x != -1 && pDist[tmpPair[1].x] * sqrtf(pCost[threadIdx.x] * pCost[tmpPair[1].x]) > pDist[threadIdx.x]) {//E                
					pDist[threadIdx.x] = pDist[tmpPair[1].x] * sqrtf(pCost[threadIdx.x] * pCost[tmpPair[1].x]);
					i++;
				}
				if (tmpPair[1].y != -1 && pDist[tmpPair[1].y] * sqrtf(pCost[threadIdx.x] * pCost[tmpPair[1].y]) > pDist[threadIdx.x]) {//W                
					pDist[threadIdx.x] = pDist[tmpPair[1].y] * sqrtf(pCost[threadIdx.x] * pCost[tmpPair[1].y]);
					i++;
				}
				if (tmpPair[2].x != -1 && pDist[tmpPair[2].x] * powf(pCost[threadIdx.x] * pCost[tmpPair[2].x], d_diag) > pDist[threadIdx.x]) {//NW                
					pDist[threadIdx.x] = pDist[tmpPair[2].x] * powf(pCost[threadIdx.x] * pCost[tmpPair[2].x], d_diag);
					i++;
				}
				if (tmpPair[2].y != -1 && pDist[tmpPair[2].y] * powf(pCost[threadIdx.x] * pCost[tmpPair[2].y], d_diag) > pDist[threadIdx.x]) {//NE               
					pDist[threadIdx.x] = pDist[tmpPair[2].y] * powf(pCost[threadIdx.x] * pCost[tmpPair[2].y], d_diag);
					i++;
				}
				if (tmpPair[3].x != -1 && pDist[tmpPair[3].x] * powf(pCost[threadIdx.x] * pCost[tmpPair[3].x], d_diag) > pDist[threadIdx.x]) {//SW               
					pDist[threadIdx.x] = pDist[tmpPair[3].x] * powf(pCost[threadIdx.x] * pCost[tmpPair[3].x], d_diag);
					i++;
				}
				if (tmpPair[3].y != -1 && pDist[tmpPair[3].y] * powf(pCost[threadIdx.x] * pCost[tmpPair[3].y], d_diag) > pDist[threadIdx.x]) {//SE               
					pDist[threadIdx.x] = pDist[tmpPair[3].y] * powf(pCost[threadIdx.x] * pCost[tmpPair[3].y], d_diag);
					i++;
				}
			}


			//Calculate HjWij
			pDist[threadIdx.x] *= pBenefit[threadIdx.x];
			__syncthreads();

			//Sum all HjWij and petal count
			//https://docs.nvidia.com/cuda/samples/6_Advanced/reduction/doc/reduction.pdf
			for (i = d_firstReduction; i > 0; i /= 2) {
				if (threadIdx.x < i && threadIdx.x + i < d_nPetals) {
					pDist[threadIdx.x] += pDist[threadIdx.x + i];
					pCount[threadIdx.x] += pCount[threadIdx.x + i];
				}
				__syncthreads();
			}

			//Final focal cell calculations
			if (threadIdx.x == 0) {
				pDist[0] += focal.x;

				//pDist[0] *= (pCount[0] < (float)d_nPetals) ? ((float)d_nPetals) / (pCount[0] + 1.0f) : 1.0f;
				//pDist[0] *= d_multFocal ? powf(focal.x, d_focalPower) : 1.0f;
				//pDist[0] *= d_mbvData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x];
				pDist[0] = powf(pDist[0], d_sumPower);
				pDist[0] *= (focal.x * d_mbvData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x]);
			}

			__syncthreads();

			//If petal == -1 then all threads calculate their original value
			if (petal == -1) orig = pDist[0];

			//Else if petal == nPetals then thread 0 adds delta to focal cell in outData
			else if (threadIdx.x == 0 && petal == d_nPetals) {
				d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] += (pDist[0] - orig);
			}

			//Else if this is the current petal then save delta and swap back cost and benefit values
			else if (threadIdx.x == petal && petal < d_nPetals && activeCells > 0) {
				delta = (pDist[0] - orig) / activeCells;
				t.x = pBenefit[threadIdx.x];
				t.y = pCost[threadIdx.x];
				pBenefit[threadIdx.x] = tmpBen;
				pCost[threadIdx.x] = tmpCst;
				tmpBen = t.x;
				tmpCst = t.y;
			}
			__syncthreads();
		}//for petal

		//Add delta to all petal cells in outData
		tmpPair[0] = d_petalData[4 * blockDim.x + threadIdx.x];
		for (i = 0; i < tmpPair[0].x; i++) {
			tmpPair[1] = d_petalData[(i + 5) * blockDim.x + threadIdx.x];
			if (tmpPair[1].x != -1 && tmpPair[1].y != -1) {
				t = tex2D<float4>(tex4Obj, float(int(blockIdx.x) + tmpPair[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + tmpPair[1].y - d_focalOffset));
				if (t.y > 0.0f) {
						atomicAdd(&(d_outData[(outOffset + blockIdx.y + tmpPair[1].y - d_focalOffset) * gridDim.x + (blockIdx.x + tmpPair[1].x - d_focalOffset)]), delta);
				}
			}
		}
	}
}

#else
//NOTE Do NOT use causes non - reproducable errors
__global__ void BenefitCBA_kernel(int inOffset, int outOffset, float *d_mbvData, int2 *d_petalData, float *d_outData)
{
	//Focal cell
	float4 focal = tex2D(tex4Ref, blockIdx.x, inOffset + blockIdx.y);
	if (focal.y == 0.0f)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = 0.0f;
	else if (focal.x == d_noData || focal.y == d_noData || focal.z == d_noData || focal.w == d_noData)
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] = d_noData;
	else {

		//Declare thread local storage
		int i, petal;
		float Hj = 0.0f, orig = 0.0f, delta = 0.0f;
		int2 tmpPair[4];
		float4 costNSWE{ -1.0f, -1.0f, -1.0f, -1.0f }; //also used for fetching cell values
		float4 costNDSD{ -1.0f, -1.0f, -1.0f, -1.0f };

		//Address shared memory
		extern __shared__ float	pCost[];
		float *pDist = &pCost[blockDim.x];
		float *pAltHab = &pCost[2 * blockDim.x];
		float *pAltCst = &pCost[3 * blockDim.x];
		pCost[threadIdx.x] = 0.0f;
		pDist[threadIdx.x] = 0.0f;
		pAltHab[threadIdx.x] = 0.0f;
		pAltCst[threadIdx.x] = 0.0f;

		//Calculate petal values from cells
		tmpPair[0] = d_petalData[4 * blockDim.x + threadIdx.x];      //petal cell count, valid cell count = 0
		tmpPair[0].y = 0;
		for (i = 0; i < tmpPair[0].x; i++) {
			tmpPair[1] = d_petalData[(i + 5) * blockDim.x + threadIdx.x];//plus 5 to jump neighbours storage
			if (tmpPair[1].x != -1 && tmpPair[1].y != -1) {
				costNSWE = tex2D(tex4Ref, float(int(blockIdx.x) + tmpPair[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + tmpPair[1].y - d_focalOffset));
				if (costNSWE.y > d_fltEps) {
					Hj += costNSWE.x;
					pCost[threadIdx.x] += costNSWE.y;
					pAltHab[threadIdx.x] += costNSWE.z;
					pAltCst[threadIdx.x] += costNSWE.w;
					++(tmpPair[0].y);
				}
			}
		}

		//Normalise petal cost and set petal counter
		if (tmpPair[0].y > 0) {
			pCost[threadIdx.x] = powf(pCost[threadIdx.x] / (float)tmpPair[0].y, 2.0f * sqrtf((float)tmpPair[0].y * d_oneOnPi));
			pAltCst[threadIdx.x] = powf(pAltCst[threadIdx.x] / (float)tmpPair[0].y, 2.0f * sqrtf((float)tmpPair[0].y * d_oneOnPi));
		}

		//Get petal neighbour indices
		tmpPair[0] = d_petalData[0 * blockDim.x + threadIdx.x];//N, S
		tmpPair[1] = d_petalData[1 * blockDim.x + threadIdx.x];//W, E
		tmpPair[2] = d_petalData[2 * blockDim.x + threadIdx.x];//NW, NE
		tmpPair[3] = d_petalData[3 * blockDim.x + threadIdx.x];//SW, SE

															   //For each petal swap values if petal. -1 is original, petal = d_nPetals is focal cell
		for (petal = -1; petal <= d_nPetals; petal++) {

			//Reset i
			i = 1;

			//If this is the current petal swap cost and benefit values with alt
			if (petal == threadIdx.x) {
				costNSWE.y = pCost[threadIdx.x];
				pCost[threadIdx.x] = pAltCst[threadIdx.x];
				pAltCst[threadIdx.x] = costNSWE.y;
			}

			//if petal == nPetals swap focal values
			if (petal == d_nPetals) {
				focal.x = focal.z;
				focal.y = powf(focal.w, 2.0f * sqrtf(d_oneOnPi));
			}

			//Calculate petal edge costs
			__syncthreads();
			costNSWE.x = tmpPair[0].x != -1 ? sqrtf(pCost[threadIdx.x] * pCost[tmpPair[0].x]) : -1.0f;
			costNSWE.y = tmpPair[0].y != -1 ? sqrtf(pCost[threadIdx.x] * pCost[tmpPair[0].y]) : -1.0f;
			costNSWE.z = tmpPair[1].x != -1 ? sqrtf(pCost[threadIdx.x] * pCost[tmpPair[1].x]) : -1.0f;
			costNSWE.w = tmpPair[1].y != -1 ? sqrtf(pCost[threadIdx.x] * pCost[tmpPair[1].y]) : -1.0f;
			costNDSD.x = tmpPair[2].x != -1 ? powf(pCost[threadIdx.x] * pCost[tmpPair[2].x], d_diag) : -1.0f;
			costNDSD.y = tmpPair[2].y != -1 ? powf(pCost[threadIdx.x] * pCost[tmpPair[2].y], d_diag) : -1.0f;
			costNDSD.z = tmpPair[3].x != -1 ? powf(pCost[threadIdx.x] * pCost[tmpPair[3].x], d_diag) : -1.0f;
			costNDSD.w = tmpPair[3].y != -1 ? powf(pCost[threadIdx.x] * pCost[tmpPair[3].y], d_diag) : -1.0f;

			//Calculate permeability for focal cell neighbours first
			if (threadIdx.x < 8)
				pDist[threadIdx.x] = threadIdx.x % 2 == 0 ? powf(pCost[threadIdx.x] * focal.y, d_diag) : sqrtf(pCost[threadIdx.x] * focal.y);
			else
				pDist[threadIdx.x] = 0.0f;

			//Loop until all permeabilities are maximised (i==0) then multiply by Hj 
			i = 1;
			while (__syncthreads_or(i)) {
				i = 0;
				if (tmpPair[0].x != -1 && pDist[tmpPair[0].x] * costNSWE.x > pDist[threadIdx.x]) {//N edge
					pDist[threadIdx.x] = pDist[tmpPair[0].x] * costNSWE.x;
					i = 1;
				}
				if (tmpPair[0].y != -1 && pDist[tmpPair[0].y] * costNSWE.y > pDist[threadIdx.x]) {//S edge
					pDist[threadIdx.x] = pDist[tmpPair[0].y] * costNSWE.y;
					i = 1;
				}
				if (tmpPair[1].x != -1 && pDist[tmpPair[1].x] * costNSWE.z > pDist[threadIdx.x]) {//W edge
					pDist[threadIdx.x] = pDist[tmpPair[1].x] * costNSWE.z;
					i = 1;
				}
				if (tmpPair[1].y != -1 && pDist[tmpPair[1].y] * costNSWE.w > pDist[threadIdx.x]) {//E edge
					pDist[threadIdx.x] = pDist[tmpPair[1].y] * costNSWE.w;
					i = 1;
				}
				if (tmpPair[2].x != -1 && pDist[tmpPair[2].x] * costNDSD.x > pDist[threadIdx.x]) {//NW edge
					pDist[threadIdx.x] = pDist[tmpPair[2].x] * costNDSD.x;
					i = 1;
				}
				if (tmpPair[2].y != -1 && pDist[tmpPair[2].y] * costNDSD.y > pDist[threadIdx.x]) {//NE edge
					pDist[threadIdx.x] = pDist[tmpPair[2].y] * costNDSD.y;
					i = 1;
				}
				if (tmpPair[3].x != -1 && pDist[tmpPair[3].x] * costNDSD.z > pDist[threadIdx.x]) {//SW edge
					pDist[threadIdx.x] = pDist[tmpPair[3].x] * costNDSD.z;
					i = 1;
				}
				if (tmpPair[3].y != -1 && pDist[tmpPair[3].y] * costNDSD.w > pDist[threadIdx.x]) {//SE edge
					pDist[threadIdx.x] = pDist[tmpPair[3].y] * costNDSD.w;
					i = 1;
				}
			}
			pDist[threadIdx.x] *= (petal == threadIdx.x ? pAltHab[threadIdx.x] : Hj);
			__syncthreads();

			//Sum all HjWij and petal count
			for (i = d_firstReduction; i > 0; i /= 2) {
				if (threadIdx.x < i && threadIdx.x + i < d_nPetals) {
					pDist[threadIdx.x] += pDist[threadIdx.x + i];
				}
				__syncthreads();
			}

			//Final focal cell calcs
			if (threadIdx.x == 0) {
				pDist[0] += focal.x;
				pDist[0] = powf(pDist[0], d_sumPower);
				pDist[0] *= (focal.x * d_mbvData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x]);
			}
			__syncthreads();

			//If petal == -1 all threads calculate original value
			if (petal == -1) orig = pDist[0];

			//If petal == nPetals thread 0 writes out focal cell
			else if (threadIdx.x == 0 && petal == d_nPetals) {
				d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x] += (pDist[0] - orig);
			}

			//If this is the current petal, save delta and swap back cost value
			else if (threadIdx.x == petal && petal < d_nPetals) {
				delta = (pDist[0] - orig) / d_petalData[4 * blockDim.x + threadIdx.x].x;
				costNSWE.y = pCost[threadIdx.x];
				pCost[threadIdx.x] = pAltCst[threadIdx.x];
				pAltCst[threadIdx.x] = costNSWE.y;
			}
			__syncthreads();
		}//for petal

		 //Add delta to each cell in the petal
		tmpPair[0] = d_petalData[4 * blockDim.x + threadIdx.x];
		for (i = 0; i < tmpPair[0].x; i++) {
			tmpPair[1] = d_petalData[(i + 5) * blockDim.x + threadIdx.x];
			if (tmpPair[1].x != -1 && tmpPair[1].y != -1) {
				costNSWE = tex2D(tex4Ref, float(int(blockIdx.x) + tmpPair[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + tmpPair[1].y - d_focalOffset));
				if (costNSWE.y > d_fltEps) {
					atomicAdd(&(d_outData[(outOffset + blockIdx.y + tmpPair[1].y - d_focalOffset) * gridDim.x + (blockIdx.x + tmpPair[1].x - d_focalOffset)]), delta);
				}
			}
		}
	}
}
#endif

////////////////////////////////////////////////////////////
//Benefit CBA small
////////////////////////////////////////////////////////////
int CUDABenefitCBA_S(CBAParams &p)
{
	cudaError_t cudaStatus;
	msgText("Performing Benefits CBA");

	//Host and Device data
	float4 *h_inData, *d_inData;
	float *h_mbvData, *d_mbvData;
	float *h_outData, *d_outData;
	int2 *d_petalData;

	//Local variables
	size_t d_pitch;
	dim3 gridSize{ (unsigned int)p.nCols, 1U, 1U };
	uint maxKernelTime = h_maxKernelTime;
	uint kernelTime = maxKernelTime;
	uint i;

	//Host 6 x 4bytes x nCells, Device 6 x 4bytes x nCells
	//Malloc pinned host memory
	msgText("Allocating host memory");
	CUDA(cudaMallocHost(&h_inData, p.nCells * sizeof(float4)));
	CUDA(cudaMallocHost(&h_mbvData, p.nCells * sizeof(float)));
	CUDA(cudaMallocHost(&h_outData, p.nCells * sizeof(float)));

	//Malloc device data
	msgText("Allocating device memory");
	CUDA(cudaMalloc(&d_petalData, sizeof(int2) * p.petalRows * p.petalCols));
	CUDA(cudaMallocPitch(&d_inData, &d_pitch, p.nCols * sizeof(float4), p.nRows));
	CUDA(cudaMalloc(&d_mbvData, p.nCells * sizeof(float)));
	CUDA(cudaMalloc(&d_outData, p.nCells * sizeof(float)));

	//Read h_inData from disk
	msgText("Reading input data to host memory");
	if (p.useAltGrids) {
		for (i = 0; i < p.nCells; i++) {
			p.habInFS.read((char*)&(h_inData[i].x), sizeof(float));
			p.prmInFS.read((char*)&(h_inData[i].y), sizeof(float));
			p.altHabInFN.read((char*)&(h_inData[i].z), sizeof(float));
			p.altPrmInFS.read((char*)&(h_inData[i].w), sizeof(float));
			h_outData[i] = 0.0f;
		}
	}else{
		for (i = 0; i < p.nCells; i++) {
			p.habInFS.read((char*)&(h_inData[i].x), sizeof(float));
			p.prmInFS.read((char*)&(h_inData[i].y), sizeof(float));
			h_inData[i].z = p.altHab;
			h_inData[i].w = p.altPrm;
			h_outData[i] = 0.0f;
		}
	}
	//read mbv data
	p.mbvInFS.read((char*)h_mbvData, p.nCells * sizeof(float));
	
	//Mask all grids nodata
	//TODO Probably not needed since kernel checks for nodata and neg values
	//for (i = 0; i < p.nCells; i++) {
	//	if (	h_inData[i].x == p.noData || h_inData[i].x < 0.0f ||
	//			h_inData[i].y == p.noData || h_inData[i].y < 0.0f ||
	//			h_inData[i].z == p.noData || h_inData[i].z < 0.0f ||
	//			h_inData[i].w == p.noData || h_inData[i].w < 0.0f ||
	//			h_mbvData[i] == p.noData || h_mbvData[i] < 0.0f)
	//	{
	//		h_inData[i].x = p.noData;
	//		h_inData[i].y = p.noData;
	//		h_inData[i].z = p.noData;
	//		h_inData[i].w = p.noData;
	//		h_mbvData[i] = p.noData;
	//	}
	//}

	//Move input data to device
	msgText("Moving input data to device memory");
	CUDA(cudaMemcpy(d_petalData, p.petalPtr.get(), sizeof(int2) * p.petalRows * p.petalCols, cudaMemcpyHostToDevice));
	CUDA(cudaMemcpy2D(d_inData, d_pitch, h_inData, p.nCols * sizeof(float4), p.nCols * sizeof(float4), p.nRows, cudaMemcpyHostToDevice));
	CUDA(cudaMemcpy(d_mbvData, h_mbvData, sizeof(float) * p.nCells, cudaMemcpyHostToDevice));
	CUDA(cudaMemcpy(d_outData, h_outData, sizeof(float) * p.nCells, cudaMemcpyHostToDevice));
	CUDA(cudaFreeHost(h_inData));
	CUDA(cudaFreeHost(h_mbvData));

	//Copy constants to device
	msgText("Setting device parameters");
	CUDA(cudaMemcpyToSymbol(d_nPetals, &(p.nPetals), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_focalOffset, &(p.fOffset), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_multFocal, &(p.multFocal), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_firstReduction, &(p.firstReduction), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_sumPower, &(p.sumPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_focalPower, &(p.focalPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(p.noData), sizeof(float)));

	//Texture object replaces texture reference 
	struct cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypePitch2D;
	resDesc.res.pitch2D.desc = cudaCreateChannelDesc<float4>();
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

	cudaTextureObject_t tex4Obj = 0;
	CUDA(cudaCreateTextureObject(&tex4Obj, &resDesc, &texDesc, NULL));

	//Profiled Kernel Call Loop
	msgText("Processing data on device");
	Profiler profiler(1000000), kernelTimer(1000000);
	profiler.Start();
	for (i = 0; i < p.nRows; i += gridSize.y) {
		//Set gridSize.y
		gridSize.y = max(gridSize.y * maxKernelTime / kernelTime, 1);
		gridSize.y = min(gridSize.y, p.nRows - i);

		//Parallel CBA kernel call
		kernelTimer.Start();
		BenefitCBA_kernel << <gridSize, p.petalCols, p.petalCols * 4 * sizeof(float) >> >(i, i, d_mbvData, d_petalData, d_outData, tex4Obj);
		CUDA(cudaDeviceSynchronize());
		kernelTime = (uint)kernelTimer.Stop();
		msgProgress("Percent complete: ", i * 100 / p.nRows);
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
	CUDA(cudaDestroyTextureObject(tex4Obj));
	CUDA(cudaFreeHost(h_outData));
	CUDA(cudaFree(d_petalData));
	CUDA(cudaFree(d_inData));
	CUDA(cudaFree(d_mbvData));
	CUDA(cudaFree(d_outData));

	//Check CUDA errors
	cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDABenefitCBA_S() Complete!");
	return int(cudaStatus);
}

////////////////////////////////////////////////////////////
//Benefit CBA large
////////////////////////////////////////////////////////////
int CUDABenefitCBA_L(CBAParams &p)
{
	cudaError_t cudaStatus;
	msgText("Performing Benefits CBA");

	//Host and Device data
	float4 *h_inData, *d_inData;
	float4 *h_inBuf, *d_inBuf;
	float *h_mbvData, *d_mbvData;
	float *h_mbvBuf, *d_mbvBuf;
	float *d_outData;
	float *h_outBuf, *d_outBuf;
	int2 *d_petalData;

	size_t d_pitch;
	uint nBufferRows = p.fOffset * 2;
	uint nDataRows = p.fOffset * 4;
	uint nFirstRows = nDataRows - p.fOffset;
	ull nReads = 0;
	dim3 gridSize{ (unsigned int) p.nCols, 1U, 1U };
	uint maxKernelTime = h_maxKernelTime;
	uint kernelTime = maxKernelTime;
	uint i, j;

	//Malloc pinned host memory
	msgText("Allocating host memory");
	CUDA(cudaMallocHost(&h_inData, p.nCols * nDataRows * sizeof(float4)));
	CUDA(cudaMallocHost(&h_inBuf, p.nCols * nBufferRows * sizeof(float4)));
	CUDA(cudaMallocHost(&h_mbvData, p.nCols * nDataRows * sizeof(float)));
	CUDA(cudaMallocHost(&h_mbvBuf, p.nCols * nBufferRows * sizeof(float)));
	CUDA(cudaMallocHost(&h_outBuf, p.nCols * nDataRows * sizeof(float)));

	//Malloc device data
	msgText("Allocating device memory");
	CUDA(cudaMalloc(&d_petalData, sizeof(int2) * p.petalRows * p.petalCols));
	CUDA(cudaMallocPitch(&d_inData, &d_pitch, p.nCols * sizeof(float4), nDataRows));
	CUDA(cudaMallocPitch(&d_inBuf, &d_pitch, p.nCols * sizeof(float4), nBufferRows));
	CUDA(cudaMalloc(&d_mbvData, p.nCols * nDataRows * sizeof(float)));
	CUDA(cudaMalloc(&d_mbvBuf, p.nCols * nBufferRows * sizeof(float)));
	CUDA(cudaMalloc(&d_outData, p.nCols * nDataRows * sizeof(float)));
	CUDA(cudaMalloc(&d_outBuf, p.nCols * nBufferRows * sizeof(float)));

	//Read h_inData from disk
	msgText("Reading input data to host memory");
	if (p.useAltGrids) {
		for (i = 0; i < nDataRows * p.nCols; i++) {
			h_inData[i].x = 0;
			h_inData[i].y = 0;
			h_inData[i].z = 0;
			h_inData[i].w = 0;
			h_mbvData[i] = 0;
			if (i >= p.fOffset * p.nCols) {
				p.habInFS.read((char*)&(h_inData[i].x), sizeof(float));
				p.prmInFS.read((char*)&(h_inData[i].y), sizeof(float));
				p.altHabInFN.read((char*)&(h_inData[i].z), sizeof(float));
				p.altPrmInFS.read((char*)&(h_inData[i].w), sizeof(float));
				p.mbvInFS.read((char*)&(h_mbvData[i]), sizeof(float));
			}
		}
	}
	else {
		for (i = 0; i < nDataRows * p.nCols; i++) {
			h_inData[i].x = 0;
			h_inData[i].y = 0;
			h_inData[i].z = 0;
			h_inData[i].w = 0;
			h_mbvData[i] = 0;
			if (i >= p.fOffset * p.nCols) {
				p.habInFS.read((char*)&(h_inData[i].x), sizeof(float));
				p.prmInFS.read((char*)&(h_inData[i].y), sizeof(float));
				h_inData[i].z = h_inData[i].x < 0.0f ? p.noData : p.altHab;
				h_inData[i].w = h_inData[i].y < 0.0f ? p.noData : p.altPrm;
				p.mbvInFS.read((char*)&(h_mbvData[i]), sizeof(float));
			}
		}
	}

	//Copy data to device
	msgText("Moving input data to device memory");
	CUDA(cudaMemcpy(d_petalData, p.petalPtr.get(), sizeof(int2) * p.petalRows * p.petalCols, cudaMemcpyHostToDevice));
	CUDA(cudaMemcpy2D(d_inData, d_pitch, h_inData, p.nCols * sizeof(float4), p.nCols * sizeof(float4), nDataRows, cudaMemcpyHostToDevice));
	CUDA(cudaMemcpy(d_mbvData, h_mbvData, p.nCols * nDataRows * sizeof(float), cudaMemcpyHostToDevice));
	CUDA(cudaMemset(d_outData, 0, p.nCols * nDataRows * sizeof(float)));
	CUDA(cudaMemset(d_outBuf, 0, p.nCols * nBufferRows * sizeof(float)));
	CUDA(cudaFreeHost(h_inData));
	CUDA(cudaFreeHost(h_mbvData));

	//Copy constants to device
	msgText("Setting device parameters");
	CUDA(cudaMemcpyToSymbol(d_nPetals, &(p.nPetals), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_focalOffset, &(p.fOffset), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_multFocal, &(p.multFocal), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_firstReduction, &(p.firstReduction), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_sumPower, &(p.sumPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_focalPower, &(p.focalPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(p.noData), sizeof(float)));

	//Texture object replaces texture reference 
	struct cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypePitch2D;
	resDesc.res.pitch2D.desc = cudaCreateChannelDesc<float4>();
	resDesc.res.pitch2D.devPtr = d_inData;
	resDesc.res.pitch2D.width = p.nCols;
	resDesc.res.pitch2D.height = nDataRows;
	resDesc.res.pitch2D.pitchInBytes = d_pitch;

	struct cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.filterMode = cudaFilterModePoint;
	texDesc.addressMode[0] = cudaAddressModeBorder;
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.readMode = cudaReadModeElementType;

	cudaTextureObject_t tex4Obj = 0;
	CUDA(cudaCreateTextureObject(&tex4Obj, &resDesc, &texDesc, NULL));

	//Profiled Kernel Call Loop
	msgText("Processing data on device");
	Profiler profiler(1000000), kernelTimer(1000000);
	profiler.Start();
	for (i = 0; i < p.nRows; i += gridSize.y) {
		//Set gridSize.y
		gridSize.y = max(min(gridSize.y * maxKernelTime / kernelTime, nBufferRows), 1);
		gridSize.y = min(gridSize.y, p.nRows - i);

		//Parallel CBA kernel call
		kernelTimer.Start();
		BenefitCBA_kernel << <gridSize, p.petalCols, p.petalCols * 4 * sizeof(float) >> > (p.fOffset, p.fOffset, d_mbvData, d_petalData, d_outData, tex4Obj);
		CUDA(cudaDeviceSynchronize());
		kernelTime = (uint)kernelTimer.Stop();
		msgProgress("Percent complete: ", i * 100 / p.nRows);

		//Clear gridSize.y rows in the host buffer then read in any new data
		nReads = (nFirstRows + i) * p.nCols;
		for (j = 0; j < p.nCols * gridSize.y; j++) {
			h_inBuf[j].x = 0.0f;
			h_inBuf[j].y = 0.0f;
			h_inBuf[j].z = 0.0f;
			h_inBuf[j].w = 0.0f;
			h_mbvBuf[j] = 0.0f;
			if (nReads + j < p.nCells) {
				p.habInFS.read((char*)&(h_inBuf[j].x), sizeof(float));
				p.prmInFS.read((char*)&(h_inBuf[j].y), sizeof(float));
				p.mbvInFS.read((char*)&h_mbvBuf[j], sizeof(float));
				if (p.useAltGrids) {
					p.altHabInFN.read((char*)&(h_inBuf[j].z), sizeof(float));
					p.altPrmInFS.read((char*)&(h_inBuf[j].w), sizeof(float));
				} else {
					h_inBuf[j].z = h_inBuf[j].x < 0.0f ? p.noData : p.altHab;
					h_inBuf[j].w = h_inBuf[j].y < 0.0f ? p.noData : p.altPrm;
				}
			}
		}
		//Copy host buffer to device buffer then shuffle into inData
		CUDA(cudaMemcpy2D(d_inBuf, d_pitch, h_inBuf, p.nCols * sizeof(float4), p.nCols * sizeof(float4), gridSize.y, cudaMemcpyHostToDevice));
		CUDA(cudaMemcpy(d_mbvBuf, h_mbvBuf, p.nCols * gridSize.y * sizeof(float), cudaMemcpyHostToDevice));
		ShuffleUp(d_inData, d_inBuf, nDataRows, gridSize.y, int(d_pitch / sizeof(float4)), true);
		ShuffleUp(d_mbvData, d_mbvBuf, nDataRows, gridSize.y, p.nCols, true);

		//Move device output buffer to host then write to disk
		if (i < p.fOffset) {//at start so don't write out the head of outData
			if (i + gridSize.y > p.fOffset) {
				int nRows = gridSize.y - p.fOffset + i;
				CUDA(cudaMemcpy(h_outBuf, &(d_outData[(p.fOffset - i) * p.nCols]), p.nCols * nRows * sizeof(float), cudaMemcpyDeviceToHost));
				p.cxtOutFS.write((const char*)h_outBuf, p.nCols * nRows * sizeof(float));
			}
			ShuffleUp(d_outData, d_outBuf, nDataRows, gridSize.y, p.nCols, true);
		}
		else if (i + gridSize.y < p.nRows) {//middle so write out the first nBufferRows in outData
			CUDA(cudaMemcpy(h_outBuf, d_outData, p.nCols * gridSize.y * sizeof(float), cudaMemcpyDeviceToHost));
			p.cxtOutFS.write((const char*)h_outBuf, p.nCols * gridSize.y * sizeof(float));
			ShuffleUp(d_outData, d_outBuf, nDataRows, gridSize.y, p.nCols, true);
		}
		else {//end so write out the nBufferRows and the tail of outData
			CUDA(cudaMemcpy(h_outBuf, d_outData, p.nCols * (gridSize.y + p.fOffset) * sizeof(float), cudaMemcpyDeviceToHost));
			p.cxtOutFS.write((const char*)h_outBuf, p.nCols * (gridSize.y + p.fOffset) * sizeof(float));
		}
		CUDA(cudaDeviceSynchronize());
	}

	profiler.Stop();
	msgText("\rPercent complete: 100");
	msgText(("Processing time: " + toStr(profiler.Total())).c_str());

	//Free pinned host and device memory
	msgText("Freeing host and device memory");
	CUDA(cudaDestroyTextureObject(tex4Obj));
	CUDA(cudaFreeHost(h_inBuf));
	CUDA(cudaFreeHost(h_mbvBuf));
	CUDA(cudaFreeHost(h_outBuf));
	CUDA(cudaFree(d_petalData));
	CUDA(cudaFree(d_inData));
	CUDA(cudaFree(d_inBuf));
	CUDA(cudaFree(d_mbvData));
	CUDA(cudaFree(d_mbvBuf));
	CUDA(cudaFree(d_outData));
	CUDA(cudaFree(d_outBuf));

	//Check CUDA errors
	cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDABenefitCBA_L() Complete!");
	return int(cudaStatus);
}

#endif