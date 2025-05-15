/*
REMP_CBA.cu - CUDA CBA functions for performing REMP analysis
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


/*
* Resolving differences between _S and _L versions

CUDAOccupancyCBA_TS_S()

for (uint t = 0; t < p.rempTimeSteps; t++)
	h_inData[i].w = p.mpcInitFactor
	if (p.haveMPCIn && t == 0)
		p.mpcInFS.read((char*)&(h_inData[i].w), sizeof(float));

	for (uint remp_i = 0; remp_i < p.rempIterations; remp_i++) {

		OccupancyCBA_kernel()

		p.mpcDen = max(h_outData[j].y, p.mpcDen);

		if (remp_i + 1 < p.rempIterations)
			//Should this happen for all iterations since last iteration does h_outData[j].y / h_inData[j].w
			h_inData[j].w = h_outData[j].y > fltEps ? h_outData[j].y / p.mpcDen : h_outData[j].y;

		if (remp_i + 1 < p.rempIterations && p.rempWriteEachIteration) {
			p.mpcOutFS.open(GridFN(p.mpcOutFNs[t] + "_r_i" + toStr(remp_i)), std::ios::out | std::ios::binary);
			p.mpcOutFS.write((const char*)&(h_outData[j].y), sizeof(float));

	if (p.haveMPCOut) {
		p.mpcOutFS.open(GridFN(p.mpcOutFNs[t]), std::ios::out | std::ios::binary);
		h_outData[j].y = h_inData[j].w > fltEps ? h_outData[j].y / h_inData[j].w : h_outData[j].y;
		p.mpcOutFS.write((const char*)&(h_outData[j].y), sizeof(float));

CUDAOccupancyCBA_TS_S()


* 
*/




#ifndef _REMP_CBA_KERNEL_CU_
#define _REMP_CBA_KERNEL_CU_
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
__global__ void OccupancyCBA_kernel(int inOffset, int outOffset, int2 *d_petalData, float2 *d_outData, cudaTextureObject_t tex4Obj)
{
	//Focal cell
	float4 focal = tex2D<float4>(tex4Obj, blockIdx.x, blockIdx.y + inOffset);
	 if (focal.x == d_noData || focal.y == d_noData) {
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].x = d_noData;
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].y = d_noData;
	} else if (focal.y <= 0.0f) {
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].x = 0.0f;
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].y = 0.0f;
	} else {
		//Declare thread local storage
		int i;
		int2 intPairs[4];
		float4 t;

		//Address shared memory
		extern __shared__ float pOcc[];         //Occupancy
		float *pMpc = &pOcc[blockDim.x];		//mpc
		float *pCost = &pOcc[2 * blockDim.x];   //cost
		float *pDist = &pOcc[3 * blockDim.x];   //Dist
		float *pCount = &pOcc[4 * blockDim.x];  //count row
		pOcc[threadIdx.x] = 0.0f;
		pMpc[threadIdx.x] = 0.0f;
		pCost[threadIdx.x] = 0.0f;
		pDist[threadIdx.x] = 0.0f;
		pCount[threadIdx.x] = 0.0f;

		//Aggregate cell values into petals
		intPairs[0] = d_petalData[4 * blockDim.x + threadIdx.x];//petal cell count, valid cell count = 0
		intPairs[0].y = 0;
		for (i = 0; i < intPairs[0].x; i++)     {
			intPairs[1] = d_petalData[(i + 5) * blockDim.x + threadIdx.x];//plus 5 to jump neighbours storage
			if (intPairs[1].x != -1 && intPairs[1].y != -1) {
				t = tex2D<float4>(tex4Obj, float(int(blockIdx.x) + intPairs[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + intPairs[1].y - d_focalOffset));
				if (t.x >= 0.0f && t.y > 0.0f) {
					pOcc[threadIdx.x] += (t.x * t.z);//Occ
					pMpc[threadIdx.x] += (t.x * t.w);//Mpc

					//Use average or maximum permeability in petal
					if (d_useMaxPetalPerm)
						//pCost[threadIdx.x] = __max(pCost[threadIdx.x], t.y); identifier "__max" is undefined
						pCost[threadIdx.x] = max(pCost[threadIdx.x], t.y);
					else
						pCost[threadIdx.x] += t.y;
	
					++(intPairs[0].y);
				}
			}
		}

		//Normalise petal cost and set petal counter
		if (intPairs[0].y > 0){
			if (d_useMaxPetalPerm)
				pCost[threadIdx.x] = powf(pCost[threadIdx.x], 2.0f * sqrtf(float(intPairs[0].y) * d_oneOnPi));
			else
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

		//Multiply Occ and MPC by permeability
		pOcc[threadIdx.x] *= pDist[threadIdx.x];
		pMpc[threadIdx.x] *= pDist[threadIdx.x];
		__syncthreads();

		//Sum all Occ, MPC and petal count
		for (i = d_firstReduction; i > 0; i /= 2) {
			if (threadIdx.x < i && threadIdx.x + i < d_nPetals) {
				pOcc[threadIdx.x] += pOcc[threadIdx.x + i];
				pMpc[threadIdx.x] += pMpc[threadIdx.x + i];
				pCount[threadIdx.x] += pCount[threadIdx.x + i];
			}
			__syncthreads();
		}
						
		//Final focal cell calculations then write Occ and MPC to outData
		if (threadIdx.x == 0) {
			pOcc[0] += (focal.x * focal.z);
			pMpc[0] += (focal.x * focal.w);
			pOcc[0] *= (pCount[0] < float(d_nPetals) ? float(d_nPetals) / (pCount[0] + 1.0f) : 1.0f);
			pMpc[0] *= (pCount[0] < float(d_nPetals) ? float(d_nPetals) / (pCount[0] + 1.0f) : 1.0f);
			d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].x = pOcc[0];
			d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].y = pMpc[0];
		}
	}
}
#else
//NOTE Do NOT use causes non - reproducable errors
__global__ void OccupancyCBA_kernel(int inOffset, int outOffset, int2 *d_petalData, float2 *d_outData)
{
	//Focal cell
	float4 focal = tex2D(tex4Ref, blockIdx.x, blockIdx.y + inOffset);
	if (focal.y == 0.0f) {
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].x = 0.0f;
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].y = 0.0f;
	}
	else if (focal.x == d_noData || focal.y == d_noData) {
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].x = d_noData;
		d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].y = d_noData;
	}
	else {
		//Registers
		int i;
		int2 intPairs[4];
		float4 costNSWE{ -1.0f, -1.0f, -1.0f, -1.0f }; //also used for fetching cell values
		float4 costNDSD{ -1.0f, -1.0f, -1.0f, -1.0f };

		//Shared memory offsets
		extern __shared__ float pOcc[];         //Occupancy
		float *pMpc = &pOcc[blockDim.x];		//mpc
		float *pDist = &pOcc[2 * blockDim.x];   //Dist
		float *pCount = &pOcc[3 * blockDim.x];
		pOcc[threadIdx.x] = 0.0f;
		pMpc[threadIdx.x] = 0.0f;
		pDist[threadIdx.x] = 0.0f;
		pCount[threadIdx.x] = 0.0f;

		//Calc shared array values
		intPairs[0] = d_petalData[4 * blockDim.x + threadIdx.x];      //petal cell count, valid cell count = 0
		intPairs[0].y = 0;
		for (i = 0; i < intPairs[0].x; i++) {
			intPairs[1] = d_petalData[(i + 5) * blockDim.x + threadIdx.x];//plus 5 to jump neighbours storage
			if (intPairs[1].x != -1 && intPairs[1].y != -1) {
				costNSWE = tex2D(tex4Ref, float(int(blockIdx.x) + intPairs[1].x - d_focalOffset), float(inOffset + int(blockIdx.y) + intPairs[1].y - d_focalOffset));
				if (costNSWE.y > d_fltEps) {//test only perm grid as hab grid may be 0.0f
					pOcc[threadIdx.x] += (costNSWE.x * costNSWE.z);//Occ
					pMpc[threadIdx.x] += (costNSWE.x * costNSWE.w);//Mpc
					pDist[threadIdx.x] += costNSWE.y;
					++(intPairs[0].y);
				}
			}
		}
		if (intPairs[0].y > 0) {
			pDist[threadIdx.x] = powf(pDist[threadIdx.x] / (float)intPairs[0].y, 2.0f * sqrtf((float)intPairs[0].y * d_oneOnPi));
			pCount[threadIdx.x] = 1.0f;
		}

		//Get neighbour petals
		intPairs[0] = d_petalData[0 * blockDim.x + threadIdx.x];//N, S
		intPairs[1] = d_petalData[1 * blockDim.x + threadIdx.x];//W, E
		intPairs[2] = d_petalData[2 * blockDim.x + threadIdx.x];//NW, NE
		intPairs[3] = d_petalData[3 * blockDim.x + threadIdx.x];//SW, SE

															   //Calculate petal edge costs
		__syncthreads();
		costNSWE.x = intPairs[0].x != -1 ? sqrtf(pDist[threadIdx.x] * pDist[intPairs[0].x]) : -1.0f;
		costNSWE.y = intPairs[0].y != -1 ? sqrtf(pDist[threadIdx.x] * pDist[intPairs[0].y]) : -1.0f;
		costNSWE.z = intPairs[1].x != -1 ? sqrtf(pDist[threadIdx.x] * pDist[intPairs[1].x]) : -1.0f;
		costNSWE.w = intPairs[1].y != -1 ? sqrtf(pDist[threadIdx.x] * pDist[intPairs[1].y]) : -1.0f;
		costNDSD.x = intPairs[2].x != -1 ? powf(pDist[threadIdx.x] * pDist[intPairs[2].x], d_diag) : -1.0f;
		costNDSD.y = intPairs[2].y != -1 ? powf(pDist[threadIdx.x] * pDist[intPairs[2].y], d_diag) : -1.0f;
		costNDSD.z = intPairs[3].x != -1 ? powf(pDist[threadIdx.x] * pDist[intPairs[3].x], d_diag) : -1.0f;
		costNDSD.w = intPairs[3].y != -1 ? powf(pDist[threadIdx.x] * pDist[intPairs[3].y], d_diag) : -1.0f;

		//Calculate permeability for focal cell neighbours first
		if (threadIdx.x < 8)
			pDist[threadIdx.x] = threadIdx.x % 2 == 0 ? powf(pDist[threadIdx.x] * focal.y, d_diag) : sqrtf(pDist[threadIdx.x] * focal.y);
		else
			pDist[threadIdx.x] = 0.0f;

		//Loop until all permeabilities are maximised (i==0) then multiply by Hj 
		i = 1;
		while (__syncthreads_or(i)) {
			i = 0;
			if (intPairs[0].x != -1 && pDist[intPairs[0].x] * costNSWE.x > pDist[threadIdx.x]) {//N edge
				pDist[threadIdx.x] = pDist[intPairs[0].x] * costNSWE.x;
				i = 1;
			}
			if (intPairs[0].y != -1 && pDist[intPairs[0].y] * costNSWE.y > pDist[threadIdx.x]) {//S edge
				pDist[threadIdx.x] = pDist[intPairs[0].y] * costNSWE.y;
				i = 1;
			}
			if (intPairs[1].x != -1 && pDist[intPairs[1].x] * costNSWE.z > pDist[threadIdx.x]) {//W edge
				pDist[threadIdx.x] = pDist[intPairs[1].x] * costNSWE.z;
				i = 1;
			}
			if (intPairs[1].y != -1 && pDist[intPairs[1].y] * costNSWE.w > pDist[threadIdx.x]) {//E edge
				pDist[threadIdx.x] = pDist[intPairs[1].y] * costNSWE.w;
				i = 1;
			}
			if (intPairs[2].x != -1 && pDist[intPairs[2].x] * costNDSD.x > pDist[threadIdx.x]) {//NW edge
				pDist[threadIdx.x] = pDist[intPairs[2].x] * costNDSD.x;
				i = 1;
			}
			if (intPairs[2].y != -1 && pDist[intPairs[2].y] * costNDSD.y > pDist[threadIdx.x]) {//NE edge
				pDist[threadIdx.x] = pDist[intPairs[2].y] * costNDSD.y;
				i = 1;
			}
			if (intPairs[3].x != -1 && pDist[intPairs[3].x] * costNDSD.z > pDist[threadIdx.x]) {//SW edge
				pDist[threadIdx.x] = pDist[intPairs[3].x] * costNDSD.z;
				i = 1;
			}
			if (intPairs[3].y != -1 && pDist[intPairs[3].y] * costNDSD.w > pDist[threadIdx.x]) {//SE edge
				pDist[threadIdx.x] = pDist[intPairs[3].y] * costNDSD.w;
				i = 1;
			}
		}

		//Calc Hjwij or 0 for threads greater than nPetals
		pOcc[threadIdx.x] *= pDist[threadIdx.x];
		pMpc[threadIdx.x] *= pDist[threadIdx.x];
		__syncthreads();

		//Reduce all pDist
		for (i = d_firstReduction; i > 0; i /= 2) {
			if (threadIdx.x < i && threadIdx.x + i < d_nPetals) {
				pOcc[threadIdx.x] += pOcc[threadIdx.x + i];
				pMpc[threadIdx.x] += pMpc[threadIdx.x + i];
				pCount[threadIdx.x] += pCount[threadIdx.x + i];
			}
			__syncthreads();
		}

		//Write out Occ and Mpc to global
		if (threadIdx.x == 0) {
			pOcc[0] += (focal.x * focal.z);
			pMpc[0] += (focal.x * focal.w);
			pOcc[0] *= (pCount[0] < float(d_nPetals) ? float(d_nPetals) / (pCount[0] + 1.0f) : 1.0f);
			pMpc[0] *= (pCount[0] < float(d_nPetals) ? float(d_nPetals) / (pCount[0] + 1.0f) : 1.0f);
			d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].x = powf(pOcc[0], d_sumPower);
			d_outData[(outOffset + blockIdx.y) * gridDim.x + blockIdx.x].y = powf(pMpc[0], d_sumPower);
		}
	}
}

#endif
////////////////////////////////////////////////////////////
//REMP time series CBA small
////////////////////////////////////////////////////////////
int CUDAOccupancyCBA_TS_S(CBAParams &p)
{
	cudaError_t cudaStatus;
	msgText("Performing REMP CBA");

	//Host and Device data
	float4 *h_inData, *d_inData;
	float2 *h_outData, *d_outData;
	int2 *d_petalData;

	std::unique_ptr<float[]> h_extInData;
	std::unique_ptr<float[]> h_nhaInData;

	//Local variables
	size_t d_pitch;
	dim3 gridSize{ (unsigned int)p.nCols, 1U, 1U };
	uint maxKernelTime = h_maxKernelTime; // 100000;
	uint kernelTime = maxKernelTime;
	uint nSharedSlots = 5;
	uint i;

	//Host 7 x 4bytes x nCells, Device 6 x 4bytes x nCells
	//Malloc pinned host memory
	msgText("Allocating host memory");
	CUDA(cudaMallocHost(&h_inData, p.nCells * sizeof(float4)));
	CUDA(cudaMallocHost(&h_outData, p.nCells * sizeof(float2)));
	if (p.haveExtIn) h_extInData = std::make_unique<float[]>(p.nCells);
	if (p.haveNHAIn) h_nhaInData = std::make_unique<float[]>(p.nCells);

	//Malloc device data
	msgText("Allocating device memory");
	CUDA(cudaMalloc(&d_petalData, sizeof(int2) *p. petalRows * p.petalCols));
	CUDA(cudaMallocPitch(&d_inData, &d_pitch, p.nCols * sizeof(float4), p.nRows));
	CUDA(cudaMalloc(&d_outData, p.nCells * sizeof(float2)));

	//Copy petal data to device
	CUDA(cudaMemcpy(d_petalData, p.petalPtr.get(), sizeof(int2) * p.petalRows * p.petalCols, cudaMemcpyHostToDevice));

	//Copy constants to device
	msgText("Setting device parameters");
	CUDA(cudaMemcpyToSymbol(d_nPetals, &(p.nPetals), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_focalOffset, &(p.fOffset), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_multFocal, &(p.multFocal), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_firstReduction, &(p.firstReduction), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_sumPower, &(p.sumPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_focalPower, &(p.focalPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(p.noData), sizeof(float)));

	CUDA(cudaMemcpyToSymbol(d_useMaxPetalPerm, &(p.useMaxPetalPerm), sizeof(bool)));

	//Texture Reference
	//msgText("Texture Reference");
	//tex4Ref.filterMode = cudaFilterModePoint;
	//tex4Ref.addressMode[0] = cudaAddressModeBorder;
	//tex4Ref.addressMode[1] = cudaAddressModeBorder;
	//cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	//CUDA(cudaBindTexture2D (NULL, &tex4Ref, d_inData, &channelDesc, p.nCols, p.nRows, d_pitch));
	
	//Replace Texture Reference with Texture Object
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
	for (uint t = 0; t < p.rempTimeSteps; t++)
	{
		msgText(("Time Step: " + toStr(t)).c_str());
		//Open hab and prm grids for this t and read data
		p.habInFS.open(GridFN(p.habInFNs[t]), std::ios::in | std::ios::binary);
		p.prmInFS.open(GridFN(p.prmInFNs[t]), std::ios::in | std::ios::binary);
		if (!p.habInFS.is_open()) return msgErrorStr(-3, p.habInFNs[t]);
		if (!p.prmInFS.is_open()) return msgErrorStr(-3, p.prmInFNs[t]);
		for (i = 0; i < p.nCells; i++) {
			p.habInFS.read((char*)&(h_inData[i].x), sizeof(float));
			p.prmInFS.read((char*)&(h_inData[i].y), sizeof(float));
			h_inData[i].w = p.mpcInitFactor;//MPC does not get coupled
		}
		p.habInFS.close();
		p.prmInFS.close();

		//Seed occ factor if un-coupled or first time step
		if (!p.rempCoupled || t == 0) {
			for (i = 0; i < p.nCells; i++) {
				h_inData[i].z = p.occInitFactor;
			}
		}

		//Read occ grid if using and first time step
		if (p.haveOccIn && t == 0){
			p.occInFS.open(GridFN(p.occInFN), std::ios::in | std::ios::binary);
			if (!p.occInFS.is_open()) return msgErrorStr(-3, p.occInFN);
			for (i = 0; i < p.nCells; i++) {
				p.occInFS.read((char*)&(h_inData[i].z), sizeof(float));		
			}
			p.occInFS.close();
		}

		//Read mpc grid if using and first time step
		if (p.haveMPCIn && t == 0){
			p.mpcInFS.open(GridFN(p.mpcInFN), std::ios::in | std::ios::binary);
			if (!p.mpcInFS.is_open()) return msgErrorStr(-3, p.mpcInFN);
			for (i = 0; i < p.nCells; i++) {
				p.mpcInFS.read((char*)&(h_inData[i].w), sizeof(float));	
			}
			p.mpcInFS.close();
		}

		//Read in nha or extinction data if it's being used
		if (p.haveNHAIn) {
			p.nhaInFS.open(GridFN(p.nhaInFNs[t]), std::ios::in | std::ios::binary);
			if (!p.nhaInFS.is_open()) return msgErrorStr(-3, p.nhaInFNs[t]);
			p.nhaInFS.read((char*)(h_nhaInData.get()), p.nCells * sizeof(float));
			p.nhaInFS.close();
		}
		else if (p.haveExtIn) {
			p.extInFS.open(GridFN(p.extInFNs[t]), std::ios::in | std::ios::binary);
			if (!p.extInFS.is_open()) return msgErrorStr(-3, p.extInFNs[t]);
			p.extInFS.read((char*)(h_extInData.get()), p.nCells * sizeof(float));
			p.extInFS.close();
		}

		//Iterate CBA
		for (uint remp_i = 0; remp_i < p.rempIterations; remp_i++) {
			msgText(("CBA Iteration: " + toStr(remp_i)).c_str());
			
			//Copy data to device
			CUDA(cudaMemcpy2D(d_inData, d_pitch, h_inData, p.nCols * sizeof(float4), p.nCols * sizeof(float4), p.nRows, cudaMemcpyHostToDevice));

			//Profiled Kernel Call Loop
			msgText("Processing data on device");
			for (i = 0; i < p.nRows; i += gridSize.y) {
				gridSize.y = min(p.nRows - i, max(gridSize.y * maxKernelTime / kernelTime, 1));
				kernelTimer.Start();
				OccupancyCBA_kernel << <gridSize, p.petalCols, p.petalCols * nSharedSlots * sizeof(float) >> >(i, i, d_petalData, d_outData, tex4Obj);
				CUDA(cudaDeviceSynchronize());
				kernelTime = int(kernelTimer.Stop());
				msgProgress("Percent complete: ", i * 100 / p.nRows);
			}
			msgText("\rPercent complete: 100");

			//Copy data from device to host
			CUDA(cudaMemcpy(h_outData, d_outData, p.nCells * sizeof(float2), cudaMemcpyDeviceToHost));
			CUDA(cudaMemset(d_outData, 0, p.nCells * sizeof(float)));

			//Calculate occupancy factor grid using NHA, extinction grid or e/c parameters
			p.mpcDen = 0.0f;
			//20190322 Added S/(S+lambda/NHA) when NHA provided
			if (p.haveNHAIn) {
				for (int j = 0; j < p.nCells; j++) {
					h_inData[j].z = float((double(h_outData[j].x) * double(h_nhaInData[j])) / ((double(h_outData[j].x) * double(h_nhaInData[j])) + p.lambda));
					p.mpcDen = max(h_outData[j].y, p.mpcDen);
				}
			}
			else if (p.haveExtIn) {
				for (int j = 0; j < p.nCells; j++) {
					h_inData[j].z = h_outData[j].x > fltEps ? float(p.c * double(h_outData[j].x) / (p.c * double(h_outData[j].x) + double(h_extInData[j]))) : h_outData[j].x;
					p.mpcDen = max(h_outData[j].y, p.mpcDen);
				}
			}
			else {
				for (int j = 0; j < p.nCells; j++) {
					h_inData[j].z = h_outData[j].x > fltEps ? float(p.c * double(h_outData[j].x) / (p.c * double(h_outData[j].x) + p.e)) : h_outData[j].x;
					p.mpcDen = max(h_outData[j].y, p.mpcDen);
				}
			}

			//Calculate mpc factor grid if not last run
			if (remp_i + 1 < p.rempIterations)
				for (int j = 0; j < p.nCells; j++)
					h_inData[j].w = h_outData[j].y > fltEps ? h_outData[j].y / p.mpcDen : h_outData[j].y;

			if (remp_i + 1 < p.rempIterations && p.rempWriteEachIteration) {
				if (p.haveOccOut) {
					//Output Intermediate Si grid if writing each iteration
					p.occSiOutFS.open(GridFN(p.occOutFNs[t] + "_Si_i" + toStr(remp_i)), std::ios::out | std::ios::binary);
					if (p.occSiOutFS.is_open()) {
						for (int j = 0; j < p.nCells; j++) p.occSiOutFS.write((const char*)&(h_outData[j].x), sizeof(float));
						p.occSiOutFS.close();
						CopyFS(HeaderFN(p.habInFNs[t]), HeaderFN(p.occOutFNs[t] + "_Si_i" + toStr(remp_i)), true);
					}
					else msgErrorStr(-4, GridFN(p.occOutFNs[t] + "_Si_i" + toStr(remp_i)));
	
					//Output Intermediate Pu grid if writing each iteration
					p.occPuOutFS.open(GridFN(p.occOutFNs[t] + "_Pu_i" + toStr(remp_i)), std::ios::out | std::ios::binary);
					if (p.occPuOutFS.is_open()) {
						for (int j = 0; j < p.nCells; j++) p.occPuOutFS.write((const char*)&(h_inData[j].z), sizeof(float));
						p.occPuOutFS.close();
						CopyFS(HeaderFN(p.habInFNs[t]), HeaderFN(p.occOutFNs[t] + "_Pu_i" + toStr(remp_i)), true);
					}
					else msgErrorStr(-4, GridFN(p.occOutFNs[t] + "_Pu_i" + toStr(remp_i)));
				}
				//Output Intermediate r grid (MPC outGrid) file if writing each iteration
				if (p.haveMPCOut) {
					p.mpcOutFS.open(GridFN(p.mpcOutFNs[t] + "_r_i" + toStr(remp_i)), std::ios::out | std::ios::binary);
					if (p.mpcOutFS.is_open()) {
						for (int j = 0; j < p.nCells; j++) p.mpcOutFS.write((const char*)&(h_outData[j].y), sizeof(float));
						p.mpcOutFS.close();
						CopyFS(HeaderFN(p.habInFNs[t]), HeaderFN(p.mpcOutFNs[t] + "_r_i" + toStr(remp_i)), true);
					}
					else msgErrorStr(-4, GridFN(p.mpcOutFNs[t] + "_r_i" + toStr(remp_i)));
				}
			}
		}

		//Open and write Occ and MPC output grids for this t
		if (p.haveOccOut) {

			p.occSiOutFS.open(GridFN(p.occOutFNs[t] + "_Si"), std::ios::out | std::ios::binary);
			if (p.occSiOutFS.is_open()) {
				for (int j = 0; j < p.nCells; j++) p.occSiOutFS.write((const char*)&(h_outData[j].x), sizeof(float));
				p.occSiOutFS.close();
				CopyFS(HeaderFN(p.habInFNs[t]), HeaderFN(p.occOutFNs[t] + "_Si"), true);
			}
			else msgErrorStr(-4, GridFN(p.occOutFNs[t] + "_Si"));


			p.occPuOutFS.open(GridFN(p.occOutFNs[t]), std::ios::out | std::ios::binary);
			if (p.occPuOutFS.is_open()) {
				for (int j = 0; j < p.nCells; j++) p.occPuOutFS.write((const char*)&(h_inData[j].z), sizeof(float));
				p.occPuOutFS.close();
			}
			else msgErrorStr(-4, GridFN(p.occOutFNs[t]));
		}
		if (p.haveMPCOut) {
			p.mpcOutFS.open(GridFN(p.mpcOutFNs[t]), std::ios::out | std::ios::binary);
			if (p.mpcOutFS.is_open()) {
				for (int j = 0; j < p.nCells; j++) {
					h_outData[j].y = h_inData[j].w > fltEps ? h_outData[j].y / h_inData[j].w : h_outData[j].y;
					p.mpcOutFS.write((const char*)&(h_outData[j].y), sizeof(float));
				}
				p.mpcOutFS.close();
			}
			else msgErrorStr(-4, GridFN(p.mpcOutFNs[t]));
		}
	}//end for t in rempTimeSteps
	
	profiler.Stop();
	msgText(("Processing time: " + toStr(profiler.Total())).c_str());

	//Free pinned host and device memory
	msgText("Freeing host and device memory");
	CUDA(cudaDestroyTextureObject(tex4Obj));
	CUDA(cudaFreeHost(h_inData));
	CUDA(cudaFreeHost(h_outData));
	CUDA(cudaFree(d_petalData));
	CUDA(cudaFree(d_inData));
	CUDA(cudaFree(d_outData));
	
	cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDA_REMP_TS_CBA_S() Complete!");
	return int(cudaStatus);
}

////////////////////////////////////////////////////////////
//REMP time series CBA large
////////////////////////////////////////////////////////////
int CUDAOccupancyCBA_TS_L(CBAParams &p)
{
	cudaError_t cudaStatus;
	msgText("Performing REMP CBA");

	//Host and Device data
	float4 *h_inData, *d_inData;
	float4 *h_inBuf, *d_inBuf;
	float2 *h_outBuf, *d_outBuf;
	int2 *d_petalData;

	std::unique_ptr<float[]> h_extInBuf;
	std::unique_ptr<float[]> h_nhaInBuf;

	//Local variables
	size_t d_pitch;
	uint nBufferRows = p.fOffset * 2;
	uint nDataRows = p.fOffset * 4;
	uint nFirstRows = nDataRows - p.fOffset;
	ull nReads = 0;
	dim3 gridSize{ (unsigned int)p.nCols, 1U, 1U };
	uint maxKernelTime = h_maxKernelTime;
	uint kernelTime = maxKernelTime;
	uint nSharedSlots = 5;
	uint i, j;

	std::string mpcTmpFN;

	//Malloc pinned host memory
	msgText("Allocating host memory");
	CUDA(cudaMallocHost(&h_inData, p.nCols * nDataRows * sizeof(float4)));
	CUDA(cudaMallocHost(&h_inBuf, p.nCols * nBufferRows * sizeof(float4)));
	CUDA(cudaMallocHost(&h_outBuf, p.nCols * nBufferRows * sizeof(float2)));
	if (p.haveExtIn) h_extInBuf = std::make_unique<float[]>(nDataRows * p.nCols);
	if (p.haveNHAIn) h_nhaInBuf = std::make_unique<float[]>(nDataRows * p.nCols);

	//Malloc device data
	msgText("Allocating device memory");
	CUDA(cudaMalloc(&d_petalData, sizeof(int2) * p.petalRows * p.petalCols));
	CUDA(cudaMallocPitch(&d_inData, &d_pitch, p.nCols * sizeof(float4), nDataRows));
	CUDA(cudaMallocPitch(&d_inBuf, &d_pitch, p.nCols * sizeof(float4), nBufferRows));
	CUDA(cudaMalloc(&d_outBuf, p.nCols * nBufferRows * sizeof(float2)));

	//Copy petal data to device
	CUDA(cudaMemcpy(d_petalData, p.petalPtr.get(), sizeof(int2) * p.petalRows * p.petalCols, cudaMemcpyHostToDevice));

	//Copy constants to device
	msgText("Setting device parameters");
	CUDA(cudaMemcpyToSymbol(d_nPetals, &(p.nPetals), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_focalOffset, &(p.fOffset), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_multFocal, &(p.multFocal), sizeof(bool)));
	CUDA(cudaMemcpyToSymbol(d_firstReduction, &(p.firstReduction), sizeof(int)));
	CUDA(cudaMemcpyToSymbol(d_sumPower, &(p.sumPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_focalPower, &(p.focalPower), sizeof(float)));
	CUDA(cudaMemcpyToSymbol(d_noData, &(p.noData), sizeof(float)));

	CUDA(cudaMemcpyToSymbol(d_useMaxPetalPerm, &(p.useMaxPetalPerm), sizeof(bool)));


	//Texture Reference
	//msgText("Setting texture reference");
	//tex4Ref.filterMode = cudaFilterModePoint;
	//tex4Ref.addressMode[0] = cudaAddressModeBorder;
	//tex4Ref.addressMode[1] = cudaAddressModeBorder;
	//cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	//CUDA(cudaBindTexture2D(NULL, &tex4Ref, d_inData, &channelDesc, p.nCols, nDataRows, d_pitch));

	//Replace Texture Reference with Texture Object
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

	//Profiled Time series loop
	msgText("Processing time series");
	Profiler profiler(1000000), kernelTimer(1000000);
	profiler.Start();
	for (unsigned int t = 0; t < p.rempTimeSteps; t++) {
		msgText(("Time Step: " + toStr(t)).c_str());
		msgText("");

		//Iterate for occupancy
		for (unsigned int remp_i = 0; remp_i < p.rempIterations; remp_i++) {
			msgText(("CBA Iteration: " + toStr(remp_i)).c_str());
			msgText("");

			//Set leading rows to 0
			for (i = 0; i < p.fOffset * p.nCols; i++)
				h_inData[i].x = h_inData[i].y = h_inData[i].z = h_inData[i].w = 0.0f;

			if (p.haveNHAIn) for (i = 0; i < p.fOffset * p.nCols; i++) h_nhaInBuf[i] = 0.0f;
			if (p.haveExtIn) for (i = 0; i < p.fOffset * p.nCols; i++) h_extInBuf[i] = 0.0f;

			//Open hab and prm grids for this t and read data
			p.habInFS.open(GridFN(p.habInFNs[t]), std::ios::in | std::ios::binary);
			p.prmInFS.open(GridFN(p.prmInFNs[t]), std::ios::in | std::ios::binary);
			if (!p.habInFS.is_open()) return msgErrorStr(-3, p.habInFNs[t]);
			if (!p.prmInFS.is_open()) return msgErrorStr(-3, p.prmInFNs[t]);

			//Open nha or ext grids for this t if provided and read data
			if (p.haveNHAIn) {
				p.nhaInFS.open(GridFN(p.nhaInFNs[t]), std::ios::in | std::ios::binary);
				if (!p.nhaInFS.is_open()) return msgErrorStr(-3, p.nhaInFNs[t]);
			}
			else if (p.haveExtIn) {
				p.extInFS.open(GridFN(p.extInFNs[t]), std::ios::in | std::ios::binary);
				if (!p.extInFS.is_open()) return msgErrorStr(-3, p.extInFNs[t]);
			}

			//Occ and MPC input grids
			if (remp_i == 0 && t == 0 && p.haveOccIn)
				p.occInFS.open(GridFN(p.occInFN), std::ios::in | std::ios::binary);
			if (remp_i == 0 && t == 0 && p.haveMPCIn)
				p.mpcInFS.open(GridFN(p.mpcInFN), std::ios::in | std::ios::binary);
			
			if (remp_i == 0 && t > 0 && p.rempCoupled)
				p.occInFS.open(GridFN(p.occOutFNs[t - 1]), std::ios::in | std::ios::binary);
			
			if (remp_i > 0 && p.haveOccOut)
				p.occInFS.open(p.occPuOutFN, std::ios::in | std::ios::binary);
			if (remp_i > 0 && p.haveMPCOut)
				p.mpcInFS.open(p.mpcOutFN, std::ios::in | std::ios::binary);

			//Open Occ and MPC output grids 
			if (p.haveOccOut) {
				if (remp_i < p.rempIterations - 1) {
					p.occSiOutFN = p.occOutFNs[t] + "_Si_i" + toStr(remp_i);
					p.occPuOutFN = p.occOutFNs[t] + "_Pu_i" + toStr(remp_i);
					p.occSiOutFS.open(GridFN(p.occSiOutFN), std::ios::out | std::ios::binary);
					p.occPuOutFS.open(GridFN(p.occPuOutFN), std::ios::out | std::ios::binary);
				}
				else {
					p.occSiOutFN = p.occOutFNs[t] + "_Si";
					p.occSiOutFS.open(GridFN(p.occSiOutFN), std::ios::out | std::ios::binary);
					p.occPuOutFS.open(GridFN(p.occOutFNs[t]), std::ios::out | std::ios::binary);
				}
			}
			if (p.haveMPCOut) {
				if (remp_i < p.rempIterations - 1) {
					p.mpcOutFN = p.mpcOutFNs[t] + "_i" + toStr(remp_i);
					p.mpcOutFS.open(GridFN(p.mpcOutFN), std::ios::out | std::ios::binary);
				}
				else {
					p.mpcOutFS.open(GridFN(p.mpcOutFNs[t]), std::ios::out | std::ios::binary);
				}
				mpcTmpFN = p.mpcOutFN + "_tmp";
				p.mpcTmpFS.open(GridFN(mpcTmpFN), std::ios::out | std::ios::binary);
			}

			//Read start of hab and perm
			for (i = p.fOffset * p.nCols; i < nDataRows * p.nCols; i++) {
				p.habInFS.read((char*)&(h_inData[i].x), sizeof(float));
				p.prmInFS.read((char*)&(h_inData[i].y), sizeof(float));
				h_inData[i].z = p.occInitFactor;
				h_inData[i].w = p.mpcInitFactor;
			}

			//Read start of nha or ext grid if using
			//Updated to use S/(S+lambda/NHA)
			if (p.nhaInFS.is_open()) 
				p.nhaInFS.read((char*)&(h_nhaInBuf[p.fOffset * p.nCols]), nDataRows * p.nCols * sizeof(float));
			else if (p.extInFS.is_open()) 
				p.extInFS.read((char*)&(h_extInBuf[p.fOffset * p.nCols]), nDataRows * p.nCols * sizeof(float));

			//Read start of occ and mpc grids if open
			if (p.occInFS.is_open()) 
				for (i = p.fOffset * p.nCols; i < nDataRows * p.nCols; i++) 
					p.occInFS.read((char*)&(h_inData[i].z), sizeof(float));
			
			if (p.mpcInFS.is_open()) 
				for (i = p.fOffset * p.nCols; i < nDataRows * p.nCols; i++) 
					p.mpcInFS.read((char*)&(h_inData[i].w), sizeof(float));

			//Copy data to device
			CUDA(cudaMemcpy2D(d_inData, d_pitch, h_inData, p.nCols * sizeof(float4), p.nCols * sizeof(float4), nDataRows, cudaMemcpyHostToDevice));
			p.mpcDen = 0.0f;

			//Profiled Kernel Call Loop
			msgText("Processing data on device");
			for (i = 0; i < p.nRows; i += gridSize.y) {
				msgProgress("Percent complete: ", i * 100 / p.nRows);

				//Set kernel gridSize.y
				gridSize.y = max(min(gridSize.y * maxKernelTime / kernelTime, nBufferRows), 1);
				if (i + gridSize.y >= p.nRows) gridSize.y = p.nRows - i;

				//Parallel CBA kernel call
				kernelTimer.Start();
				OccupancyCBA_kernel << <gridSize, p.petalCols, p.petalCols * nSharedSlots * sizeof(float) >> > (p.fOffset, 0, d_petalData, d_outBuf, tex4Obj);
				CUDA(cudaDeviceSynchronize());
				kernelTime = (uint)kernelTimer.Stop();

				//Clear gridSize.y rows in the host buffer then read in any new data
				nReads = (nFirstRows + i) * p.nCols;
				for (j = 0; j < p.nCols * gridSize.y; j++) {
					h_inBuf[j].x = 0.0;
					h_inBuf[j].y = 0.0;
					h_inBuf[j].z = p.occInitFactor;
					h_inBuf[j].w = p.mpcInitFactor;
					if (p.haveExtIn) h_extInBuf[j] = 0.0f;
					if (p.haveNHAIn) h_nhaInBuf[j] = 0.0f;
					if (nReads + j < p.nCells) {
						p.habInFS.read((char*)&(h_inBuf[j].x), sizeof(float));
						p.prmInFS.read((char*)&(h_inBuf[j].y), sizeof(float));
						if (p.occInFS.is_open()) p.occInFS.read((char*)&(h_inBuf[j].z), sizeof(float));
						if (p.mpcInFS.is_open()) p.mpcInFS.read((char*)&(h_inBuf[j].w), sizeof(float));
						if (p.nhaInFS.is_open()) p.nhaInFS.read((char*)&(h_nhaInBuf[j]), sizeof(float));
						if (p.extInFS.is_open()) p.extInFS.read((char*)&(h_extInBuf[j]), sizeof(float));
					}
				}
				//Copy host buffer to device buffer then shuffle into inData
				CUDA(cudaMemcpy2D(d_inBuf, d_pitch, h_inBuf, p.nCols * sizeof(float4), p.nCols * sizeof(float4), gridSize.y, cudaMemcpyHostToDevice));
				ShuffleUp(d_inData, d_inBuf, nDataRows, gridSize.y, int(d_pitch / sizeof(float4)), true);

				//Copy and clear output device buffer
				CUDA(cudaMemcpy(h_outBuf, d_outBuf, p.nCols * gridSize.y * sizeof(float2), cudaMemcpyDeviceToHost));
				CUDA(cudaMemset(d_outBuf, 0, p.nCols * nBufferRows * sizeof(float2)));

				//Calculate occupancy factor grid and get max MPC
				//20190322 Added Si/(Si+lambda/NHA) when NHA provided
				if (p.haveNHAIn) {
					for (j = 0; j < p.nCols * gridSize.y; j++) {
						if (p.occSiOutFS.is_open()) p.occSiOutFS.write((const char*)&(h_outBuf[j].x), sizeof(float));
						//h_outBuf[j].x = h_outBuf[j].x > fltEps && h_nhaInBuf[j] > fltEps ? float(double(h_outBuf[j].x) / (double(h_outBuf[j].x) + (p.lambda / double(h_nhaInBuf[j])))) : h_outBuf[j].x;
						h_outBuf[j].x *= h_nhaInBuf[j];
						h_outBuf[j].x = float(double(h_outBuf[j].x) / (double(h_outBuf[j].x) + p.lambda));
						p.mpcDen = max(h_outBuf[j].y, p.mpcDen);
						if (p.occPuOutFS.is_open()) p.occPuOutFS.write((const char*)&(h_outBuf[j].x), sizeof(float));
						if (p.mpcTmpFS.is_open()) p.mpcTmpFS.write((const char*)&(h_outBuf[j].y), sizeof(float));
					}
				}
				else if (p.haveExtIn) {
					for (j = 0; j < p.nCols * gridSize.y; j++) {
						if (p.occSiOutFS.is_open()) p.occSiOutFS.write((const char*)&(h_outBuf[j].x), sizeof(float));
						h_outBuf[j].x = h_outBuf[j].x > fltEps ? float(p.c * double(h_outBuf[j].x) / (p.c * double(h_outBuf[j].x) + double(h_extInBuf[j]))) : h_outBuf[j].x;
						p.mpcDen = max(h_outBuf[j].y, p.mpcDen);
						if (p.occPuOutFS.is_open()) p.occPuOutFS.write((const char*)&(h_outBuf[j].x), sizeof(float));
						if (p.mpcTmpFS.is_open()) p.mpcTmpFS.write((const char*)&(h_outBuf[j].y), sizeof(float));
					}
				}
				else {
					for (j = 0; j < p.nCols * gridSize.y; j++) {
						if (p.occSiOutFS.is_open()) p.occSiOutFS.write((const char*)&(h_outBuf[j].x), sizeof(float));
						h_outBuf[j].x = h_outBuf[j].x > fltEps ? float(p.c * double(h_outBuf[j].x) / (p.c * double(h_outBuf[j].x) + p.e)) : h_outBuf[j].x;
						p.mpcDen = max(h_outBuf[j].y, p.mpcDen);
						if (p.occPuOutFS.is_open()) p.occPuOutFS.write((const char*)&(h_outBuf[j].x), sizeof(float));
						if (p.mpcTmpFS.is_open()) p.mpcTmpFS.write((const char*)&(h_outBuf[j].y), sizeof(float));
					}
				}
				CUDA(cudaDeviceSynchronize());
			}
			msgText("\rPercent complete: 100");

			//Calculate MPC grid as MPC/maxMPC or MPC/mpcIn if last iteration
			if (p.haveMPCOut) {
				p.mpcTmpFS.close();
				p.mpcTmpFS.open(GridFN(mpcTmpFN), std::ios::in | std::ios::binary);
				float mpcNum, mpcVal;
				p.mpcInFS.seekg(0);
				for (j = 0; j < p.nCells; j++) {
					p.mpcTmpFS.read((char *)&mpcNum, sizeof(float));
					if (remp_i == p.rempIterations - 1)
						p.mpcInFS.read((char *)&(p.mpcDen), sizeof(float));
					mpcVal = p.mpcDen > fltEps ? mpcNum / p.mpcDen : mpcNum;
					p.mpcOutFS.write((const char*)&mpcVal, sizeof(float));
				}
			}

			//Close all open files and write output headers
			p.habInFS.close();
			p.prmInFS.close();
			if (p.nhaInFS.is_open()) p.nhaInFS.close();
			if (p.extInFS.is_open()) p.extInFS.close();
			if (p.occInFS.is_open()) p.occInFS.close();
			if (p.mpcInFS.is_open()) p.mpcInFS.close();
			if (p.occPuOutFS.is_open()) p.occPuOutFS.close();
			if (p.occSiOutFS.is_open()) p.occSiOutFS.close();
			if (p.mpcOutFS.is_open()) p.mpcOutFS.close();
			if (p.mpcTmpFS.is_open()) p.mpcTmpFS.close();
		}// end for remp_i in rempIterations
	}//end for t in rempTimeSteps

	profiler.Stop();
	msgText(("Processing time: " + toStr(profiler.Total())).c_str());

	//Free pinned host and device memory and ext buffer
	msgText("Freeing host and device memory");
	CUDA(cudaDestroyTextureObject(tex4Obj));
	CUDA(cudaFreeHost(h_inData));
	CUDA(cudaFreeHost(h_inBuf));
	CUDA(cudaFreeHost(h_outBuf));
	CUDA(cudaFree(d_inData));
	CUDA(cudaFree(d_inBuf));
	CUDA(cudaFree(d_outBuf));
	CUDA(cudaFree(d_petalData));

	cudaStatus = cudaGetLastError();
	CUDA(cudaDeviceReset());
	msgText((std::string("Device status ") + cudaGetErrorString(cudaStatus)).c_str());
	msgText("CUDA_REMP_TS_CBA_L() Complete!");
	return int(cudaStatus);
}

/**/
#endif