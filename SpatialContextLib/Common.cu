/*
Common.cu - Common CUDA device variables and functions
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


#ifndef _CUDA_CBA_COMMON_CU_
#define _CUDA_CBA_COMMON_CU_
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <iostream>
#include "Common.cuh"
#include "Parameters.h"
#include "FileUtilsLib.h"

//CUDA error checking wrapper
inline void CUDACheck(cudaError_t code, const char *file, int line, bool abort) {
	if (code != cudaSuccess) {
		char msg[256];
		sprintf(msg, "CUDA: %s  %s  %d\n", cudaGetErrorString(code), file, line);
		msgText(msg);
		if (abort) {
			cudaDeviceReset();
			throw std::runtime_error(msg);
		}
	}
}

unsigned int h_maxKernelTime = 10000;

//Define constant device variables here to avoid multiple def errors in VS17
//#if defined _WIN32 && _MSC_VER < 1920
#if defined _WIN32
//Common Host Variables
const float fltEps = FLT_MIN;//1.175494351e-38F

 //Common Defined Device Constants
 __constant__ __device__ float	d_oneOnPi =	M_1_PI;//0.318309886183790671538;
 __constant__ __device__ float	d_diag = M_SQRT1_2;//0.707106781186547524401;
 __constant__ __device__ float	d_fltEps = FLT_MIN;//1.175494351e-38F
 __constant__ __device__ float	d_sqrt2 = M_SQRT2;//1.4142135623730950488016887242097

 //Common Derived Device Constants with default values
 __constant__ __device__ bool	d_multFocal = true;
 __constant__ __device__ bool	d_quickEHA = false;
 __constant__ __device__ int	d_nPetals = 0;
 __constant__ __device__ int	d_focalOffset = 0;
 __constant__ __device__ int	d_firstReduction = 1;//int(powf(2, ceil(log2(float(p.petalCols))))) / 2
 __constant__ __device__ uint	d_nCols = 1;
 __constant__ __device__ uint	d_nColsAligned = 32;
 __constant__ __device__ uint	d_nRows = 1;
 __constant__ __device__ uint	d_nCells = 1;
 __constant__ __device__ float	d_sumPower = 1.0f;
 __constant__ __device__ float	d_focalPower = 1.0f;
 __constant__ __device__ float	d_noData = -9999.0f;
 __constant__ __device__ float	d_maxCond = 1.0f;
 __constant__ __device__ float	d_maxPerm = 1.0f;

 //Occupancy
 //__constant__ __device__ double	d_occC = 1.0f;
 //__constant__ __device__ double	d_occE = 1.0f;
 __constant__ __device__ bool	d_useMaxPetalPerm = false;

 //Resilience
 __constant__ __device__ uint d_nTxLyrs = 0;
 __constant__ __device__ uint d_txDataOffset = 0;
 __constant__ __device__ uint d_nDispersals = 0;

 //Dispersal probability
 __constant__ __device__ bool d_doAveragePetalHabitat = false;
 __constant__ __device__ bool d_doDispersalProbability = false;
 __constant__ __device__ float d_dispersalYears = 50;

 //Gaussian dispersal
 __constant__ __device__ uint d_nGenerations = 1;
 __constant__ __device__ float d_habFactor = 1.0f;
 __constant__ __device__ float d_gaussConst = 1.0f;
 __constant__ __device__ float d_cellSize = 1.0f;

 //GDM Sim to class
 __constant__ __device__ uint d_nClasses = 0;


 //Texture references
 //texture<float, cudaTextureType2D, cudaReadModeElementType> tex1Ref;
 //texture<float2, cudaTextureType2D, cudaReadModeElementType> tex2Ref;
 //texture<float4, cudaTextureType2D, cudaReadModeElementType> tex4Ref;
#endif

static 	cudaDeviceProp deviceProp; 

//Common Device Functions
bool CUDASetDevice(int d)
{
	try {
		int devCount;
		CUDA(cudaGetDeviceCount(&devCount));
		if (devCount - 1 < d) return false;
		CUDA(cudaSetDevice(d));
		return true;
	}
	catch (std::exception) {
		return false;
	}
}

bool CUDAGetGlobalMem(size_t &globMem, int d)
{
	try {
		int devCount;
		CUDA(cudaGetDeviceCount(&devCount));
		if (devCount - 1 < d) return false;
		CUDA(cudaGetDeviceProperties(&deviceProp, d));
		globMem = deviceProp.totalGlobalMem;
		return true;
	}
	catch (std::exception) {
		return false;
	}
}

bool CUDAGetMemInfo(size_t &devFree, size_t &devTotal, int d)
{
	try {
		int devCount;
		CUDA(cudaGetDeviceCount(&devCount));
		if (devCount - 1 < d) return false;
		CUDA(cudaMemGetInfo(&devFree, &devTotal));
		return true;
	}
	catch (std::exception) {
		return false;
	}
}

//Need to buffer when nBands x nCells x 4bytes > device free mem or 2GB if x32
//Throws Unable to get device memory exception
bool BufferRequired(ull nBands, ull nCells, int d) {
	size_t freeMem, totalMem, spareMem = size_t(1 << 28);//268,435,456 bytes left spare
	//if (!CUDAGetMemInfo(freeMem, totalMem, d))throw std::exception("Unable to get device memory");// return msgErrorStr(-10);
	if (!CUDAGetMemInfo(freeMem, totalMem, d))throw std::exception();// no std::exception(char *) or (std::string) constructor
	return
		(nBands * nCells * sizeof(float) > freeMem - spareMem) ||
		(int(sizeof(size_t)) == 4 && nBands * nCells * sizeof(float) > INT_MAX);
}

//Shuffle up kernel
template<typename T>
__global__ void ShuffleUp_kernel(T *dest, T *buf, int destRows, int bufRows, int width)
{
	uint r, c = blockIdx.x * blockDim.x + threadIdx.x;
	if (c < width) {
		for (r = 0; r < destRows - bufRows; r++)
			dest[r * width + c] = dest[(r + bufRows) * width + c];
		for (r = 0; r < bufRows; r++)
			dest[(destRows - bufRows + r) * width + c] = buf[r * width + c];
	}
}

//Shuffle up kernel wrapper
template<typename T>
bool ShuffleUp(T *dest, T *buf, int destRows, int bufRows, int width, bool useKernal)
{
	uint blockDimx = 512;
	uint gridDimx = (width + blockDimx - 1) / blockDimx;
	ShuffleUp_kernel<<<gridDimx, blockDimx>>>(dest, buf, destRows, bufRows, width);
	return true;
}

//Shuffle tx data up kernel for resilience
template<typename T>
__global__ void ShuffleTxUp_kernel(T *dest, T *buf, int destRows, int bufRows, int width)
{
	uint t, r, c = blockIdx.x * blockDim.x + threadIdx.x;
	if (c < width) {
		for (t = 0; t < d_nTxLyrs; t++) {
			uint txDataOffset = t * destRows * width;
			uint txBufOffset = t * bufRows * width;
			for (r = 0; r < destRows - bufRows; r++)
				dest[txDataOffset + r * width + c] = dest[txDataOffset + (r + bufRows) * width + c];
			for (r = 0; r < bufRows; r++)
				dest[txDataOffset + (destRows - bufRows + r) * width + c] = buf[txBufOffset + r * width + c];
		}
	}
}

//Shuffle tx data up kernel wrapper
template<typename T>
bool ShuffleTxUp(T *dest, T *buf, int destRows, int bufRows, int width, bool useKernal)
{
	uint blockDimx = 512;
	uint gridDimx = (width + blockDimx - 1) / blockDimx;
	ShuffleTxUp_kernel << <gridDimx, blockDimx >> >(dest, buf, destRows, bufRows, width);
	return true;
}

//Shuffle Up, extern templates need expected type declarations 
//template<typename T> bool ShuffleUp(T *dest, T *buf, int destRows, int bufRows, int width, bool useKernal = true);
template bool ShuffleUp(float *, float *, int, int, int, bool);
template bool ShuffleUp(float2 *, float2 *, int, int, int, bool);
template bool ShuffleUp(float4 *, float4 *, int, int, int, bool);

//Shuffle tx data up for resilience
//template<typename T> bool ShuffleTxUp(T *dest, T *buf, int destRows, int bufRows, int width, bool useKernal = true);
template bool ShuffleTxUp(float *, float *, int, int, int, bool);



#endif