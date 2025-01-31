/*
Common.cuh - Common CUDA device variables and functions
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


#ifndef _CUDA_CBA_COMMON_H_
#define _CUDA_CBA_COMMON_H_
#define _USE_MATH_DEFINES

#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <climits>
#include <cfloat>

//typedef unsigned int uint;
using uint = unsigned int;
using ull = unsigned long long;
//CUDA error checking macro
#define CUDA(call) { CUDACheck((call), __FILE__, __LINE__); }
inline void CUDACheck(cudaError_t code, const char *file, int line, bool abort = true);

//PreCalcEdge kernels disabled until resolved
//#define USE_PRECALC_EDGE_KERNELS 1

//Define constant device variables in Common.cu to avoid multiple def errors in VS17
//#if defined _WIN32 && _MSC_VER < 1920
#if defined _WIN32
//Common Host Variables
extern unsigned int h_maxKernelTime;
extern const float fltEps;// = FLT_MIN;//1.175494351e-38F

//Common Defined Device Constants
 extern __constant__ __device__ float	d_oneOnPi;
 extern __constant__ __device__ float	d_diag;
 extern __constant__ __device__ float	d_fltEps;
 extern __constant__ __device__ float	d_sqrt2;

 //Common runtime derived Device Constants
 extern __constant__ __device__ bool	d_multFocal;
 extern __constant__ __device__ bool	d_quickEHA;
 extern __constant__ __device__ int		d_nPetals;
 extern __constant__ __device__ int		d_focalOffset;
 extern __constant__ __device__ int		d_firstReduction;
 extern __constant__ __device__ uint	d_nCols;
 extern __constant__ __device__ uint	d_nColsAligned;
 extern __constant__ __device__ uint	d_nRows;
 extern __constant__ __device__ uint	d_nCells;
 extern __constant__ __device__ float	d_sumPower;
 extern __constant__ __device__ float	d_focalPower;
 extern __constant__ __device__ float	d_noData;
 extern __constant__ __device__ float	d_maxCond;
 extern __constant__ __device__ float	d_maxPerm;

 //Occupancy
 //extern __constant__ __device__ double	d_occC;
 //extern __constant__ __device__ double	d_occE;
 extern __constant__ __device__ bool		d_useMaxPetalPerm;

 //Resilience
 extern __constant__ __device__ uint d_nTxLyrs;
 extern __constant__ __device__ uint d_txDataOffset;
 extern __constant__ __device__ uint d_nDispersals;

 //Testing dispersal probability
 extern __constant__ __device__ bool d_doAveragePetalHabitat;
 extern __constant__ __device__ bool d_doDispersalProbability;
 extern __constant__ __device__ float d_dispersalYears;

 //Testing gaussian dispersal
 extern __constant__ __device__ uint d_nGenerations;
 extern __constant__ __device__ float d_habFactor;
 extern __constant__ __device__ float d_gaussConst;
 extern __constant__ __device__ float d_cellSize;

 //GDM Sim to class
 extern __constant__ __device__ uint d_nClasses;

 //Texture references
 //extern texture<float, cudaTextureType2D, cudaReadModeElementType> tex1Ref;
 //extern texture<float2, cudaTextureType2D, cudaReadModeElementType> tex2Ref;
 //extern texture<float4, cudaTextureType2D, cudaReadModeElementType> tex4Ref;

#else
//Defining constant device variables here seems fine for GCC and _MSC_VER > 1920

//Common Host Variables
extern unsigned int h_maxKernelTime;// = 10000;
const float fltEps = FLT_MIN;//1.175494351e-38F

//Common Defined Device Constants
__constant__ __device__ float	d_oneOnPi =	M_1_PI;//0.318309886183790671538;
__constant__ __device__ float	d_diag = M_SQRT1_2;//0.707106781186547524401;
__constant__ __device__ float	d_fltEps = FLT_MIN;//1.175494351e-38F
__constant__ __device__ float	d_sqrt2 = M_SQRT2;//1.4142135623730950488016887242097

//Common Derived Device Constants with default values
__constant__ __device__ bool	d_multFocal = true;
__constant__ __device__ bool	d_quickEHA = false;
__constant__ __device__ int		d_nPetals = 0;
__constant__ __device__ int		d_focalOffset = 0;
__constant__ __device__ int		d_firstReduction = 1;//int(powf(2, ceil(log2(float(p.petalCols))))) / 2
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



//Layered textures and surfaces for tx grids
//texture<float, cudaTextureType2DLayered> txSrcLyrs, txDstLyrs;
//surface<float, cudaSurfaceType2DLayered> txSrcSurf, txDstSurf;


//Common Device Functions
bool CUDASetDevice(int d = 0);
bool CUDAGetGlobalMem(size_t &globMem, int d = 0);
bool CUDAGetMemInfo(size_t &devFree, size_t &devTotal, int d = 0);
bool BufferRequired(ull nBands, ull nCells, int d = 0);

// //Shuffle Up, extern templates need expected type declarations 
template<typename T> bool ShuffleUp(T *dest, T *buf, int destRows, int bufRows, int width, bool useKernal = true);
// template bool ShuffleUp(float *, float *, int, int, int, bool);
// template bool ShuffleUp(float2 *, float2 *, int, int, int, bool);
// template bool ShuffleUp(float4 *, float4 *, int, int, int, bool);

// //Shuffle tx data up for resilience
template<typename T> bool ShuffleTxUp(T *dest, T *buf, int destRows, int bufRows, int width, bool useKernal = true);
// template bool ShuffleTxUp(float *, float *, int, int, int, bool);

#endif