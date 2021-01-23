#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

using namespace std;

typedef unsigned int uint32_t;
typedef unsigned short int uint16_t;
typedef float float32_t;
typedef double float64_t;

// ====================================================
// Util functions to check CUDA GPU compatibility
// ====================================================
bool cudaCheck(cudaError_t cudaStatus);
int checkCompatibility(int gpu_id);
void cudaCheckLastError();

// ====================================================
// Function to preprocess data for GPU
// ====================================================
void preprocessDataForGPU(uint32_t* data, int NUM_COMPARTMENTS, uint32_t* compartmentsPerBlock, uint32_t* offsetPerBlock, int NUM_BLOCKS);

// ====================================================
// CUDA Kernels for Ax operation
// ====================================================
__global__ void multiply_Ax_ICpart(
    uint32_t*  voxelIDs,
    uint32_t*  fiberIDs,
    uint16_t*  orienIDs,
    float32_t* lengths,
    uint32_t*  segmentsPerBlock,
    uint32_t*  offsetPerBlock,
    float32_t* lut,
    float64_t* x,
    float64_t* y);

__global__ void multiply_Ax_ECpart(
        uint32_t*  voxelIDs,
        uint16_t*  orienIDs,
        uint32_t*  segmentsPerBlock,
        uint32_t*  offsetPerBlock,
        float32_t* lut,
        float64_t* x,
        float64_t* y);

__global__ void multiply_Ax_ISOpart(
    float32_t* lut,
    float64_t* x,
    float64_t* y);

// ====================================================
// CUDA Kernels for A'y operation
// ====================================================
__global__ void multiply_Aty_ICpart(
    uint32_t*  TvoxelIC,
    uint32_t*  TfiberIC,
    uint16_t*  TorienIC,
    float32_t* TlengthIC,
    uint32_t*  compartmentsPerBlock,
    uint32_t*  offsetPerBlock,
    float32_t* lut,
    float64_t* x,
    float64_t* y);

__global__ void multiply_Aty_ECpart(
    uint32_t*  voxelEC,
    uint16_t*  orienEC,
    uint32_t*  segmentsPerBlock,
    uint32_t*  offsetPerBlock,
    float32_t* lut,
    float64_t* x,
    float64_t* y);

__global__ void multiply_Aty_ISOpart(
    float* lut,
    double* x,
    double* y);

// ====================================================
// Constant global values in the GPU
// ====================================================
__constant__ int NUM_VOXELS;
__constant__ int NUM_FIBERS;
__constant__ int NUM_PEAKS;
__constant__ int NUM_ORIENTATIONS;
__constant__ int NUM_SAMPLES;
__constant__ int NUM_DIAMETERS;
__constant__ int NUM_ZEPPELINS;
__constant__ int NUM_BALLS;
__constant__ int NUM_ROWS;        
__constant__ int NUM_COLS;      
__constant__ int SIZE_LUTIC;      
__constant__ int SIZE_LUTEC;     
__constant__ int SIZE_LUTISO;

// ====================================================
// Pointers to A (IC part) in the GPU
// ====================================================
static uint32_t*  gpu_voxelIC;
static uint32_t*  gpu_fiberIC;
static uint16_t*  gpu_orienIC;
static float32_t* gpu_lengthIC;
static uint32_t*  gpu_segmentsPerBlockIC;
static uint32_t*  gpu_offsetPerBlockIC;

// ====================================================
// Pointers to A' (IC part) in the GPU
// ====================================================
static uint32_t*  gpu_TvoxelIC;
static uint32_t*  gpu_TfiberIC;
static uint16_t*  gpu_TorienIC;
static float32_t* gpu_TlengthIC;
static uint32_t*  gpu_TfibersPerBlockIC;
static uint32_t*  gpu_ToffsetPerBlockIC;

// ====================================================
// Pointers to A (EC part) in the GPU
// ====================================================
static uint32_t* gpu_voxelEC;
static uint16_t* gpu_orienEC;
static uint32_t* gpu_segmentsPerBlockEC;
static uint32_t* gpu_offsetPerBlockEC;

// ====================================================
// Pointers to LUTs in the GPU
// ====================================================
static float32_t* gpu_lutIC;
static float32_t* gpu_lutEC;
static float32_t* gpu_lutISO;

// ====================================================
// Pointers to x and y in the GPU
// ====================================================
static float64_t* gpu_x;
static float64_t* gpu_y;

// ============================================================================
// This class creates an instance of the LinearOperator in GPU memory
// ============================================================================
class CudaLinearOperator {

    // constant values in CPU
    int nrows;
    int ncols;
    int nvoxels;
    int nfibers;
    int nsegments;

    // CUDA GPU status
    bool cudaStatus;
    int  cudaError;

    public:
        CudaLinearOperator(
            // pointers to IC data in CPU memory
            uint32_t* voxelIC,
            uint32_t* fiberIC,
            uint16_t* orienIC,
            float*    lengthIC,
            float*    lutIC,
            // pointers to EC data in CPU memory
            uint32_t* voxelEC,
            uint16_t* orienEC,
            float*    lutEC,
            // pointer to ISO data in CPU memory
            float*    lutISO,
            // dataset constant values
            int nsegments,
            int nvoxels,      
            int nfibers,      
            int npeaks,
            int norientations,
            int nsamples,     
            int ndiameters,   
            int nzeppelins,   
            int nballs,
            // flag to ensure we create the operator only one time
            int fcall,
            // id of the selected CUDA gpu
            int gpu_id);

        ~CudaLinearOperator();

        int  getCudaStatus() { return (int)cudaStatus; }
        void setTransposeData(uint32_t*  voxelIDs, uint32_t*  fiberIDs, uint16_t*  orienIDs, float32_t* lengths);
        void destroy();

        void  dot(float64_t* v_in, float64_t* v_out);
        void Tdot(float64_t* v_in, float64_t* v_out);
};