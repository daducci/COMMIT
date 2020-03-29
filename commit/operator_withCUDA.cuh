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

bool cudaCheck(cudaError_t cudaStatus);
void preprocessDataForGPU(uint32_t* data, int NUM_COMPARTMENTS, uint32_t* compartmentsPerBlock, uint32_t* offsetPerBlock, int NUM_BLOCKS);

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

class CudaLinearOperator {

    // pointers to IC data in GPU memory
    uint32_t*  voxelIC;
    uint32_t*  fiberIC;
    uint16_t*  orienIC;
    float32_t* lengthIC;
    uint32_t*  segmentsPerBlockIC;
    uint32_t*  offsetPerBlockIC;

    // pointers to IC data (transpose) in GPU memory
    uint32_t*  TvoxelIC;
    uint32_t*  TfiberIC;
    uint16_t*  TorienIC;
    float32_t* TlengthIC;
    uint32_t*  TfibersPerBlockIC;
    uint32_t*  ToffsetPerBlockIC;

    // pointers to EC data in GPU memory
    uint32_t* voxelEC;
    uint16_t* orienEC;
    uint32_t* segmentsPerBlockEC;
    uint32_t* offsetPerBlockEC;

    // pointers to LUTs in GPU memory
    float32_t* lutIC;
    float32_t* lutEC;
    float32_t* lutISO;

    // pointers to vector x and y
    float64_t* x;
    float64_t* y;

    // constant values in GPU
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

    // constant values in CPU
    int nrows;
    int ncols;
    int nvoxels;
    int nfibers;
    int nsegments;

    // CUDA GPU status
    bool cudaStatus;

    public:
        CudaLinearOperator(
            uint32_t* voxelIC,
            uint32_t* fiberIC,
            uint16_t* orienIC,
            float*    lengthIC,
            float*    lutIC,
        
            uint32_t* voxelEC,
            uint16_t* orienEC,
            float*    lutEC,
        
            float*    lutISO,
        
            int nsegments,
            int nvoxels,      
            int nfibers,      
            int npeaks,       
            int norientations,
            int nsamples,     
            int ndiameters,   
            int nzeppelins,   
            int nballs);

        ~CudaLinearOperator();

        bool getCudaStatus() const { return cudaStatus; }
        void setTransposeData(uint32_t*  voxelIDs, uint32_t*  fiberIDs, uint16_t*  orienIDs, float32_t* lengths);

        void  dot(float64_t* v_in, float64_t* v_out);
        void Tdot(float64_t* v_in, float64_t* v_out);
};