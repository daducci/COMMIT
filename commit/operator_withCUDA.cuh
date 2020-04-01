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
bool checkCompatibility(size_t required_mem, int gpu_id = 0) {
    int num_gpus;
    cudaError_t cudaStatus;
    
    //printf("-> Checking availability of CUDA:\n");
    cudaStatus = cudaGetDeviceCount(&num_gpus);

    if (num_gpus <= 0 || num_gpus <= gpu_id) {
        printf("\t* the selected GPU does not exist or is not detected \n");
        return false;
    }

    if(cudaStatus == cudaSuccess){
        cudaDeviceProp gpu_properties;
        cudaGetDeviceProperties(&gpu_properties, gpu_id);

        printf("\t* checking availability of CUDA ... [ OK ]\n");
        printf("\t* number of CUDA GPUs detected: %d\n", num_gpus);
        printf("\t* using GPU %s with ID %d... \n", gpu_properties.name, gpu_id);

        if (required_mem <= gpu_properties.totalGlobalMem) {
            printf("\t* using %f GB of total %f GB... [ OK ]\n", required_mem*1e-9, gpu_properties.totalGlobalMem*1e-9);
        }
        else {
            printf("\t* using %f GB of total %f GB... [ ERROR ]: dictionary too big for GPU memory\n", required_mem*1e-9, gpu_properties.totalGlobalMem*1e-9);
        }

        if(gpu_properties.major >= 5){
            printf("\t* compute capability: %d.%d [ OK ]\n", gpu_properties.major, gpu_properties.minor);
        }
        else{
            printf("\t* compute capability: %d.%d [ ERROR ]. GPU compute capability must be at least 5.0\n", gpu_properties.major, gpu_properties.minor);
            return false;
        }

        return true;
    }
    else{
        printf("\t* checking availability of CUDA ... [ ERROR ]: CUDA is not available or GPU is not CUDA compatible\n");
        return false;
    }
}

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

// pointers to IC data in GPU memory
static uint32_t*  gpu_voxelIC;
static uint32_t*  gpu_fiberIC;
static uint16_t*  gpu_orienIC;
static float32_t* gpu_lengthIC;
static uint32_t*  gpu_segmentsPerBlockIC;
static uint32_t*  gpu_offsetPerBlockIC;

// pointers to IC data (transpose) in GPU memory
static uint32_t*  gpu_TvoxelIC;
static uint32_t*  gpu_TfiberIC;
static uint16_t*  gpu_TorienIC;
static float32_t* gpu_TlengthIC;
static uint32_t*  gpu_TfibersPerBlockIC;
static uint32_t*  gpu_ToffsetPerBlockIC;

// pointers to EC data in GPU memory
static uint32_t* gpu_voxelEC;
static uint16_t* gpu_orienEC;
static uint32_t* gpu_segmentsPerBlockEC;
static uint32_t* gpu_offsetPerBlockEC;

// pointers to LUTs in GPU memory
static float32_t* gpu_lutIC;
static float32_t* gpu_lutEC;
static float32_t* gpu_lutISO;

// pointers to vector x and y
static float64_t* gpu_x;
static float64_t* gpu_y;

class CudaLinearOperator {

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
            int nballs,
        
            bool fcall);

        ~CudaLinearOperator();

        int  getCudaStatus() { return (int)cudaStatus; }
        void setTransposeData(uint32_t*  voxelIDs, uint32_t*  fiberIDs, uint16_t*  orienIDs, float32_t* lengths);
        void destroy();

        void  dot(float64_t* v_in, float64_t* v_out);
        void Tdot(float64_t* v_in, float64_t* v_out);
};