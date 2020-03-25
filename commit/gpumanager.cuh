#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

//#define __dual__ __host__ __device__

using namespace std;

typedef unsigned int uint32_t;
typedef unsigned short int uint16_t;
typedef float float32_t;
typedef double float64_t;

bool cudaCheck(cudaError_t cudaStatus);
void preprocessDataForGPU(uint32_t* data, int NUM_COMPARTMENTS, uint32_t* compartmentsPerBlock, uint32_t* offsetPerBlock, int NUM_BLOCKS);

/*class segment_t {
    public:
    
    // pointer to the GPU memory where the array is stored
    uint32_t voxelID;
    uint32_t fiberID;
    uint16_t orienID;
    float length;
  
    __dual__  segment();
    __dual__ ~segment();
};//*/

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

class CudaLinearOperator {

    // pointers to IC data in GPU memory
    uint32_t*  voxelIC;
    uint32_t*  fiberIC;
    uint16_t*  orienIC;
    float32_t* lengthIC;

    // pointers to IC data (transpose) in GPU memory
    uint32_t*  voxelICt;
    uint32_t*  fiberICt;
    uint16_t*  orienICt;
    float32_t* lengthICt;
    uint32_t* fibersPerBlockICt;
    uint32_t* offsetPerBlockICt;

    // auxiliar arrays for GPU
    uint32_t* segmentsPerBlockIC;
    uint32_t* offsetPerBlockIC;
    uint32_t* segmentsPerBlockEC;
    uint32_t* offsetPerBlockEC;

    // pointers to EC data in GPU memory
    uint32_t*  voxelEC;
    uint16_t*  orienEC;

    // pointers to LUTs in GPU memory
    float32_t* lutIC;
    float32_t* lutEC;
    float32_t* lutISO;

    // pointers to vector x and y
    float64_t* x;
    float64_t* y;

    // dimensions of the operator
    int nrows;
    int ncols;
    int nvoxels;
    int nfibers;

    public:
        CudaLinearOperator(
            uint32_t* voxelIC,
            uint32_t* fiberIC,
            uint16_t* orienIC,
            float32_t*    lengthIC,
            float32_t*    lutIC,
        
            uint32_t* voxelEC,
            uint16_t* orienEC,
            float32_t*    lutEC,
        
            float32_t*    lutISO,
        
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

        void setTransposeData(uint32_t* voxelIDs, uint32_t* fiberIDs, uint16_t* orienIDs, float32_t* lengths);
        void multiplyByX(float64_t* x, float64_t* y);
        void multiplyByY(float64_t* y, float64_t* x);
};