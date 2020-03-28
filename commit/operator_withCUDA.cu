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
            int nballs)
        {
            this->nvoxels = nvoxels;
            this->nfibers = nfibers;
            this->nrows = nvoxels * nsamples;
            this->ncols = nfibers*ndiameters + npeaks*nzeppelins + nvoxels*nballs;
            int size_lutic  = ndiameters*norientations*nsamples;
            int size_lutec  = nzeppelins*norientations*nsamples;
            int size_lutiso = nballs*nsamples;
            bool status;
        
            uint32_t* segmentsPerBlock = (uint32_t*) malloc(nvoxels*sizeof(uint32_t));
            uint32_t* offsetPerBlock   = (uint32_t*) malloc(nvoxels*sizeof(uint32_t));
        
            // copy constant values to GPU
            printf("\t* constant global values ... ");
            status = true;
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_VOXELS,       &nvoxels,       sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_FIBERS,       &nfibers,       sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_PEAKS,        &npeaks,        sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_ORIENTATIONS, &norientations, sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_SAMPLES,      &nsamples,      sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_DIAMETERS,    &ndiameters,    sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_ZEPPELINS,    &nzeppelins,    sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_BALLS,        &nballs,        sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_ROWS,         &nrows,         sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(NUM_COLS,         &ncols,         sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTIC,       &size_lutic,    sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTEC,       &size_lutec,    sizeof(int)) );
            status = status && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTISO,      &size_lutiso,   sizeof(int)) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
        
            // alloc memory in GPU for vectors x and y
            printf("\t* memory for vectors x and y ... ");
            status = true;
            status = status && cudaCheck( cudaMalloc((void**)&(this->x), ncols*sizeof(float64_t)) );
            status = status && cudaCheck( cudaMalloc((void**)&(this->y), nrows*sizeof(float64_t)) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            // alloc GPU memory for segments
            printf("\t* memory for LUT (IC part) ... ");
            status = true;
            status = status && cudaCheck( cudaMalloc((void**)&(this->lutIC), size_lutic*sizeof(float32_t)) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* copying LUT in GPU (IC part) ... ");
            status = true;
            status = status && cudaCheck( cudaMemcpy(this->lutIC, lutIC, size_lutic*sizeof(float32_t), cudaMemcpyHostToDevice) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* allocating memory for LUT in GPU (EC part) ... ");
            status = cudaCheck( cudaMalloc((void**)&(this->lutEC), size_lutec*sizeof(float32_t)) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* copying LUT in GPU (EC part) ... ");
            status = cudaCheck( cudaMemcpy(this->lutEC, lutEC, size_lutec*sizeof(float32_t), cudaMemcpyHostToDevice) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* allocating memory for LUT in GPU (ISO part) ... ");
            status = cudaCheck( cudaMalloc((void**)&(this->lutISO), size_lutiso*sizeof(float32_t)) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* copying LUT in GPU (ISO part) ... ");
            status = cudaCheck( cudaMemcpy(this->lutISO, lutISO, size_lutiso*sizeof(float32_t), cudaMemcpyHostToDevice) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* preprocessing data for GPU ... ");
            preprocessDataForGPU(voxelIC, nsegments, segmentsPerBlock, offsetPerBlock, nvoxels);
            printf("\n");
        
            printf("\t* fiber segments memory allocation ... ");
            status = true;
            status = status && cudaCheck( cudaMalloc((void**)&(this->voxelIC),  nsegments*sizeof(uint32_t))  );
            status = status && cudaCheck( cudaMalloc((void**)&(this->fiberIC),  nsegments*sizeof(uint32_t))  );
            status = status && cudaCheck( cudaMalloc((void**)&(this->orienIC),  nsegments*sizeof(uint16_t))  );
            status = status && cudaCheck( cudaMalloc((void**)&(this->lengthIC), nsegments*sizeof(float32_t)) );
            status = status && cudaCheck( cudaMalloc((void**)&(this->segmentsPerBlockIC), nvoxels*sizeof(uint32_t)) );
            status = status && cudaCheck( cudaMalloc((void**)&(this->offsetPerBlockIC),   nvoxels*sizeof(uint32_t)) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* transfering fiber segments ... ");
            status = true;
            status = status && cudaCheck( cudaMemcpy(this->voxelIC,  voxelIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(this->fiberIC,  fiberIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(this->orienIC,  orienIC,  nsegments*sizeof(uint16_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(this->lengthIC, lengthIC, nsegments*sizeof(float32_t), cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(this->segmentsPerBlockIC, segmentsPerBlock, nvoxels*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(this->offsetPerBlockIC,   offsetPerBlock,   nvoxels*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            // ---------------------------------------- EC DATA ---------------------------------------- //
            printf("\t* allocating memory for operator A in GPU (EC part) ... ");
            status = true;
            status = status && cudaCheck( cudaMalloc((void**)&(this->voxelEC),  npeaks*sizeof(uint32_t)) );
            status = status && cudaCheck( cudaMalloc((void**)&(this->orienEC),  npeaks*sizeof(uint16_t)) );
            status = status && cudaCheck( cudaMalloc((void**)&(this->segmentsPerBlockEC), nvoxels*sizeof(uint32_t))  );
            status = status && cudaCheck( cudaMalloc((void**)&(this->offsetPerBlockEC),   nvoxels*sizeof(uint32_t))  );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* preprocessing EC data for GPU ... ");
            preprocessDataForGPU(voxelEC, npeaks, segmentsPerBlock, offsetPerBlock, nvoxels);
            printf("\n");
        
            printf("\t* copying operator A to GPU (EC part) ... ");
            status = true;
            status = status && cudaCheck( cudaMemcpy(this->voxelEC,            voxelEC,              npeaks*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(this->orienEC,            orienEC,              npeaks*sizeof(uint16_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(this->segmentsPerBlockEC, segmentsPerBlock,     nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(this->offsetPerBlockEC,   offsetPerBlock,       nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            free(segmentsPerBlock);
            free(offsetPerBlock);
        }
        
        ~CudaLinearOperator(){
            cudaFree(voxelIC);
            cudaFree(fiberIC);
            cudaFree(orienIC);
            cudaFree(lengthIC);
            cudaFree(lutIC);
            cudaFree(segmentsPerBlockIC);
            cudaFree(offsetPerBlockIC);
            
            cudaFree(voxelEC);
            cudaFree(orienEC);
            cudaFree(lutEC);
            cudaFree(segmentsPerBlockEC);
            cudaFree(offsetPerBlockEC);
        
            cudaFree(lutISO);
        
            cudaFree(voxelICt);
            cudaFree(fiberICt);
            cudaFree(orienICt);
            cudaFree(lengthICt);
            cudaFree(fibersPerBlockICt);
            cudaFree(offsetPerBlockICt);
        
            cudaFree(x);
            cudaFree(y);
        
            printf("\t* reseting GPU ... ");
            bool status = true;
            status = status && cudaCheck( cudaDeviceReset() );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        }

        void setTransposeData(
            uint32_t*  voxelIDs,
            uint32_t*  fiberIDs,
            uint16_t*  orienIDs,
            float32_t* lengths,
            int nsegments)
        {
            bool status;
            uint32_t*  fibersPerBlock = (uint32_t*) malloc(nfibers*sizeof(uint32_t));
            uint32_t*  offsetPerBlock = (uint32_t*) malloc(nfibers*sizeof(uint32_t));
        
            preprocessDataForGPU(fiberIDs, nsegments, fibersPerBlock, offsetPerBlock, nfibers);
        
            printf("\t* extra memory for operator A' ... ");
            status = true;
            status = status && cudaCheck( cudaMalloc((void**)&(voxelICt),  nsegments*sizeof(uint32_t))  );
            status = status && cudaCheck( cudaMalloc((void**)&(fiberICt),  nsegments*sizeof(uint32_t))  );
            status = status && cudaCheck( cudaMalloc((void**)&(orienICt),  nsegments*sizeof(uint16_t))  );
            status = status && cudaCheck( cudaMalloc((void**)&(lengthICt), nsegments*sizeof(float32_t)) );
            status = status && cudaCheck( cudaMalloc((void**)&(fibersPerBlockICt), nfibers*sizeof(uint32_t)) );
            status = status && cudaCheck( cudaMalloc((void**)&(offsetPerBlockICt), nfibers*sizeof(uint32_t)) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            printf("\t* transfering memory for operator A' ... ");
            status = true;
            status = status && cudaCheck( cudaMemcpy(voxelICt,  voxelIDs, nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(fiberICt,  fiberIDs, nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(orienICt,  orienIDs, nsegments*sizeof(uint16_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(lengthICt, lengths,  nsegments*sizeof(float32_t), cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(fibersPerBlockICt, fibersPerBlock, nfibers*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            status = status && cudaCheck( cudaMemcpy(offsetPerBlockICt, offsetPerBlock, nfibers*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
            if (status) printf("[ OK ]\n");
            else        printf("[ ERROR ]\n");
        
            free(fibersPerBlock);
            free(offsetPerBlock);
        }

        void multiplyByX(float64_t* x, float64_t* y){
            // Copy vector x to the GPU
            cudaMemcpy(this->x, x, ncols*sizeof(double), cudaMemcpyHostToDevice);

            // Multiply IC part in the GPU
            multiply_Ax_ICpart<<<nvoxels, 1024>>>(voxelIC, fiberIC, orienIC, lengthIC, segmentsPerBlockIC, offsetPerBlockIC, lutIC, this->x, this->y);

            //cudaCheckKernel();

            // Multiply EC part in the GPU
            multiply_Ax_ECpart<<<nvoxels, 512>>>(voxelEC, orienEC, segmentsPerBlockEC, offsetPerBlockEC, lutEC, this->x, this->y);

            //cudaCheckKernel();

            // Multiply ISO part in the GPU
            multiply_Ax_ISOpart<<<nvoxels, 512>>>(lutISO, this->x, this->y);

            //cudaCheckKernel();

            // Copy back result to CPU
            cudaMemcpy(y, this->y, nrows*sizeof(double), cudaMemcpyDeviceToHost);
        }

        void multiplyByY(float64_t* v_in, float64_t* v_out){
        
            // Copy vector y to the GPU
            //cudaCheck( cudaMemset(gpu_x, 0, NUM_COLS*sizeof(float64_t)) );
            //cudaCheck( cudaMemcpy(gpu_x, x, NUM_COLS*sizeof(double), cudaMemcpyHostToDevice) );
            cudaCheck( cudaMemcpy(y, v_in, nrows*sizeof(double), cudaMemcpyHostToDevice) );
        
            // Multiply IC part in the GPU
            multiply_Aty_ICpart<<<nfibers, 512>>>(voxelICt, fiberICt, orienICt, lengthICt, fibersPerBlockICt, offsetPerBlockICt, lutIC, x, y);
        
            //cudaCheckKernel();//*/
        
            // Multiply EC part in the GPU
            multiply_Aty_ECpart<<<nvoxels, 512>>>(voxelEC, orienEC, segmentsPerBlockEC, offsetPerBlockEC, lutEC, x, y);
        
            //cudaCheckKernel();
        
            // Multiply ISO part in the GPU
            multiply_Aty_ISOpart<<<nvoxels, 512>>>(lutISO, x, y);
        
            //cudaCheckKernel();//*/
        
            // Copy back result to CPU
            cudaCheck( cudaMemcpy(v_out, x, ncols*sizeof(double), cudaMemcpyDeviceToHost) );
                
            /*printf("\n\n VECTOR X EC PART:\n");
            for(int i = NUM_FIBERS*NUM_RESFUNCIC; i < NUM_FIBERS*NUM_RESFUNCIC+20; i++)
                printf("%lf ", x[i]);
            printf("\n\n");//*/
        }
};

bool cudaCheck(cudaError_t cudaStatus){
    return cudaStatus == cudaSuccess;
}

void preprocessDataForGPU(uint32_t* data, int NUM_COMPARTMENTS, uint32_t* compartmentsPerBlock, uint32_t* offsetPerBlock, int NUM_BLOCKS){

    // fill arrays with zeros
    memset(compartmentsPerBlock, 0, NUM_BLOCKS * sizeof(uint32_t));
    memset(offsetPerBlock,       0, NUM_BLOCKS * sizeof(uint32_t));

    // count compartments per block
    for(int i = 0; i < NUM_COMPARTMENTS; i++)
        compartmentsPerBlock[data[i]]++;

    // calculate offset per block
    offsetPerBlock[0] = 0;
    for(int i = 1; i < NUM_BLOCKS; i++)
        offsetPerBlock[i] = offsetPerBlock[i-1] + compartmentsPerBlock[i-1];
}

/*CudaLinearOperator::CudaLinearOperator(
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
    int nballs)
{
    this->nvoxels = nvoxels;
    this->nfibers = nfibers;
    this->nrows = nvoxels * nsamples;
    this->ncols = nfibers*ndiameters + npeaks*nzeppelins + nvoxels*nballs;
    int size_lutic  = ndiameters*norientations*nsamples;
    int size_lutec  = nzeppelins*norientations*nsamples;
    int size_lutiso = nballs*nsamples;
    bool status;

    uint32_t* segmentsPerBlock = (uint32_t*) malloc(nvoxels*sizeof(uint32_t));
    uint32_t* offsetPerBlock   = (uint32_t*) malloc(nvoxels*sizeof(uint32_t));

    // copy constant values to GPU
    printf("\t* constant global values ... ");
    status = true;
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_VOXELS,       &nvoxels,       sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_FIBERS,       &nfibers,       sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_PEAKS,        &npeaks,        sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_ORIENTATIONS, &norientations, sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_SAMPLES,      &nsamples,      sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_DIAMETERS,    &ndiameters,    sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_ZEPPELINS,    &nzeppelins,    sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_BALLS,        &nballs,        sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_ROWS,         &nrows,         sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(NUM_COLS,         &ncols,         sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTIC,       &size_lutic,    sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTEC,       &size_lutec,    sizeof(int)) );
    status = status && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTISO,      &size_lutiso,   sizeof(int)) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");


    // alloc memory in GPU for vectors x and y
    printf("\t* memory for vectors x and y ... ");
    status = true;
    status = status && cudaCheck( cudaMalloc((void**)&(this->x), ncols*sizeof(float64_t)) );
    status = status && cudaCheck( cudaMalloc((void**)&(this->y), nrows*sizeof(float64_t)) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    // alloc GPU memory for segments
    printf("\t* memory for LUT (IC part) ... ");
    status = true;
    status = status && cudaCheck( cudaMalloc((void**)&(this->lutIC), size_lutic*sizeof(float32_t)) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* copying LUT in GPU (IC part) ... ");
    status = true;
    status = status && cudaCheck( cudaMemcpy(this->lutIC, lutIC, size_lutic*sizeof(float32_t), cudaMemcpyHostToDevice) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* allocating memory for LUT in GPU (EC part) ... ");
    status = cudaCheck( cudaMalloc((void**)&(this->lutEC), size_lutec*sizeof(float32_t)) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* copying LUT in GPU (EC part) ... ");
    status = cudaCheck( cudaMemcpy(this->lutEC, lutEC, size_lutec*sizeof(float32_t), cudaMemcpyHostToDevice) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* allocating memory for LUT in GPU (ISO part) ... ");
    status = cudaCheck( cudaMalloc((void**)&(this->lutISO), size_lutiso*sizeof(float32_t)) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* copying LUT in GPU (ISO part) ... ");
    status = cudaCheck( cudaMemcpy(this->lutISO, lutISO, size_lutiso*sizeof(float32_t), cudaMemcpyHostToDevice) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* preprocessing data for GPU ... ");
    preprocessDataForGPU(voxelIC, nsegments, segmentsPerBlock, offsetPerBlock, nvoxels);
    printf("\n");

    printf("\t* fiber segments memory allocation ... ");
    status = true;
    status = status && cudaCheck( cudaMalloc((void**)&(this->voxelIC),  nsegments*sizeof(uint32_t))  );
    status = status && cudaCheck( cudaMalloc((void**)&(this->fiberIC),  nsegments*sizeof(uint32_t))  );
    status = status && cudaCheck( cudaMalloc((void**)&(this->orienIC),  nsegments*sizeof(uint16_t))  );
    status = status && cudaCheck( cudaMalloc((void**)&(this->lengthIC), nsegments*sizeof(float32_t)) );
    status = status && cudaCheck( cudaMalloc((void**)&(this->segmentsPerBlockIC), nvoxels*sizeof(uint32_t)) );
    status = status && cudaCheck( cudaMalloc((void**)&(this->offsetPerBlockIC),   nvoxels*sizeof(uint32_t)) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* transfering fiber segments ... ");
    status = true;
    status = status && cudaCheck( cudaMemcpy(this->voxelIC,  voxelIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(this->fiberIC,  fiberIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(this->orienIC,  orienIC,  nsegments*sizeof(uint16_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(this->lengthIC, lengthIC, nsegments*sizeof(float32_t), cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(this->segmentsPerBlockIC, segmentsPerBlock, nvoxels*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(this->offsetPerBlockIC,   offsetPerBlock,   nvoxels*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    // ---------------------------------------- EC DATA ---------------------------------------- //
    printf("\t* allocating memory for operator A in GPU (EC part) ... ");
    status = true;
    status = status && cudaCheck( cudaMalloc((void**)&(this->voxelEC),  npeaks*sizeof(uint32_t)) );
    status = status && cudaCheck( cudaMalloc((void**)&(this->orienEC),  npeaks*sizeof(uint16_t)) );
    status = status && cudaCheck( cudaMalloc((void**)&(this->segmentsPerBlockEC), nvoxels*sizeof(uint32_t))  );
    status = status && cudaCheck( cudaMalloc((void**)&(this->offsetPerBlockEC),   nvoxels*sizeof(uint32_t))  );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* preprocessing EC data for GPU ... ");
    preprocessDataForGPU(voxelEC, npeaks, segmentsPerBlock, offsetPerBlock, nvoxels);
    printf("\n");

    printf("\t* copying operator A to GPU (EC part) ... ");
    status = true;
    status = status && cudaCheck( cudaMemcpy(this->voxelEC,            voxelEC,              npeaks*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(this->orienEC,            orienEC,              npeaks*sizeof(uint16_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(this->segmentsPerBlockEC, segmentsPerBlock,     nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(this->offsetPerBlockEC,   offsetPerBlock,       nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    free(segmentsPerBlock);
    free(offsetPerBlock);
}*/

/*CudaLinearOperator::~CudaLinearOperator(){
    cudaFree(voxelIC);
    cudaFree(fiberIC);
    cudaFree(orienIC);
    cudaFree(lengthIC);
    cudaFree(lutIC);
    cudaFree(segmentsPerBlockIC);
    cudaFree(offsetPerBlockIC);
    
    cudaFree(voxelEC);
    cudaFree(orienEC);
    cudaFree(lutEC);
    cudaFree(segmentsPerBlockEC);
    cudaFree(offsetPerBlockEC);

    cudaFree(lutISO);

    cudaFree(voxelICt);
    cudaFree(fiberICt);
    cudaFree(orienICt);
    cudaFree(lengthICt);
    cudaFree(fibersPerBlockICt);
    cudaFree(offsetPerBlockICt);

    cudaFree(x);
    cudaFree(y);

    printf("\t* reseting GPU ... ");
    bool status = true;
    status = status && cudaCheck( cudaDeviceReset() );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");
}*/

/*void CudaLinearOperator::setTransposeData(
    uint32_t*  voxelIDs,
    uint32_t*  fiberIDs,
    uint16_t*  orienIDs,
    float32_t* lengths,
    int nsegments)
{
    bool status;
    uint32_t*  fibersPerBlock = (uint32_t*) malloc(nfibers*sizeof(uint32_t));
    uint32_t*  offsetPerBlock = (uint32_t*) malloc(nfibers*sizeof(uint32_t));

    preprocessDataForGPU(fiberIDs, nsegments, fibersPerBlock, offsetPerBlock, nfibers);

    printf("\t* extra memory for operator A' ... ");
    status = true;
    status = status && cudaCheck( cudaMalloc((void**)&(voxelICt),  nsegments*sizeof(uint32_t))  );
    status = status && cudaCheck( cudaMalloc((void**)&(fiberICt),  nsegments*sizeof(uint32_t))  );
    status = status && cudaCheck( cudaMalloc((void**)&(orienICt),  nsegments*sizeof(uint16_t))  );
    status = status && cudaCheck( cudaMalloc((void**)&(lengthICt), nsegments*sizeof(float32_t)) );
    status = status && cudaCheck( cudaMalloc((void**)&(fibersPerBlockICt), nfibers*sizeof(uint32_t)) );
    status = status && cudaCheck( cudaMalloc((void**)&(offsetPerBlockICt), nfibers*sizeof(uint32_t)) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    printf("\t* transfering memory for operator A' ... ");
    status = true;
    status = status && cudaCheck( cudaMemcpy(voxelICt,  voxelIDs, nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(fiberICt,  fiberIDs, nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(orienICt,  orienIDs, nsegments*sizeof(uint16_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(lengthICt, lengths,  nsegments*sizeof(float32_t), cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(fibersPerBlockICt, fibersPerBlock, nfibers*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    status = status && cudaCheck( cudaMemcpy(offsetPerBlockICt, offsetPerBlock, nfibers*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    if (status) printf("[ OK ]\n");
    else        printf("[ ERROR ]\n");

    free(fibersPerBlock);
    free(offsetPerBlock);
}*/

__global__ void multiply_Ax_ICpart(
    uint32_t*  voxelIDs,
    uint32_t*  fiberIDs,
    uint16_t*  orienIDs,
    float32_t* lengths,
    uint32_t*  segmentsPerBlock,
    uint32_t*  offsetPerBlock,
    float32_t* lut,
    float64_t* x,
    float64_t* y){

    __shared__ float64_t shmem[1024];

    uint32_t bid = blockIdx.x;
    uint32_t tid = threadIdx.x;
    uint32_t gid = threadIdx.x / 512;
    uint32_t sid = threadIdx.x - 512*gid;

    shmem[tid] = 0.0;

    if(sid >= NUM_SAMPLES) return;

    uint32_t offset = offsetPerBlock[bid] + (segmentsPerBlock[bid]/2)*gid;
    uint32_t nsegments = segmentsPerBlock[bid]/2 + (segmentsPerBlock[bid]%2)*gid;

    //segment_t* segment = segments + offset;
    uint32_t*  voxel  = voxelIDs + offset;
    uint32_t*  fiber  = fiberIDs + offset;
    uint16_t*  orien  = orienIDs + offset;
    float32_t* length = lengths  + offset;

    float64_t sum = 0.0;

    for(int i = 0; i < nsegments; i++){
        int offset_lut = (*orien)*NUM_SAMPLES + sid;

        float64_t aux = 0.0;
        for(int j = 0; j < NUM_DIAMETERS; j++){
            aux += (double)(lut[offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES])*x[(*fiber) + j*NUM_FIBERS];
            //aux += tex1Dfetch(tex_lutIC, offset_lut + j*num_orientations*num_samples) * x[(*fiber) + j*num_fibers];
        }

        sum += aux * (*length);

        fiber++;
        orien++;
        length++;
    }

    shmem[tid] = sum;
    __syncthreads();

    if(tid < NUM_SAMPLES)
        y[(*voxel)*NUM_SAMPLES + sid] = sum + shmem[tid+512];
}

__global__ void multiply_Ax_ECpart(
    uint32_t*  voxelIDs,
    uint16_t*  orienIDs,
    uint32_t*  segmentsPerBlock,
    uint32_t*  offsetPerBlock,
    float32_t* lut,
    float64_t* x,
    float64_t* y)
{
    uint32_t bid = blockIdx.x;
    uint32_t tid = threadIdx.x;

    if(tid >= NUM_SAMPLES) return;

    uint32_t offset  = offsetPerBlock[bid];
    uint32_t nsegments = segmentsPerBlock[bid];

    //compartmentEC_t* excomp = excomps + offset;
    uint32_t* voxel = voxelIDs + offset;
    uint16_t* orien = orienIDs + offset;

    uint32_t target = NUM_FIBERS*NUM_DIAMETERS + offset;

    float64_t sum = 0.0;
    for(int i = 0; i < nsegments; i++){
        uint32_t offset_lut = (*orien)*NUM_SAMPLES + tid;

        for(int j = 0; j < NUM_ZEPPELINS; j++)
            sum += (double)(lut[offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES])*x[target + j*NUM_PEAKS + i];
            //sum += tex1Dfetch(tex_lutEC, offset_lut + j*num_orientations*num_samples) * x[target + j*num_excomps + i];

        orien++;
    }

    y[(*voxel)*NUM_SAMPLES + tid] += sum;
}

__global__ void multiply_Ax_ISOpart(
    float32_t* lut,
    float64_t* x,
    float64_t* y)
{
    uint32_t bid = blockIdx.x;
    uint32_t tid = threadIdx.x;

    if(tid >= NUM_SAMPLES) return;

    uint32_t target = NUM_FIBERS*NUM_DIAMETERS + NUM_PEAKS*NUM_ZEPPELINS + bid;

    float64_t sum = 0.0;
    for(int j = 0; j < NUM_BALLS; j++)
        sum += (double)(lut[j*NUM_SAMPLES + tid])*x[target + j*NUM_VOXELS];
        //sum += (double)(tex1Dfetch(tex_lutISO, j*num_samples + tid))*x[target + j*num_voxels];
        

    y[bid*NUM_SAMPLES + tid] += sum;
}

__global__ void multiply_Aty_ICpart(
    uint32_t*  voxelICt,
    uint32_t*  fiberICt,
    uint16_t*  orienICt,
    float32_t* lengthICt,
    uint32_t*  compartmentsPerBlock,
    uint32_t*  offsetPerBlock,
    float32_t* lut,
    float64_t* x,
    float64_t* y)
{
    __shared__ float64_t shmem[512];

    uint32_t bid = blockIdx.x;
    uint32_t tid = threadIdx.x;

    shmem[tid] = 0.0;

    if(tid >= NUM_SAMPLES) return;

    /*if(bid == 0 && tid == 0){
    for(int i = 0; i < 10; i++){
    printf("%d %d %d %f\n", voxelICt[i], fiberICt[i], orientICt[i], lengthICt[i]);
    }
    }
    else if(bid != 0) return;
    //__syncthreads();//*/

    uint32_t offset = offsetPerBlock[bid];
    uint32_t nsegments = offset + compartmentsPerBlock[bid];

    //segment_t* segment = segments + offset;
    uint32_t*  voxel  = voxelICt  + offset;
    uint32_t*  fiber  = fiberICt  + offset;
    uint16_t*  orien  = orienICt  + offset;
    float32_t* length = lengthICt + offset;
    //uint fiber = segment->fiber;

    for(int j = 0; j < NUM_DIAMETERS; j++){
        int offset_lut = j*NUM_ORIENTATIONS*NUM_SAMPLES + tid;

        float64_t sum = 0.0;
        //segment = segments + offset;
        voxel  = voxelICt  + offset;
        orien  = orienICt  + offset;
        length = lengthICt + offset;
        for(int i = offset; i < nsegments; i++){
            sum += ((float64_t)(*length)) *( (float64_t) lut[offset_lut + (*orien)*NUM_SAMPLES] )* y[(*voxel)*NUM_SAMPLES + tid];
            //sum += ((float64_t)(*length)) *( (float64_t) tex1Dfetch(tex_lutIC, offset_lut + (*orient)*num_samples) )* y[(*voxel)*num_samples + tid];
            //segment++;
            voxel++;
            //fiber++;
            orien++;
            length++;
        }

        shmem[tid] = sum;
        __syncthreads();

        if(tid < 256) shmem[tid] += shmem[tid + 256]; __syncthreads();
        if(tid < 128) shmem[tid] += shmem[tid + 128]; __syncthreads();
        if(tid <  64) shmem[tid] += shmem[tid +  64]; __syncthreads();
        if(tid <  32) shmem[tid] += shmem[tid +  32]; __syncthreads();
        if(tid <  16) shmem[tid] += shmem[tid +  16]; __syncthreads();
        if(tid <   8) shmem[tid] += shmem[tid +   8]; __syncthreads();
        if(tid <   4) shmem[tid] += shmem[tid +   4]; __syncthreads();
        //if(tid <   2) shmem[tid] += shmem[tid +   2]; __syncthreads();

        if(tid == 0) x[j*NUM_FIBERS + (*fiber)] = shmem[0] + shmem[1] + shmem[2] + shmem[3];

        __syncthreads();
    }
}

__global__ void multiply_Aty_ECpart(
    uint32_t*  voxelEC,
    uint16_t*  orienEC,
    uint32_t*  segmentsPerBlock,
    uint32_t*  offsetPerBlock,
    float32_t* lut,
    float64_t* x,
    float64_t* y)
{
    __shared__ float64_t shmem[512];

    uint32_t bid = blockIdx.x;
    uint32_t tid = threadIdx.x;

    shmem[tid] = 0.0;

    if(tid >= NUM_SAMPLES) return;

    uint32_t offset  = offsetPerBlock[bid];
    uint32_t ncompartments = segmentsPerBlock[bid] + offset;

    //compartmentEC_t* peak = peaks + offset;
    uint32_t* voxel = voxelEC + offset;
    uint16_t* orien = orienEC + offset;

    for(int j = 0; j < NUM_ZEPPELINS; j++){        
        uint32_t offset_lut = j*NUM_ORIENTATIONS*NUM_SAMPLES + tid;

        //peak = peaks + offset;
        voxel = voxelEC + offset;
        orien = orienEC + offset;
        for(int i = offset; i < ncompartments; i++){
            //shmem[tid] =( (float64_t)tex1Dfetch(tex_lutEC, (*orient)*num_samples + offset_lut) )* y[(*voxel)*num_samples + tid];
            shmem[tid] =( (float64_t)(lut[(*orien)*NUM_SAMPLES + offset_lut] ))* y[(*voxel)*NUM_SAMPLES + tid];
            __syncthreads();

            //if(bid == 0){
            //printf("%lf\n", lut[(peak->orientation)*num_samples + lut_offset] * y[(peak->voxel)*num_samples + tid]);

            if(tid < 256) shmem[tid] += shmem[tid + 256]; __syncthreads();
            if(tid < 128) shmem[tid] += shmem[tid + 128]; __syncthreads();
            if(tid <  64) shmem[tid] += shmem[tid +  64]; __syncthreads();
            if(tid <  32) shmem[tid] += shmem[tid +  32]; __syncthreads();
            if(tid <  16) shmem[tid] += shmem[tid +  16]; __syncthreads();
            if(tid <   8) shmem[tid] += shmem[tid +   8]; __syncthreads();
            if(tid <   4) shmem[tid] += shmem[tid +   4]; __syncthreads();
            if(tid <   2) shmem[tid] += shmem[tid +   2]; __syncthreads();

            if(tid == 0) x[NUM_FIBERS*NUM_DIAMETERS + j*NUM_PEAKS + i] = shmem[0] + shmem[1];
            //}

            //peak++;
            voxel++;
            orien++;
            __syncthreads();
        }
    }
} //*/

__global__ void multiply_Aty_ISOpart(float* lut, double* x, double* y){
    __shared__ double shmem[512];

    uint bid = blockIdx.x;
    uint tid = threadIdx.x;
    uint offset = NUM_FIBERS*NUM_DIAMETERS + NUM_PEAKS*NUM_ZEPPELINS + bid;

    shmem[tid] = 0.0;

    if(tid >= NUM_SAMPLES) return;

    for(int j = 0; j < NUM_BALLS; j++){
        shmem[tid] =( (float64_t) lut[j*NUM_SAMPLES + tid] )* y[bid*NUM_SAMPLES + tid];
        //shmem[tid] =( (float64_t) tex1Dfetch(tex_lutISO, j*num_samples + tid) )* y[bid*num_samples + tid];
        __syncthreads();

        if(tid < 256) shmem[tid] += shmem[tid + 256]; __syncthreads();
        if(tid < 128) shmem[tid] += shmem[tid + 128]; __syncthreads();
        if(tid <  64) shmem[tid] += shmem[tid +  64]; __syncthreads();
        if(tid <  32) shmem[tid] += shmem[tid +  32]; __syncthreads();
        if(tid <  16) shmem[tid] += shmem[tid +  16]; __syncthreads();
        if(tid <   8) shmem[tid] += shmem[tid +   8]; __syncthreads();
        if(tid <   4) shmem[tid] += shmem[tid +   4]; __syncthreads(); 

        if(tid == 0)
            x[offset + j*NUM_VOXELS] = shmem[0] + shmem[1] + shmem[2] + shmem[3];
    }
}//*/

/*void CudaLinearOperator::multiplyByX(float64_t* x, float64_t* y){

    // Copy vector x to the GPU
    cudaMemcpy(this->x, x, ncols*sizeof(double), cudaMemcpyHostToDevice);

    // Multiply IC part in the GPU
    multiply_Ax_ICpart<<<nvoxels, 1024>>>(voxelIC, fiberIC, orienIC, lengthIC, segmentsPerBlockIC, offsetPerBlockIC, lutIC, this->x, this->y);

    //cudaCheckKernel();

    // Multiply EC part in the GPU
    multiply_Ax_ECpart<<<nvoxels, 512>>>(voxelEC, orienEC, segmentsPerBlockEC, offsetPerBlockEC, lutEC, this->x, this->y);

    //cudaCheckKernel();

    // Multiply ISO part in the GPU
    multiply_Ax_ISOpart<<<nvoxels, 512>>>(lutISO, this->x, this->y);

    //cudaCheckKernel();

    // Copy back result to CPU
    cudaMemcpy(y, this->y, nrows*sizeof(double), cudaMemcpyDeviceToHost);
}*/

/*void CudaLinearOperator::multiplyByY(float64_t* v_in, float64_t* v_out){
        
    // Copy vector y to the GPU
    //cudaCheck( cudaMemset(gpu_x, 0, NUM_COLS*sizeof(float64_t)) );
    //cudaCheck( cudaMemcpy(gpu_x, x, NUM_COLS*sizeof(double), cudaMemcpyHostToDevice) );
    cudaCheck( cudaMemcpy(y, v_in, nrows*sizeof(double), cudaMemcpyHostToDevice) );

    // Multiply IC part in the GPU
    multiply_Aty_ICpart<<<nfibers, 512>>>(voxelICt, fiberICt, orienICt, lengthICt, fibersPerBlockICt, offsetPerBlockICt, lutIC, x, y);

    //cudaCheckKernel();//*/

    // Multiply EC part in the GPU
    multiply_Aty_ECpart<<<nvoxels, 512>>>(voxelEC, orienEC, segmentsPerBlockEC, offsetPerBlockEC, lutEC, x, y);

    //cudaCheckKernel();

    // Multiply ISO part in the GPU
    multiply_Aty_ISOpart<<<nvoxels, 512>>>(lutISO, x, y);

    //cudaCheckKernel();//*/

    // Copy back result to CPU
    cudaCheck( cudaMemcpy(v_out, x, ncols*sizeof(double), cudaMemcpyDeviceToHost) );
        
    /*printf("\n\n VECTOR X EC PART:\n");
    for(int i = NUM_FIBERS*NUM_RESFUNCIC; i < NUM_FIBERS*NUM_RESFUNCIC+20; i++)
        printf("%lf ", x[i]);
    printf("\n\n");//*/
}*/