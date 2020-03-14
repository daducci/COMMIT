#include "gpumanager.cuh"

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

/*
__dual__ segment::segment() {}

__dual__ segment::~segment() {}
//*/

CudaLinearOperator::CudaLinearOperator(
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
    int nrows = nvoxels * nsamples;
    int ncols = nfibers*ndiameters + npeaks*nzeppelins + nvoxels*nballs;
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

CudaLinearOperator::~CudaLinearOperator(){
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

    cudaFree(x);
    cudaFree(y);
}