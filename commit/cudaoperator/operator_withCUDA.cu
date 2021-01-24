#include "operator_withCUDA.cuh"

// ====================================================
// Textures for LUT in the GPU
// ====================================================
texture<float32_t, 1, cudaReadModeElementType> tex_lutIC;
texture<float32_t, 1, cudaReadModeElementType> tex_lutEC;
texture<float32_t, 1, cudaReadModeElementType> tex_lutISO;


int checkCompatibility(int gpuID) {
    int gpuCount;
    cudaError_t cudaStatus;
    
    cudaStatus = cudaGetDeviceCount(&gpuCount);

    if (gpuCount <= 0 || gpuID >= gpuCount || cudaStatus != cudaSuccess) return 1;

    cudaStatus = cudaSetDevice(gpuID);

    if (cudaStatus != cudaSuccess) return 2;

    cudaDeviceProp gpuProperties;
    cudaStatus = cudaGetDeviceProperties(&gpuProperties, gpuID);

    if (cudaStatus != cudaSuccess) return 3;

    printf("\t* selected GPU...       [ %s ]\n",     gpuProperties.name);
    printf("\t* total memory...       [ %.2fGB ]\n", gpuProperties.totalGlobalMem*1e-9);
    printf("\t* compute capability... [ %d.%d ]\n",  gpuProperties.major, gpuProperties.minor);

    if(gpuProperties.major < 5) return 4;

    return 0;
}

void cudaCheckLastError()
{
    cudaError_t err = cudaGetLastError();

    if(err != cudaSuccess){
        printf("CUDA Error: %s\n", cudaGetErrorString(err));
        exit(-1);
    }
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

int CudaLinearOperator::setDictionary(uint32_t* voxelIC, uint32_t* fiberIC, uint16_t* orienIC, float32_t* lengthIC, uint32_t* voxelEC, uint16_t* orienEC){
    
    cudaError_t cudaStatus;

    uint32_t* segmentsPerBlock = (uint32_t*) malloc(nvoxels*sizeof(uint32_t));
    uint32_t* offsetPerBlock   = (uint32_t*) malloc(nvoxels*sizeof(uint32_t));

    if (segmentsPerBlock == NULL || offsetPerBlock == NULL) return -1;

    preprocessDataForGPU(voxelIC, nsegments, segmentsPerBlock, offsetPerBlock, nvoxels);

    cudaStatus = cudaMalloc((void**)&gpu_segmentsPerBlockIC, nvoxels*sizeof(uint32_t));
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_offsetPerBlockIC,   nvoxels*sizeof(uint32_t));
    if (cudaStatus != cudaSuccess) return 1;

    cudaStatus = cudaMemcpy(gpu_segmentsPerBlockIC, segmentsPerBlock, nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    cudaStatus = cudaMemcpy(gpu_offsetPerBlockIC,   offsetPerBlock,   nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;

    if (npeaks > 0){
        preprocessDataForGPU(voxelEC, npeaks, segmentsPerBlock, offsetPerBlock, nvoxels);

        cudaStatus = cudaMalloc((void**)&gpu_segmentsPerBlockEC, nvoxels*sizeof(uint32_t));
        if (cudaStatus != cudaSuccess) return 1;
        cudaStatus = cudaMalloc((void**)&gpu_offsetPerBlockEC,   nvoxels*sizeof(uint32_t));
        if (cudaStatus != cudaSuccess) return 1;

        cudaStatus = cudaMemcpy(gpu_segmentsPerBlockEC, segmentsPerBlock, nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) return 2;
        cudaStatus = cudaMemcpy(gpu_offsetPerBlockEC,   offsetPerBlock,   nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) return 2;
    }

    free(segmentsPerBlock);
    free(offsetPerBlock);

    // alloc IC part of the dictionary in GPU
    cudaStatus = cudaMalloc((void**)&gpu_voxelIC,  nsegments*sizeof(uint32_t)); 
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_fiberIC,  nsegments*sizeof(uint32_t)); 
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_orienIC,  nsegments*sizeof(uint16_t)); 
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_lengthIC, nsegments*sizeof(float32_t));
    if (cudaStatus != cudaSuccess) return 1;

    // transfer IC part of the dictionary to GPU
    cudaStatus = cudaMemcpy(gpu_voxelIC,  voxelIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    cudaStatus = cudaMemcpy(gpu_fiberIC,  fiberIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    cudaStatus = cudaMemcpy(gpu_orienIC,  orienIC,  nsegments*sizeof(uint16_t),  cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    cudaStatus = cudaMemcpy(gpu_lengthIC, lengthIC, nsegments*sizeof(float32_t), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;

    if (npeaks > 0){
        // alloc EC part of the dictionary in GPU
        cudaStatus = cudaMalloc((void**)&gpu_voxelEC,  npeaks*sizeof(uint32_t));
        if (cudaStatus != cudaSuccess) return 1;
        cudaStatus = cudaMalloc((void**)&gpu_orienEC,  npeaks*sizeof(uint16_t));
        if (cudaStatus != cudaSuccess) return 1;

        // transfer EC part of the dictionary to GPU
        cudaStatus = cudaMemcpy(gpu_voxelEC,  voxelEC,  npeaks*sizeof(uint32_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) return 2;
        cudaStatus = cudaMemcpy(gpu_orienEC,  orienEC,  npeaks*sizeof(uint16_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) return 2;
    }

    return 0;
}

int CudaLinearOperator::setTransposeDictionary(uint32_t* TvoxelIC, uint32_t* TfiberIC, uint16_t* TorienIC, float32_t* TlengthIC){
    
    cudaError_t cudaStatus;

    uint32_t*  fibersPerBlock = (uint32_t*) malloc(nfibers*sizeof(uint32_t));
    uint32_t*  offsetPerBlock = (uint32_t*) malloc(nfibers*sizeof(uint32_t));
    if(fibersPerBlock == NULL || offsetPerBlock == NULL) return -1;

    preprocessDataForGPU(TfiberIC, nsegments, fibersPerBlock, offsetPerBlock, nfibers);

    cudaStatus = cudaMalloc((void**)&gpu_TfibersPerBlockIC, nfibers*sizeof(uint32_t));
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_ToffsetPerBlockIC, nfibers*sizeof(uint32_t));
    if (cudaStatus != cudaSuccess) return 1;

    cudaStatus = cudaMemcpy(gpu_TfibersPerBlockIC, fibersPerBlock, nfibers*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    cudaStatus = cudaMemcpy(gpu_ToffsetPerBlockIC, offsetPerBlock, nfibers*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;

    free(fibersPerBlock);
    free(offsetPerBlock);

    cudaStatus = cudaMalloc((void**)&gpu_TvoxelIC,  nsegments*sizeof(uint32_t)) ;
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_TfiberIC,  nsegments*sizeof(uint32_t)) ;
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_TorienIC,  nsegments*sizeof(uint16_t)) ;
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_TlengthIC, nsegments*sizeof(float32_t));
    if (cudaStatus != cudaSuccess) return 1;

    cudaStatus = cudaMemcpy(gpu_TvoxelIC,  TvoxelIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    cudaStatus = cudaMemcpy(gpu_TfiberIC,  TfiberIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    cudaStatus = cudaMemcpy(gpu_TorienIC,  TorienIC,  nsegments*sizeof(uint16_t),  cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    cudaStatus = cudaMemcpy(gpu_TlengthIC, TlengthIC, nsegments*sizeof(float32_t), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) return 2;
    
    return 0;
}

int CudaLinearOperator::setKernels(float32_t* lutIC, float32_t* lutEC, float32_t* lutISO){

    cudaError_t cudaStatus;

    if (ndiameters > 0){
        cudaStatus = cudaMalloc((void**)&gpu_lutIC, size_lutic*sizeof(float32_t));
        if (cudaStatus != cudaSuccess) return 1;
        cudaStatus = cudaMemcpy(gpu_lutIC, lutIC, size_lutic*sizeof(float32_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) return 2;

        tex_lutIC.addressMode[0] = cudaAddressModeBorder;
        tex_lutIC.addressMode[1] = cudaAddressModeBorder;
        tex_lutIC.filterMode = cudaFilterModePoint;
        tex_lutIC.normalized = false;

        cudaStatus = cudaBindTexture(NULL, tex_lutIC,  gpu_lutIC,  size_lutic*sizeof(float32_t));
        if (cudaStatus != cudaSuccess) return 3;
    }

    if (nzeppelins > 0){
        cudaStatus = cudaMalloc((void**)&gpu_lutEC,  size_lutec*sizeof(float32_t));
        if (cudaStatus != cudaSuccess) return 1;
        cudaStatus = cudaMemcpy(gpu_lutEC, lutEC, size_lutec*sizeof(float32_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) return 2;

        tex_lutEC.addressMode[0] = cudaAddressModeBorder;
        tex_lutEC.addressMode[1] = cudaAddressModeBorder;
        tex_lutEC.filterMode = cudaFilterModePoint;
        tex_lutEC.normalized = false;

        cudaStatus = cudaBindTexture(NULL, tex_lutEC,  gpu_lutEC,  size_lutec*sizeof(float32_t));
        if (cudaStatus != cudaSuccess) return 3;
    }

    if (nballs > 0){
        cudaStatus = cudaMalloc((void**)&gpu_lutISO, size_lutiso*sizeof(float32_t));
        if (cudaStatus != cudaSuccess) return 1;
        cudaStatus = cudaMemcpy(gpu_lutISO, lutISO, size_lutiso*sizeof(float32_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) return 2;

        tex_lutISO.addressMode[0] = cudaAddressModeBorder;
        tex_lutISO.addressMode[1] = cudaAddressModeBorder;
        tex_lutISO.filterMode = cudaFilterModePoint;
        tex_lutISO.normalized = false;

        cudaStatus = cudaBindTexture(NULL, tex_lutISO, gpu_lutISO, size_lutiso*sizeof(float32_t));
        if (cudaStatus != cudaSuccess) return 3;
    }

    return 0;
}

int CudaLinearOperator::setVectors(){
    
    cudaError_t cudaStatus;

    cudaStatus = cudaMalloc((void**)&gpu_x, ncols*sizeof(float64_t));
    if (cudaStatus != cudaSuccess) return 1;
    cudaStatus = cudaMalloc((void**)&gpu_y, nrows*sizeof(float64_t));
    if (cudaStatus != cudaSuccess) return 1;
    
    return 0;
}

int CudaLinearOperator::setGlobals(){
    
    cudaError_t cudaStatus;

    cudaStatus = cudaMemcpyToSymbol(NUM_VOXELS,       &nvoxels,       sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_FIBERS,       &nfibers,       sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_PEAKS,        &npeaks,        sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_ORIENTATIONS, &norientations, sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_SAMPLES,      &nsamples,      sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_DIAMETERS,    &ndiameters,    sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_ZEPPELINS,    &nzeppelins,    sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_BALLS,        &nballs,        sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_ROWS,         &nrows,         sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(NUM_COLS,         &ncols,         sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(SIZE_LUTIC,       &size_lutic,    sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(SIZE_LUTEC,       &size_lutec,    sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    cudaStatus = cudaMemcpyToSymbol(SIZE_LUTISO,      &size_lutiso,   sizeof(int));
    if (cudaStatus != cudaSuccess) return -1;
    
    return 0;
}

CudaLinearOperator::CudaLinearOperator(int nsegments, int nvoxels, int nfibers, int npeaks, int norientations, int nsamples, int ndiameters, int nzeppelins, int nballs){

    this->nsegments = nsegments;
    this->nvoxels = nvoxels;
    this->nfibers = nfibers;
    this->npeaks = npeaks;
    this->norientations = norientations;
    this->nsamples = nsamples;
    this->ndiameters = ndiameters;
    this->nzeppelins = nzeppelins;   
    this->nballs = nballs;
    this->size_lutic = ndiameters*norientations*nsamples;
    this->size_lutec = nzeppelins*norientations*nsamples;
    this->size_lutiso = nballs*nsamples;
    this->nrows = nvoxels*nsamples;
    this->ncols = nfibers*ndiameters + npeaks*nzeppelins + nvoxels*nballs;
}

CudaLinearOperator::~CudaLinearOperator() {}

int CudaLinearOperator::destroy(){
    cudaError_t cudaStatus;    

    printf("\t* deleting A...   ");
    cudaStatus = cudaFree(gpu_voxelIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_fiberIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_orienIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_lengthIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_voxelEC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_orienEC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_segmentsPerBlockIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_offsetPerBlockIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_segmentsPerBlockEC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_offsetPerBlockEC);
    if (cudaStatus != cudaSuccess) return 5;

    printf("\t* deleting A'...  ");
    cudaStatus = cudaFree(gpu_TvoxelIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_TfiberIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_TorienIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_TlengthIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_TfibersPerBlockIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_ToffsetPerBlockIC);
    if (cudaStatus != cudaSuccess) return 5;

    printf("\t* deleting x&y... ");
    cudaStatus = cudaFree(gpu_x);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_y);
    if (cudaStatus != cudaSuccess) return 5;

    printf("\t* deleting LUT... ");
    cudaStatus = cudaFree(gpu_lutIC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_lutEC);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaFree(gpu_lutISO);
    if (cudaStatus != cudaSuccess) return 5;
    cudaStatus = cudaUnbindTexture(tex_lutIC);
    if (cudaStatus != cudaSuccess) return 6;
    cudaStatus = cudaUnbindTexture(tex_lutEC);
    if (cudaStatus != cudaSuccess) return 6;
    cudaStatus = cudaUnbindTexture(tex_lutISO);
    if (cudaStatus != cudaSuccess) return 6;

    printf("\t* reseting GPU... ");
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) return 7;

    return 0;
}

void cudaCheckKernel(){
    cudaError_t cudaStatus;
    
    cudaStatus = cudaGetLastError();
	if(cudaStatus != cudaSuccess)
        fprintf(stderr, "\t* kernel launch... [ ERROR ]: %s\n\n", cudaGetErrorString(cudaStatus));
    else
        printf("\t* kernel launch... [ OK ]\n");

    cudaStatus = cudaDeviceSynchronize();
	if(cudaStatus != cudaSuccess)
        fprintf(stderr, "\t* cudaDeviceSynchronize() after launching kernel... [ ERROR ]: %d\n", cudaStatus);
    else
        printf("\t* cudaDeviceSynchronize() after launching kernel... [ OK ]\n");
}

void CudaLinearOperator::dot(float64_t* v_in, float64_t* v_out){
    
    // Copy vector x to the GPU
    cudaMemcpy(gpu_x, v_in, ncols*sizeof(double), cudaMemcpyHostToDevice);
    //cudaCheckLastError();

    // Multiply IC part in the GPU
    multiply_Ax_ICpart<<<nvoxels, 1024>>>(gpu_voxelIC, gpu_fiberIC, gpu_orienIC, gpu_lengthIC, gpu_segmentsPerBlockIC, gpu_offsetPerBlockIC, gpu_lutIC, gpu_x, gpu_y);
    //cudaCheckLastError();

    // Multiply EC part in the GPU
    multiply_Ax_ECpart<<<nvoxels, 512>>>(gpu_voxelEC, gpu_orienEC, gpu_segmentsPerBlockEC, gpu_offsetPerBlockEC, gpu_lutEC, gpu_x, gpu_y);
    //cudaCheckLastError();

    // Multiply ISO part in the GPU
    multiply_Ax_ISOpart<<<nvoxels, 512>>>(gpu_lutISO, gpu_x, gpu_y);
    //cudaCheckLastError();

    // Copy back result to CPU
    cudaMemcpy(v_out, gpu_y, nrows*sizeof(double), cudaMemcpyDeviceToHost);
    //cudaCheckLastError();
}

void CudaLinearOperator::Tdot(float64_t* v_in, float64_t* v_out){
    
    // Copy vector y to the GPU
    cudaMemcpy(gpu_y, v_in, nrows*sizeof(double), cudaMemcpyHostToDevice);
    //cudaCheckLastError();

    // Multiply IC part in the GPU
    multiply_Aty_ICpart<<<nfibers, 512>>>(gpu_TvoxelIC, gpu_TfiberIC, gpu_TorienIC, gpu_TlengthIC, gpu_TfibersPerBlockIC, gpu_ToffsetPerBlockIC, gpu_lutIC, gpu_x, gpu_y);
    //cudaCheckLastError();

    // Multiply EC part in the GPU
    multiply_Aty_ECpart<<<nvoxels, 512>>>(gpu_voxelEC, gpu_orienEC, gpu_segmentsPerBlockEC, gpu_offsetPerBlockEC, gpu_lutEC, gpu_x, gpu_y);
    //cudaCheckLastError();

    // Multiply ISO part in the GPU
    multiply_Aty_ISOpart<<<nvoxels, 512>>>(gpu_lutISO, gpu_x, gpu_y);
    //cudaCheckLastError();

    // Copy back result to CPU
    cudaMemcpy(v_out, gpu_x, ncols*sizeof(double), cudaMemcpyDeviceToHost);
    //cudaCheckLastError();
}

// ------------------------------------------------------- KERNELS ------------------------------------------------------- //
__global__ void multiply_Ax_ICpart(uint32_t*  voxelIDs,
                                   uint32_t*  fiberIDs,
                                   uint16_t*  orienIDs,
                                   float32_t* lengths,
                                   uint32_t*  segmentsPerBlock,
                                   uint32_t*  offsetPerBlock,
                                   float32_t* lut,
                                   float64_t* x,
                                   float64_t* y)
{
    __shared__ float64_t shmem[1024];

    uint32_t bid = blockIdx.x;
    uint32_t tid = threadIdx.x;
    uint32_t gid = threadIdx.x / 512;
    uint32_t sid = threadIdx.x - 512*gid;

    shmem[tid] = 0.0;

    if(sid >= NUM_SAMPLES) return;

    uint32_t offset = offsetPerBlock[bid] + (segmentsPerBlock[bid]/2)*gid;
    uint32_t nsegments = segmentsPerBlock[bid]/2 + (segmentsPerBlock[bid]%2)*gid;

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
            //aux += tex1Dfetch(tex_lutIC, offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES) * x[(*fiber) + j*NUM_FIBERS];
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

    uint32_t* voxel = voxelIDs + offset;
    uint16_t* orien = orienIDs + offset;

    uint32_t target = NUM_FIBERS*NUM_DIAMETERS + offset;

    float64_t sum = 0.0;
    for(int i = 0; i < nsegments; i++){
        uint32_t offset_lut = (*orien)*NUM_SAMPLES + tid;

        for(int j = 0; j < NUM_ZEPPELINS; j++)
            sum += (double)(lut[offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES])*x[target + j*NUM_PEAKS + i];
            //sum += tex1Dfetch(tex_lutEC, offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES) * x[target + j*NUM_PEAKS + i];

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
        //sum += (double)(tex1Dfetch(tex_lutISO, j*NUM_SAMPLES + tid))*x[target + j*NUM_VOXELS];
        

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

    uint32_t offset = offsetPerBlock[bid];
    uint32_t nsegments = offset + compartmentsPerBlock[bid];

    uint32_t*  voxel  = voxelICt  + offset;
    uint32_t*  fiber  = fiberICt  + offset;
    uint16_t*  orien  = orienICt  + offset;
    float32_t* length = lengthICt + offset;

    for(int j = 0; j < NUM_DIAMETERS; j++){
        int offset_lut = j*NUM_ORIENTATIONS*NUM_SAMPLES + tid;

        float64_t sum = 0.0;
        voxel  = voxelICt  + offset;
        orien  = orienICt  + offset;
        length = lengthICt + offset;
        for(int i = offset; i < nsegments; i++){
            sum += ((float64_t)(*length)) *( (float64_t) lut[offset_lut + (*orien)*NUM_SAMPLES] )* y[(*voxel)*NUM_SAMPLES + tid];
            //sum += ((float64_t)(*length)) *( (float64_t) tex1Dfetch(tex_lutIC, offset_lut + (*orien)*NUM_SAMPLES) )* y[(*voxel)*NUM_SAMPLES + tid];

            voxel++;
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

    uint32_t* voxel = voxelEC + offset;
    uint16_t* orien = orienEC + offset;

    for(int j = 0; j < NUM_ZEPPELINS; j++){        
        uint32_t offset_lut = j*NUM_ORIENTATIONS*NUM_SAMPLES + tid;

        voxel = voxelEC + offset;
        orien = orienEC + offset;
        for(int i = offset; i < ncompartments; i++){
            shmem[tid] =( (float64_t)(lut[(*orien)*NUM_SAMPLES + offset_lut] ))* y[(*voxel)*NUM_SAMPLES + tid];
            //shmem[tid] =( (float64_t)tex1Dfetch(tex_lutEC, (*orien)*NUM_SAMPLES + offset_lut) )* y[(*voxel)*NUM_SAMPLES + tid];
            __syncthreads();

            if(tid < 256) shmem[tid] += shmem[tid + 256]; __syncthreads();
            if(tid < 128) shmem[tid] += shmem[tid + 128]; __syncthreads();
            if(tid <  64) shmem[tid] += shmem[tid +  64]; __syncthreads();
            if(tid <  32) shmem[tid] += shmem[tid +  32]; __syncthreads();
            if(tid <  16) shmem[tid] += shmem[tid +  16]; __syncthreads();
            if(tid <   8) shmem[tid] += shmem[tid +   8]; __syncthreads();
            if(tid <   4) shmem[tid] += shmem[tid +   4]; __syncthreads();
            if(tid <   2) shmem[tid] += shmem[tid +   2]; __syncthreads();

            if(tid == 0) x[NUM_FIBERS*NUM_DIAMETERS + j*NUM_PEAKS + i] = shmem[0] + shmem[1];

            voxel++;
            orien++;
            __syncthreads();
        }
    }
}

__global__ void multiply_Aty_ISOpart(float* lut, double* x, double* y){
    __shared__ double shmem[512];

    uint bid = blockIdx.x;
    uint tid = threadIdx.x;
    uint offset = NUM_FIBERS*NUM_DIAMETERS + NUM_PEAKS*NUM_ZEPPELINS + bid;

    shmem[tid] = 0.0;

    if(tid >= NUM_SAMPLES) return;

    for(int j = 0; j < NUM_BALLS; j++){
        shmem[tid] =( (float64_t) lut[j*NUM_SAMPLES + tid] )* y[bid*NUM_SAMPLES + tid];
        //shmem[tid] =( (float64_t) tex1Dfetch(tex_lutISO, j*NUM_SAMPLES + tid) )* y[bid*NUM_SAMPLES + tid];
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
}

