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
            sum += (double)(lut[lut_offset + j*NUM_ORIENTATIONS*NUM_SAMPLES])*x[target + j*NUM_PEAKS + i];
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

void CudaLinearOperator::multiplyByX(float64_t* x, float64_t* y){

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

void CudaLinearOperator::multiplyByY(float64_t* y, float64_t* x){
        
    // Copy vector y to the GPU
    //cudaCheck( cudaMemset(gpu_x, 0, NUM_COLS*sizeof(float64_t)) );
    //cudaCheck( cudaMemcpy(gpu_x, x, NUM_COLS*sizeof(double), cudaMemcpyHostToDevice) );
    //cudaCheck( cudaMemcpy(gpu_y, y, NUM_ROWS*sizeof(double), cudaMemcpyHostToDevice) );

    // Multiply IC part in the GPU
    //multiply_Aty_ICpart<<<NUM_FIBERS, 512>>>(gpu_voxelICt, gpu_fiberICt, gpu_orientICt, gpu_lengthICt, gpu_segmentsPerBlockICt, gpu_offsetPerBlockICt, gpu_lutIC, gpu_x, gpu_y);

    //cudaCheckKernel();//*/

    // Multiply EC part in the GPU
    //multiply_Aty_ECpart<<<NUM_VOXELS, 512>>>(gpu_voxelEC, gpu_orientEC, gpu_segmentsPerBlockEC, gpu_offsetPerBlockEC, gpu_lutEC, gpu_x, gpu_y);

    //cudaCheckKernel();

    // Multiply ISO part in the GPU
    //multiply_Aty_ISOpart<<<NUM_VOXELS, 512>>>(gpu_lutISO, gpu_x, gpu_y);

    //cudaCheckKernel();//*/

    // Copy back result to CPU
    //cudaCheck( cudaMemcpy(x, gpu_x, NUM_COLS*sizeof(double), cudaMemcpyDeviceToHost) ); 
        
    /*printf("\n\n VECTOR X EC PART:\n");
    for(int i = NUM_FIBERS*NUM_RESFUNCIC; i < NUM_FIBERS*NUM_RESFUNCIC+20; i++)
        printf("%lf ", x[i]);
    printf("\n\n");//*/
}