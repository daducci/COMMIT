#include "operator_withCUDA.cuh"

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

bool checkCompatibility(size_t required_mem, int gpu_id) {
    int num_gpus;
    cudaError_t cudaStatus;
    
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
        printf("\t* using GPU with ID %d... [ %s ]\n", gpu_id, gpu_properties.name);

        if (required_mem <= gpu_properties.totalGlobalMem) {
            printf("\t* using %.2f GB of total %.2f GB... [ OK ]\n", required_mem*1e-9, gpu_properties.totalGlobalMem*1e-9);
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

CudaLinearOperator::CudaLinearOperator(
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

    int fcall)
{
    this->nsegments = nsegments;
    this->nvoxels   = nvoxels;
    this->nfibers   = nfibers;
    this->nrows     = nvoxels * nsamples;
    this->ncols     = nfibers*ndiameters + npeaks*nzeppelins + nvoxels*nballs;

    if (fcall == 1) {
    int size_lutic  = ndiameters*norientations*nsamples;
    int size_lutec  = nzeppelins*norientations*nsamples;
    int size_lutiso = nballs*nsamples;

    size_t required_mem = 28*(size_t)nsegments + 6.0*(size_t)npeaks + 8.0*(size_t)nfibers + 16.0*(size_t)nvoxels + 4.0*((size_t)size_lutic + (size_t)size_lutec + (size_t)size_lutiso + (size_t)this->nrows + (size_t)this->ncols);
    checkCompatibility(required_mem, 0);

    // transfer constant values to the GPU
    printf("\t* constant values ... ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_VOXELS,       &nvoxels,       sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_FIBERS,       &nfibers,       sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_PEAKS,        &npeaks,        sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_ORIENTATIONS, &norientations, sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_SAMPLES,      &nsamples,      sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_DIAMETERS,    &ndiameters,    sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_ZEPPELINS,    &nzeppelins,    sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_BALLS,        &nballs,        sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_ROWS,         &nrows,         sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(NUM_COLS,         &ncols,         sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTIC,       &size_lutic,    sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTEC,       &size_lutec,    sizeof(int)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpyToSymbol(SIZE_LUTISO,      &size_lutiso,   sizeof(int)) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    // alloc memory in GPU for vectors x and y
    printf("\t* vectors x&y ... ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_x, ncols*sizeof(float64_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_y, nrows*sizeof(float64_t)) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    // pre-process data for GPU
    printf("\t* pre-processing ... ");
    cudaStatus = true;
    uint32_t* segmentsPerBlock = (uint32_t*) malloc(nvoxels*sizeof(uint32_t));
    uint32_t* offsetPerBlock   = (uint32_t*) malloc(nvoxels*sizeof(uint32_t));

    preprocessDataForGPU(voxelIC, nsegments, segmentsPerBlock, offsetPerBlock, nvoxels);

    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_segmentsPerBlockIC, nvoxels*sizeof(uint32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_offsetPerBlockIC,   nvoxels*sizeof(uint32_t)) );

    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_segmentsPerBlockIC, segmentsPerBlock, nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_offsetPerBlockIC,   offsetPerBlock,   nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice) );

    preprocessDataForGPU(voxelEC, npeaks, segmentsPerBlock, offsetPerBlock, nvoxels);

    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_segmentsPerBlockEC, nvoxels*sizeof(uint32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_offsetPerBlockEC,   nvoxels*sizeof(uint32_t)) );

    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_segmentsPerBlockEC, segmentsPerBlock, nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_offsetPerBlockEC,   offsetPerBlock,   nvoxels*sizeof(uint32_t), cudaMemcpyHostToDevice) );

    free(segmentsPerBlock);
    free(offsetPerBlock);
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    // alloc and transfer LUTs
    printf("\t* loading LUTs ... ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_lutIC, size_lutic*sizeof(float32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_lutIC, lutIC, size_lutic*sizeof(float32_t), cudaMemcpyHostToDevice) );

    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_lutEC,  size_lutec*sizeof(float32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_lutEC, lutEC, size_lutec*sizeof(float32_t), cudaMemcpyHostToDevice) );

    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_lutISO, size_lutiso*sizeof(float32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_lutISO, lutISO, size_lutiso*sizeof(float32_t), cudaMemcpyHostToDevice) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    // configure texture for LUTs
    tex_lutIC.addressMode[0] = cudaAddressModeBorder;
    tex_lutIC.addressMode[1] = cudaAddressModeBorder;
    tex_lutIC.filterMode = cudaFilterModePoint;
    tex_lutIC.normalized = false;

    tex_lutEC.addressMode[0] = cudaAddressModeBorder;
    tex_lutEC.addressMode[1] = cudaAddressModeBorder;
    tex_lutEC.filterMode = cudaFilterModePoint;
    tex_lutEC.normalized = false;

    tex_lutISO.addressMode[0] = cudaAddressModeBorder;
    tex_lutISO.addressMode[1] = cudaAddressModeBorder;
    tex_lutISO.filterMode = cudaFilterModePoint;
    tex_lutISO.normalized = false;

    printf("\t* linking LUTs to a texture memory ... ");
    cudaStatus = cudaStatus && cudaCheck( cudaBindTexture(NULL, tex_lutIC,  gpu_lutIC,  size_lutic  * sizeof(float32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaBindTexture(NULL, tex_lutEC,  gpu_lutEC,  size_lutec  * sizeof(float32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaBindTexture(NULL, tex_lutISO, gpu_lutISO, size_lutiso * sizeof(float32_t)) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    // alloc and transfer operator A
    printf("\t* A  operator... ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_voxelIC,  nsegments*sizeof(uint32_t))  );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_fiberIC,  nsegments*sizeof(uint32_t))  );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_orienIC,  nsegments*sizeof(uint16_t))  );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_lengthIC, nsegments*sizeof(float32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_voxelEC,  npeaks*sizeof(uint32_t))     );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_orienEC,  npeaks*sizeof(uint16_t))     );

    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_voxelIC,  voxelIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_fiberIC,  fiberIC,  nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_orienIC,  orienIC,  nsegments*sizeof(uint16_t),  cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_lengthIC, lengthIC, nsegments*sizeof(float32_t), cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_voxelEC,  voxelEC,  npeaks*sizeof(uint32_t),     cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_orienEC,  orienEC,  npeaks*sizeof(uint16_t),     cudaMemcpyHostToDevice) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");
    }

}

CudaLinearOperator::~CudaLinearOperator() {}

void CudaLinearOperator::destroy(){
    bool cudaStatus;

    printf("\n-> Deleting GPU memory:\n");

    printf("\t* deleting A...   ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_voxelIC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_fiberIC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_orienIC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_lengthIC) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_voxelEC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_orienEC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_segmentsPerBlockIC) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_offsetPerBlockIC)   );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_segmentsPerBlockEC) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_offsetPerBlockEC)   );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    printf("\t* deleting A'...  ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_TvoxelIC) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_TfiberIC) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_TorienIC) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_TlengthIC) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_TfibersPerBlockIC) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_ToffsetPerBlockIC) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    printf("\t* deleting x&y... ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_x) );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_y) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    printf("\t* deleting LUT... ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_lutIC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_lutEC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaFree(gpu_lutISO) );
    cudaStatus = cudaStatus && cudaCheck( cudaUnbindTexture(tex_lutIC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaUnbindTexture(tex_lutEC)  );
    cudaStatus = cudaStatus && cudaCheck( cudaUnbindTexture(tex_lutISO) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");

    printf("\t* reseting GPU... ");
    cudaStatus = true;
    cudaStatus = cudaStatus && cudaCheck( cudaDeviceReset() );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");
}

void CudaLinearOperator::setTransposeData(uint32_t*  voxelIDs,
                                          uint32_t*  fiberIDs,
                                          uint16_t*  orienIDs,
                                          float32_t* lengths)
{
    printf("\t* A' operator... ");
    cudaStatus = true;
    uint32_t*  fibersPerBlock = (uint32_t*) malloc(nfibers*sizeof(uint32_t));
    uint32_t*  offsetPerBlock = (uint32_t*) malloc(nfibers*sizeof(uint32_t));

    if(fibersPerBlock == NULL || offsetPerBlock == NULL) printf("problemas\n");

    preprocessDataForGPU(fiberIDs, nsegments, fibersPerBlock, offsetPerBlock, nfibers);

    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_TfibersPerBlockIC, nfibers*sizeof(uint32_t)) );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_ToffsetPerBlockIC, nfibers*sizeof(uint32_t)) );

    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_TfibersPerBlockIC, fibersPerBlock, nfibers*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_ToffsetPerBlockIC, offsetPerBlock, nfibers*sizeof(uint32_t),  cudaMemcpyHostToDevice) );

    free(fibersPerBlock);
    free(offsetPerBlock);

    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_TvoxelIC,  nsegments*sizeof(uint32_t))  );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_TfiberIC,  nsegments*sizeof(uint32_t))  );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_TorienIC,  nsegments*sizeof(uint16_t))  );
    cudaStatus = cudaStatus && cudaCheck( cudaMalloc((void**)&gpu_TlengthIC, nsegments*sizeof(float32_t)) );

    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_TvoxelIC,  voxelIDs, nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_TfiberIC,  fiberIDs, nsegments*sizeof(uint32_t),  cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_TorienIC,  orienIDs, nsegments*sizeof(uint16_t),  cudaMemcpyHostToDevice) );
    cudaStatus = cudaStatus && cudaCheck( cudaMemcpy(gpu_TlengthIC, lengths,  nsegments*sizeof(float32_t), cudaMemcpyHostToDevice) );
    if (cudaStatus) printf("[ OK ]\n");
    else            printf("[ CUDA ERROR ]\n");
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
    //cudaError_t cudaStatus;
    
    // Copy vector x to the GPU
    cudaMemcpy(gpu_x, v_in, ncols*sizeof(double), cudaMemcpyHostToDevice);
    /*if (cudaStatus != cudaSuccess) printf("\t* tranfering x to GPU ... [ ERROR ]: %s\n", cudaGetErrorString(cudaStatus));
    else                           printf("\t* tranfering x to GPU ... [   OK  ]\n");//*/

    // Multiply IC part in the GPU
    multiply_Ax_ICpart<<<nvoxels, 1024>>>(gpu_voxelIC, gpu_fiberIC, gpu_orienIC, gpu_lengthIC, gpu_segmentsPerBlockIC, gpu_offsetPerBlockIC, gpu_lutIC, gpu_x, gpu_y);

    //cudaCheckKernel();

    // Multiply EC part in the GPU
    multiply_Ax_ECpart<<<nvoxels, 512>>>(gpu_voxelEC, gpu_orienEC, gpu_segmentsPerBlockEC, gpu_offsetPerBlockEC, gpu_lutEC, gpu_x, gpu_y);

    //cudaCheckKernel();

    // Multiply ISO part in the GPU
    multiply_Ax_ISOpart<<<nvoxels, 512>>>(gpu_lutISO, gpu_x, gpu_y);

    //cudaCheckKernel();

    // Copy back result to CPU
    cudaMemcpy(v_out, gpu_y, nrows*sizeof(double), cudaMemcpyDeviceToHost);
    /*if (cudaStatus != cudaSuccess) printf("\t* tranfering y to CPU ... [ ERROR ]: %s\n", cudaGetErrorString(cudaStatus));
    else                           printf("\t* tranfering y to CPU ... [   OK  ]\n");//*/
}

void CudaLinearOperator::Tdot(float64_t* v_in, float64_t* v_out){
        
    //cudaError_t cudaStatus;
    // Copy vector y to the GPU
    //cudaCheck( cudaMemset(gpu_x, 0, NUM_COLS*sizeof(float64_t)) );
    //cudaCheck( cudaMemcpy(gpu_x, x, NUM_COLS*sizeof(double), cudaMemcpyHostToDevice) );
    cudaMemcpy(gpu_y, v_in, nrows*sizeof(double), cudaMemcpyHostToDevice);
    /*if (cudaStatus != cudaSuccess) printf("\t* tranfering y to GPU ... [ ERROR ]: %s\n", cudaGetErrorString(cudaStatus));
    else                           printf("\t* tranfering y to GPU ... [   OK  ]\n");//*/

    // Multiply IC part in the GPU
    multiply_Aty_ICpart<<<nfibers, 512>>>(gpu_TvoxelIC, gpu_TfiberIC, gpu_TorienIC, gpu_TlengthIC, gpu_TfibersPerBlockIC, gpu_ToffsetPerBlockIC, gpu_lutIC, gpu_x, gpu_y);

    //cudaCheckKernel();

    // Multiply EC part in the GPU
    multiply_Aty_ECpart<<<nvoxels, 512>>>(gpu_voxelEC, gpu_orienEC, gpu_segmentsPerBlockEC, gpu_offsetPerBlockEC, gpu_lutEC, gpu_x, gpu_y);

    //cudaCheckKernel();

    // Multiply ISO part in the GPU
    multiply_Aty_ISOpart<<<nvoxels, 512>>>(gpu_lutISO, gpu_x, gpu_y);

    //cudaCheckKernel();

    // Copy back result to CPU
    cudaMemcpy(v_out, gpu_x, ncols*sizeof(double), cudaMemcpyDeviceToHost);
    /*if (cudaStatus != cudaSuccess) printf("\t* tranfering x to CPU ... [ ERROR ]: %s\n", cudaGetErrorString(cudaStatus));
    else                           printf("\t* tranfering x to CPU ... [   OK  ]\n");//*/
        
    /*printf("\n\n VECTOR X EC PART:\n");
    for(int i = NUM_FIBERS*NUM_RESFUNCIC; i < NUM_FIBERS*NUM_RESFUNCIC+20; i++)
        printf("%lf ", x[i]);
    printf("\n\n");//*/
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
            //aux += (double)(lut[offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES])*x[(*fiber) + j*NUM_FIBERS];
            aux += tex1Dfetch(tex_lutIC, offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES) * x[(*fiber) + j*NUM_FIBERS];
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
            //sum += (double)(lut[offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES])*x[target + j*NUM_PEAKS + i];
            sum += tex1Dfetch(tex_lutEC, offset_lut + j*NUM_ORIENTATIONS*NUM_SAMPLES) * x[target + j*NUM_PEAKS + i];

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
        //sum += (double)(lut[j*NUM_SAMPLES + tid])*x[target + j*NUM_VOXELS];
        sum += (double)(tex1Dfetch(tex_lutISO, j*NUM_SAMPLES + tid))*x[target + j*NUM_VOXELS];
        

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
            //summ += ((float64_t)(*length)) *( (float64_t) lut[offset_lut + (*orien)*NUM_SAMPLES] )* y[(*voxel)*NUM_SAMPLES + tid];
            sum += ((float64_t)(*length)) *( (float64_t) tex1Dfetch(tex_lutIC, offset_lut + (*orien)*NUM_SAMPLES) )* y[(*voxel)*NUM_SAMPLES + tid];

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
            //shmem[tid] =( (float64_t)(lut[(*orien)*NUM_SAMPLES + offset_lut] ))* y[(*voxel)*NUM_SAMPLES + tid];
            shmem[tid] =( (float64_t)tex1Dfetch(tex_lutEC, (*orien)*NUM_SAMPLES + offset_lut) )* y[(*voxel)*NUM_SAMPLES + tid];
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
        //shmem[tid] =( (float64_t) lut[j*NUM_SAMPLES + tid] )* y[bid*NUM_SAMPLES + tid];
        shmem[tid] =( (float64_t) tex1Dfetch(tex_lutISO, j*NUM_SAMPLES + tid) )* y[bid*NUM_SAMPLES + tid];
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

