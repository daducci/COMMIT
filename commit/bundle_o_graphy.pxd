# distutils: language = c++
# cython: language_level = 3

from libcpp cimport bool

cdef bool adapt_streamline(float [:,:] streamline, float [:,::1] to_RASMM, float [:] abc_to_RASMM, float [:,::1] to_VOXMM, float [:] abc_to_VOXMM,
                            int tempts, int pt_adapt, double m_variance, float [:, :, ::1] niiWM_img )
cdef double random_gaussian(double mov_var)
cdef void assign_random_gaussian_pt(double[:] out, int assign_ix, double mov_var)
cdef gaussian(int n, double mov_var)
cdef trk2dict_update(lut, index_list, diff_seg, float [:,:] streamlines, int* lengths, int* ptr_buff_size, double blur_core_extent, int Nx, int Ny, int Nz,
                    float Px, float Py, float Pz, int n_count,float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, float min_seg_len, float min_fiber_len, float max_fiber_len,
                    float* ptrPEAKS, double* ptrPeaksAffine, bool [:] flip_peaks, int Np, float vf_THR, float* _ptrMASK, float* _ptrISO, float* _ptrTDI,
                    float* ptrTractsAffine, unsigned short ndirs, short* ptrHashTable, 
                    unsigned char* pDict_TRK_kept, float* pDict_TRK_norm, unsigned int* pDict_IC_f, unsigned int* pDict_IC_v,
                    unsigned short* pDict_IC_o, float* pDict_IC_len, float* pDict_TRK_len, float* pDict_Tot_segm_len,
                    unsigned int*  pDict_EC_v, unsigned short* pDict_EC_o, int num_vox )