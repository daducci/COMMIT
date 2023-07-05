#!python
# cython: language_level=3, c_string_type=str, c_string_encoding=ascii, boundscheck=False, wraparound=False, profile=False

import cython
import itertools
from multiprocessing import Pool
import numpy as np
cimport numpy as np
from collections import defaultdict
from libc.math cimport sqrt, log
from libc.math cimport round as cround
from dipy.tracking.streamline import set_number_of_points
from dicelib.streamline import smooth

import os
import tqdm
from rdp import rdp

from libcpp cimport bool
from libc.stdlib cimport rand, srand, RAND_MAX
from libc.time cimport time
from libc.stdio cimport printf


from amico.util import LOG, NOTE, WARNING, ERROR

srand(time(NULL)) # seed for random generation

cdef extern from "trk2dictionary_c.cpp":
    int trk2dictionary_update_run(
        int* new_seg_len, float* ptr_start, int* lenghts, int* ptrIdx_list, int Nx, int Ny, int Nz,
        float Px, float Py, float Pz, int n_count,
        float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, float min_seg_len, float min_fiber_len, float max_fiber_len,
        float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz, float* _ptrMASK, float* _ptrISO, float* ptrTDI,
        double* ptrPeaksAffine, int nReplicas, double* ptrBlurRho, double* ptrBlurAngle, double* ptrBlurWeights, float* ptrTractsAffine,
        unsigned short ndirs, short* ptrHashTable, unsigned char* pDict_TRK_kept, float* pDict_TRK_norm, unsigned int* pDict_IC_f,
        unsigned int* pDict_IC_v, unsigned short* pDict_IC_o, float* pDict_IC_len, float* pDict_TRK_len, float* pDict_Tot_segm_len,
        unsigned int*  pDict_EC_v, unsigned short* pDict_EC_o
        ) nogil



@cython.boundscheck(True)
cdef double random_uniform():
    cdef double r = rand()
    return r / RAND_MAX


@cython.boundscheck(True)
cdef double random_gaussian(double mov_var):
    cdef double x1, x2, x3, w
    cdef double stdDev = mov_var
    w = 2.0
    while (w >= 1.0):
        x1 = 2.0 * random_uniform() - 1.0
        x2 = 2.0 * random_uniform() - 1.0
        x3 = 2.0 * random_uniform() - 1.0
        w = x1 * x1 + x2 * x2 + x3 * x3

    w = ((-2.0 * log(w)) / w) ** 0.5
    return x1 * w * stdDev


cdef void assign_random_gaussian_pt(double[:] out, int assign_ix, double mov_var):
    cdef double x1, x2, x3, w
    cdef double stdDev = mov_var
    w = 2.0
    while (w >= 1.0):
        x1 = 2.0 * random_uniform() - 1.0
        x2 = 2.0 * random_uniform() - 1.0
        x3 = 2.0 * random_uniform() - 1.0
        w = x1 * x1 + x2 * x2 + x3 * x3

    w = sqrt((-2.0 * log(w)) / w)
    out[assign_ix] = stdDev * x1 * w
    out[assign_ix + 1] = stdDev * x2 * w
    out[assign_ix + 2] = stdDev * x3 * w


cdef gaussian(int n, double mov_var):
    cdef int i
    cdef double[:] result = np.zeros(n, dtype='f8', order='C')
    for i in range(n // 2):  # Int division ensures trailing index if n is odd.
        assign_random_gaussian_pt(result, i * 2, mov_var)
    if n % 2 == 1:
        result[n - 1] = random_gaussian(mov_var)

    return result


cdef inline int randint(int lower, int upper) nogil:
    return rand() % (upper - lower + 1)


cpdef smooth_fib(streamlines, lengths, int n_count):
    cdef float [:,:] streamline_out =  np.ascontiguousarray( np.zeros( (3*10000,n_count) ).astype(np.float32) )
    trk_fiber_out = []
    lengths_out = np.zeros(n_count, dtype=np.int32)

    for f in range(n_count):
        streamline, n = smooth( streamlines[f], lengths[f], 1, 1 )
        lengths_out[f] = n
        trk_fiber_out.append( streamline )
    streamlines_out = np.vstack([s for s in trk_fiber_out])
    lengths_out = np.ascontiguousarray( lengths_out).astype(np.int32)
    return streamlines_out, lengths_out


cpdef smooth_final(streamlines, lengths, int n_count):
    cdef float [:,:] streamline_out =  np.ascontiguousarray( np.zeros( (3*10000,n_count) ).astype(np.float32) )
    trk_fiber_out = []
    lengths_out = np.zeros(n_count, dtype=np.int32)

    for f in range(n_count):
        streamline, n = smooth( streamlines[f], lengths[f], 1, 1 )
        lengths_out[f] = n
        trk_fiber_out.append( streamline )

    return trk_fiber_out



cdef float [:] apply_affine(float [:] pt, float [:,::1] M,
                                float [:] abc, float [:] pt_aff) nogil:

    pt_aff[0] = ((pt[0]*M[0,0] + pt[1]*M[1,0] + pt[2]*M[2,0]) + abc[0])
    pt_aff[1] = ((pt[0]*M[0,1] + pt[1]*M[1,1] + pt[2]*M[2,1]) + abc[1])
    pt_aff[2] = ((pt[0]*M[0,2] + pt[1]*M[1,2] + pt[2]*M[2,2]) + abc[2])

    return pt_aff



cdef bool adapt_streamline( float [:,:] streamline, float [:,::1] to_RASMM, float [:] abc_to_RASMM,
                            float [:,::1] to_VOXMM, float [:] abc_to_VOXMM, int tempts, int pt_adapt,
                            double m_variance, float [:, :, ::1] niiWM_img):
    """Compute the length of a streamline.

    Parameters
    ----------
    streamline : Nx3 numpy array
        The streamline data
    n : int
        Writes first n points of the streamline. If n<=0 (default), writes all points

    Returns
    -------
    length : double
        Length of the streamline in mm
    """

    cdef:
        int n = streamline.shape[0]
        int i,j, vox_x, vox_y, vox_z
        double [:] random_displ
        bool goodMove = False
        float length = 0.0
        int choose_pt = randint(0,n)
        float [:] temp_pt = np.zeros(3, dtype='f4', order='C')
        float [:] rand_pt = np.zeros(3, dtype='f4', order='C')

    if pt_adapt==0:
        rand_pt = streamline[choose_pt,:]
        random_displ = np.array(gaussian(3, m_variance))
        for i in xrange(tempts):
            temp_pt[:] = apply_affine(rand_pt, to_VOXMM, abc_to_VOXMM, temp_pt)
            temp_pt[0] += random_displ[0]
            temp_pt[1] += random_displ[1]
            temp_pt[2] += random_displ[2]

            vox_x = <int> cround(temp_pt[0])
            vox_y = <int> cround(temp_pt[1])
            vox_z = <int> cround(temp_pt[2])

            if ( niiWM_img[vox_x, vox_y, vox_z] != 0 ):
                streamline[choose_pt,:] = apply_affine(temp_pt, to_RASMM, abc_to_RASMM, temp_pt)
                break
        if i<tempts-1:
            goodMove = True
    else:
        for j in xrange(n):
            for i in xrange(tempts):
                random_displ = np.array(gaussian(3, m_variance))
                temp_pt[:] = apply_affine(streamline[j,:], to_VOXMM, abc_to_VOXMM, temp_pt)

                temp_pt[0] += random_displ[0]
                temp_pt[1] += random_displ[1]
                temp_pt[2] += random_displ[2]

                vox_x = <int>  cround(temp_pt[0])
                vox_y = <int>  cround(temp_pt[1])
                vox_z = <int>  cround(temp_pt[2])

                # length += sqrt( (ptr[3]-ptr[0])**2 + (ptr[4]-ptr[1])**2 + (ptr[5]-ptr[2])**2 )
                if ( niiWM_img[vox_x, vox_y, vox_z] != 0 ):
                    streamline[j] = apply_affine(temp_pt, to_RASMM, abc_to_RASMM, temp_pt)
                    break

            if i<tempts-1:
                goodMove = True
    return goodMove


cdef trk2dict_update(lut, index_list, diff_seg, float [:,:] streamlines, int* ptrlengths, int* ptr_buff_size, double blur_core_extent, int Nx, int Ny, int Nz,
                    float Px, float Py, float Pz, int n_count,float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, float min_seg_len, float min_fiber_len, float max_fiber_len,
                    float* ptrPEAKS, double* ptrPeaksAffine, bool [:] flip_peaks, int Np, float vf_THR, float* _ptrMASK, float* _ptrISO, float* _ptrTDI,
                    float* ptrTractsAffine, unsigned short ndirs, short* ptrHashTable, 
                    unsigned char* pDict_TRK_kept, float* pDict_TRK_norm, unsigned int* pDict_IC_f, unsigned int* pDict_IC_v,
                    unsigned short* pDict_IC_o, float* pDict_IC_len, float* pDict_TRK_len, float* pDict_Tot_segm_len,
                    unsigned int*  pDict_EC_v, unsigned short* pDict_EC_o, int num_vox):

    cdef:
        int i,j, pt_start, i_s
        float* ptr_start = &streamlines[0,0]
        double [:] blurRho
        double [:] blurAngle
        double [:] blurWeights
        bool [:] blurApplyTo
        int nReplicas, seg
        int [:] new_seg_len = np.zeros(len(index_list), dtype=np.int32)
        int* pt_new_seg = &new_seg_len[0]
        float blur_gauss_extent = 0.25
        float blur_spacing = 0.25
        float blur_gauss_min = 0.1
        float blur_sigma
        int [:] idx_list = np.array(index_list, np.int32)
        int* ptrIdx_list = &idx_list[0]


    if (blur_gauss_extent==0 and blur_core_extent==0) or (blur_spacing==0) :
        nReplicas = 1
        blurRho = np.array( [0.0], np.double )
        blurAngle = np.array( [0.0], np.double )
        blurWeights = np.array( [1], np.double )
    else:
        tmp = np.arange(0,blur_core_extent+blur_gauss_extent+1e-6,blur_spacing)
        tmp = np.concatenate( (tmp,-tmp[1:][::-1]) )
        x, y = np.meshgrid( tmp, tmp )
        r = np.sqrt( x*x + y*y )
        idx = (r <= blur_core_extent+blur_gauss_extent)
        blurRho = r[idx]
        blurAngle = np.arctan2(y,x)[idx]
        nReplicas = blurRho.size

        blurWeights = np.empty( nReplicas, np.double  )
        if blur_gauss_extent == 0 :
            blurWeights[:] = 1.0
        else:
            blur_sigma = blur_gauss_extent / np.sqrt( -2.0 * np.log( blur_gauss_min ) )
            for i in xrange(nReplicas):
                if blurRho[i] <= blur_core_extent :
                    blurWeights[i] = 1.0
                else:
                    blurWeights[i] = np.exp( -(blurRho[i] - blur_core_extent)**2 / (2.0*blur_sigma**2) )

    trk2dictionary_update_run(pt_new_seg, ptr_start, ptrlengths, ptrIdx_list, Nx, Ny, Nz, Px, Py, Pz, n_count,
                            fiber_shiftX, fiber_shiftY, fiber_shiftZ, min_seg_len, min_fiber_len, max_fiber_len,
                            ptrPEAKS, Np, vf_THR, -1 if flip_peaks[0] else 1, -1 if flip_peaks[1] else 1, -1 if flip_peaks[2] else 1,
                            _ptrMASK, _ptrISO, _ptrTDI, ptrPeaksAffine, nReplicas, &blurRho[0], &blurAngle[0], &blurWeights[0],
                            ptrTractsAffine, ndirs, ptrHashTable,
                            pDict_TRK_kept, pDict_TRK_norm, pDict_IC_f, pDict_IC_v, pDict_IC_o, pDict_IC_len, pDict_TRK_len,
                            pDict_Tot_segm_len, pDict_EC_v, pDict_EC_o
                            );
    # num_seg = 0
    # print(f"fib idx {int(ptrIdx_list[0])}, num of segments {int(pt_new_seg[0])}")
    # for i in range(int(pt_new_seg[0])+1):
    #     print(pDict_IC_f[i])

    # pt_start = ptr_buff_size[0]
    # for i in xrange(len(index_list)):
    #     for _ in xrange(new_seg_len[i]):
    #         pDict_IC_v[pt_start] = lut[pDict_IC_v[pt_start]]
    #         segm_dict[index_list[i]].extend([pDict_IC_f[pt_start]])
    #         pt_start += 1
    # diff_seg = diff_seg - new_seg_len.sum()
    # ptr_buff_size[0] += diff_seg
    temp_i_s = 0
    for i,k in enumerate(index_list):
        # new_s = []
        # print(f"segments for fib {k}: {new_seg_len[i]}")
        for i_s in xrange(new_seg_len[i]):
            # new_s.extend( [pDict_IC_f[i_s]] )
            # if new_seg_len[i]<10: 
            # print(f"pDict_IC_v no lut: {pDict_IC_v[temp_i_s]}")
            pDict_IC_v[temp_i_s] = lut[pDict_IC_v[temp_i_s]]
            # if new_seg_len[i]<10:
            # print(f"pDict_IC_v with lut: {pDict_IC_v[temp_i_s]}")
            temp_i_s += 1

        # segm_dict[k] = new_s
    # print(f"number old segments: {diff_seg}")
    # print(f"number new segments: {temp_i_s}\n")
    diff_seg_temp = diff_seg - temp_i_s
    ptr_buff_size[0] += diff_seg_temp


def create_prop_dict():
    prop_Dict       =   {}
    prop_Dict["MV"] =   np.arange(0, 25, 1)
    prop_Dict["A"]  =   np.arange(25, 50, 1)
    prop_Dict["B"]  =   np.arange(50, 75, 1)
    prop_Dict["K"]  =   np.arange(75, 100, 1)
    return prop_Dict


def compute_temp_schedule(opt_params):
    startTemp   = opt_params["startTemp"]
    endTemp     = opt_params["endTemp"]
    MAX_ITER_1  = opt_params["MAX_ITER_1"]
    SA_schedule = []
    it          = np.arange(MAX_ITER_1 + 1)
    SA_schedule = [-np.log(endTemp/startTemp)/i for i in it]
    # SA_schedule = np.repeat(1, MAX_ITER)
    return SA_schedule


def tractogram_2spline(input_tractogram):
    with Pool(28) as p:
        Sft_prop = list(tqdm.tqdm(p.imap(spline_repr, input_tractogram), total=len(input_tractogram)))
        return Sft_prop


def spline_repr(Fb, smooth=0.3):
    # reduced = rdp(Fb, epsilon=smooth, algo="iter")
    # if len(reduced) < 6 or np.isnan(np.sum(reduced)):
    reduced = set_number_of_points(Fb, 7)
    return reduced


def streamline2spline(input_set_streamlines, smth=0.3, parallel=False):
    if parallel:
        input_set_spline = tractogram_2spline(input_set_streamlines)
    else:
        input_set_spline = [spline_repr(Fb, smooth=smth) for Fb in tqdm.tqdm(input_set_streamlines)]
    return input_set_spline

def streamline_to_voxmm(streamlines, affine):
    streamlines_voxmm = []
    for streamline in streamlines:
        streamline_voxmm = []
        for point in streamline:
            point_voxmm = np.dot(affine, np.array([point[0], point[1], point[2], 1]))[:3]
            streamline_voxmm.append(point_voxmm)
        streamlines_voxmm.append(streamline_voxmm)
    return streamlines_voxmm


def compute_assignments(dictionary):
    # os.system( f'tck2connectome -force -symmetric -assignment_radial_search 2 -out_assignments {input_tractogram_filename}_fibers_assignment.txt {input_tractogram_filename} {gm_filename} {input_tractogram_filename}_connectome.csv' )
    # os.system( f'dice_tractogram_cluster.py {input_tractogram_filename} --save_assignments {input_tractogram_filename}_fibers_assignment.txt --atlas {gm_filename} --reference {gm_filename}' )

    file_assignments = os.path.join(dictionary['path_out'], f"{dictionary['filename_tractogram']}_clustered_thr_{float(dictionary['blur_clust_thr'][0])}_assignments.txt")
    os.system(f"dice_tractogram_cluster.py {dictionary['filename_tractogram']} --reference {dictionary['atlas']} "+
        f"--clust_threshold {dictionary['blur_clust_thr'][0]} --atlas {dictionary['atlas']} "+
        f"--n_threads {dictionary['nthreads']} --save_assignments {file_assignments} -f -v")
    
    print(f"Assignments file: {file_assignments}")
    asgn = np.loadtxt(file_assignments, dtype =float).astype(int)
    asgn = [tuple(sorted(c)) for c in asgn]
    connections = list(conn for conn,_ in itertools.groupby(asgn))
    connections_dict = dict((idx, value) for idx,value in enumerate(asgn))

    reversed_dict = {}
    for k, v in connections_dict.items():
        reversed_dict[v] = reversed_dict.get(v, []) + [k]
    # reversed_dict = defaultdict(list) #defaultdict(set) and replace append() with add() to remove duplicates
    # for key, value in connections_dict.items():
    #     reversed_dict[value].append(key)

    return reversed_dict


def assignments_to_dict(assignments_file):
    asgn = np.loadtxt(assignments_file, dtype =float).astype(int)
    asgn = [tuple(sorted(c)) for c in asgn]
    connections = list(conn for conn,_ in itertools.groupby(asgn))
    connections_dict = dict((idx, value) for idx,value in enumerate(asgn))

    reversed_dict = {}
    for k, v in connections_dict.items():
        reversed_dict[v] = reversed_dict.get(v, []) + [k]
    return reversed_dict