#!python
#cython: language_level=3, boundscheck=False, wraparound=False, profile=False

from os import cpu_count as num_cpu
from os.path import join as pjoin

from setuptools import Extension

from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import nibabel as nib

from amico.models import BaseModel, StickZeppelinBall as _StickZeppelinBall, CylinderZeppelinBall as _CylinderZeppelinBall
import amico.util as util

from dicelib.ui import setup_logger, set_verbose, ProgressBar


logger = setup_logger('models')



class StickZeppelinBall(_StickZeppelinBall):
    """Simulate the response functions according to the Stick-Zeppelin-Ball model.
    See the AMICO.model module for details.
    """
    def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
        util.set_verbose(2)
        return super().resample(in_path, idx_out, Ylm_out, doMergeB0, ndirs)


class CylinderZeppelinBall(_CylinderZeppelinBall):
    """Simulate the response functions according to the Cylinder-Zeppelin-Ball model.
    See the AMICO.model module for details.
    """
    def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
        util.set_verbose(2)
        return super().resample(in_path, idx_out, Ylm_out, doMergeB0, ndirs)


class VolumeFractions(BaseModel):
    """Implements a simple model where each compartment contributes only with
       its own volume fraction. This model has been created to test there
       ability to remove false positive fibers with COMMIT.
    """
    def __init__(self):
        self.id = 'VolumeFractions'
        self.name = 'Volume fractions'
        self.maps_name = []
        self.maps_descr = []
        self.nolut = True

    def set(self):
        return

    def get_params(self):
        params = {}
        params['id'] = self.id
        params['name'] = self.name
        return params

    def set_solver(self):
        logger.error('Not implemented')

    def generate(self, out_path, aux, idx_in, idx_out, ndirs):
        return

    def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
        if doMergeB0:
            nS = 1 + self.scheme.dwi_count
            merge_idx = np.hstack((self.scheme.b0_idx[0], self.scheme.dwi_idx))
        else:
            nS = self.scheme.nS
            merge_idx = np.arange(nS)

        KERNELS = {}
        KERNELS['model'] = self.id
        KERNELS['wmr'] = np.ones((1, ndirs, nS), dtype=np.float32)
        KERNELS['wmh'] = np.ones((0, ndirs, nS), dtype=np.float32)
        KERNELS['iso'] = np.ones((0, nS), dtype=np.float32)
        return KERNELS

    def fit(self, evaluation):
        logger.error('Not implemented')


class VolumeFractions_w_Lesions( BaseModel ) :
    """Implements a simple model where each compartment contributes only with
       its own volume fraction. This model has been created to test there
       ability to remove false positive fibers with COMMIT.
    """

    def __init__( self ) :
        self.id             = 'Lesion'
        self.name           = 'Lesion'
        self.restrictedISO  = None
        self.maps_name      = [ ]
        self.maps_descr     = [ ]

    def get_params( self ) :
        params = {}
        params['id'] = self.id
        params['name'] = self.name
        params["ISO"] = self.restrictedISO
        return params

    def set( self, ISO_map=None ) :
        self.restrictedISO = ISO_map

    def set_solver( self ) :
        logger.error( 'Not implemented' )

    def generate( self, out_path, aux, idx_in, idx_out, ndirs ) :
        return
    
    def fit(self, evaluation):
        """Placeholder implementation for the abstract method."""
        logger.error('Not implemented')

    def resample( self, in_path, idx_out, Ylm_out, doMergeB0, ndirs ) :
        if doMergeB0:
            nS = 1+self.scheme.dwi_count
            merge_idx = np.hstack((self.scheme.b0_idx[0],self.scheme.dwi_idx))
        else:
            nS = self.scheme.nS
            merge_idx = np.arange(nS)

        KERNELS = {}
        KERNELS['model'] = self.id
        KERNELS['wmr']   = np.ones( (1,ndirs,nS), dtype=np.float32 )
        KERNELS['wmh']   = np.ones( (0,ndirs,nS), dtype=np.float32 )
        if self.restrictedISO:
            KERNELS['iso']   = -1 * np.ones( (1,nS), dtype=np.float32 )
        else:
            KERNELS['iso']   = np.ones( (0,nS), dtype=np.float32 )

        return KERNELS
    
    def find_idx(self, ic_v, ic_les, dict_icf, max_size, progress_bar, chunk_num):
        cdef unsigned int [::1]ic_les_view = ic_les
        cdef long les_shape = ic_les.shape[0]
        cdef long ic_v_shape = ic_v.shape[0]
        cdef int chunk_num_mem = chunk_num
        cdef unsigned int [::1]progress_bar_mem = progress_bar


        # indices_vox = np.zeros(les_shape*ic_v_shape, dtype=np.uint32)
        indices_vox = np.zeros(max_size, dtype=np.uint32)
        cdef unsigned int [::1]indices_vox_view = indices_vox

        # indices_str = np.zeros(les_shape*ic_v_shape, dtype=np.uint32)
        indices_str = np.zeros(max_size, dtype=np.uint32)
        cdef unsigned int [::1]indices_str_view = indices_str

        cdef unsigned int[::1] ic_v_view = ic_v

        cdef int count = 0

        cdef unsigned int [::1]dict_icf_view = dict_icf

        cdef unsigned int ic_les_val = 0
        cdef size_t i, j = 0

        with nogil:
            for i in range(les_shape):
                ic_les_val = ic_les_view[i]
                for j in range(ic_v_shape):
                    if ic_v_view[j] == ic_les_val:# and ic_les_val > 0:
                        indices_vox_view[count] = ic_v_view[j]
                        indices_str_view[count] = dict_icf_view[j]
                        count += 1
                progress_bar_mem[chunk_num_mem] += 1

        return indices_vox[:count], indices_str[:count]
    
    def _postprocess(self, preproc_dict, verbose=1):
        set_verbose("commitwipmodels", verbose)

        ISO_map = np.array(preproc_dict['niiISO_img'], dtype=np.float32)
        IC_map = np.array(preproc_dict['niiIC_img'], dtype=np.float32)

        dictionary = preproc_dict['DICTIONARY']

        kept = dictionary['TRK']['kept']
        x_weights = preproc_dict["streamline_weights"][kept==1].copy()
        x_weights_scaled = x_weights.copy()
        new_weights = np.zeros_like(kept, dtype=np.float32)


        IC = IC_map[dictionary['MASK_ix'], dictionary['MASK_iy'], dictionary['MASK_iz']]
        ISO = ISO_map[dictionary['MASK_ix'], dictionary['MASK_iy'], dictionary['MASK_iz']]
        ISO_scaled = np.zeros_like(ISO) 
        ISO_scaled[ISO>0] = (IC[ISO>0] - ISO[ISO>0]) / IC[ISO>0]

        ISO_scaled_save = np.zeros_like(ISO_map)
        ISO_scaled_save[dictionary['MASK_ix'], dictionary['MASK_iy'], dictionary['MASK_iz']] = ISO_scaled
        nib.save(nib.Nifti1Image(ISO_scaled_save, preproc_dict['affine']), pjoin(preproc_dict["RESULTS_path"],'IC_lesion_scaled.nii.gz'))
        
        idx_les = np.argwhere(ISO_scaled > 0)[:,0].astype(np.uint32)

    
        result = []
        dict_idx_v = dictionary['IC']['v']
        cdef unsigned int [::1]dict_idx_v_view = dict_idx_v

        cdef unsigned int cpu_count = num_cpu()
        cdef unsigned int[:] find_idx_progress = np.zeros(cpu_count, dtype=np.uint32)

        n = idx_les.shape[0]
        c = n // cpu_count
        max_size = int(3e9/cpu_count)
        chunks = []
        for ii, jj in zip(range(0, n, c), range(c, n+1, c)):
            chunks.append((ii, jj))
        if chunks[len(chunks)-1][1] != n:
            chunks[len(chunks)-1] = (chunks[len(chunks)-1][0], n)

        logger.subinfo('Recomputing streamlines weights accounting for lesions', indent_lvl=2, indent_char='-', with_progress=True)
        with ProgressBar(multithread_progress=find_idx_progress, total=n,
                         disable=verbose<3, hide_on_exit=True, subinfo=True) as pbar:
            with ThreadPoolExecutor(max_workers=cpu_count) as executor:
                futures = [executor.submit(self.find_idx, dict_idx_v_view, idx_les[ii:jj], dictionary['IC']['fiber'], max_size, find_idx_progress, chunk_num) for chunk_num, (ii, jj) in enumerate(chunks)]
                for future in as_completed(futures):
                    result.append(future.result())

        idx_vox = []
        idx_str = []
        for r in result:
            idx_vox.extend(r[0].tolist())
            idx_str.extend(r[1].tolist())

        cdef double [::1] x_weights_view = x_weights
        cdef double [::1] x_weights_scaled_view = x_weights_scaled
        cdef double x_weight = 0
        cdef float [::1] ISO_scaled_view = ISO_scaled

        cdef size_t vox = 0
        cdef size_t str_i = 0

        for vox, str_i in zip(idx_vox, idx_str):
            x_weight = x_weights_view[str_i] * ISO_scaled_view[vox]
            if x_weight < x_weights_scaled_view[str_i]:
                x_weights_scaled_view[str_i] = x_weight

        new_weights[kept==1] = x_weights_scaled
        coeffs_format='%.5e'

        np.savetxt( pjoin(preproc_dict["RESULTS_path"],'streamline_weights_reweighted.txt'), new_weights, fmt=coeffs_format )