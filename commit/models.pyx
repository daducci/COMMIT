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


# class VolumeFractions(BaseModel):
#     """Implements a simple model where each compartment contributes only with
#        its own volume fraction. This model has been created to test there
#        ability to remove false positive fibers with COMMIT.
#     """
#     def __init__(self):
#         self.id = 'VolumeFractions'
#         self.name = 'Volume fractions'
#         self.maps_name = []
#         self.maps_descr = []
#         self.nolut = True

#     def set(self):
#         return

#     def get_params(self):
#         params = {}
#         params['id'] = self.id
#         params['name'] = self.name
#         return params

#     def set_solver(self):
#         logger.error('Not implemented')

#     def generate(self, out_path, aux, idx_in, idx_out, ndirs):
#         return

#     def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
#         if doMergeB0:
#             nS = 1 + self.scheme.dwi_count
#             merge_idx = np.hstack((self.scheme.b0_idx[0], self.scheme.dwi_idx))
#         else:
#             nS = self.scheme.nS
#             merge_idx = np.arange(nS)

#         KERNELS = {}
#         KERNELS['model'] = self.id
#         KERNELS['wmr'] = np.ones((1, ndirs, nS), dtype=np.float32)
#         KERNELS['wmh'] = np.ones((0, ndirs, nS), dtype=np.float32)
#         KERNELS['iso'] = np.ones((0, nS), dtype=np.float32)
#         return KERNELS

#     def fit(self, evaluation):
#         logger.error('Not implemented')


class ScalarMap( BaseModel ) :
    """Implements a simple model where each compartment contributes only with
       its own volume fraction. This model has been created to test there
       ability to remove false positive fibers with COMMIT.
    """

    def __init__( self ) :
        self.id             = 'ScalarMap'
        self.name           = 'ScalarMap'
        self.lesion_mask    = False
        self.maps_name      = [ ]
        self.maps_descr     = [ ]

    def get_params( self ) :
        params = {}
        params['id'] = self.id
        params['name'] = self.name
        return params

    def set( self ) :
        return

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
        if self.lesion_mask:
            KERNELS['iso']   = -1 * np.ones( (1,nS), dtype=np.float32 )
        else:
            KERNELS['iso']   = np.ones( (0,nS), dtype=np.float32 )

        return KERNELS


    def _postprocess(self, temp_data, verbose=1):
        """Rescale the streamline weights using the estimate tissue damage in all imaging voxels"""
        dictionary = temp_data['DICTIONARY']

        # compute the actual damage in each voxel
        IC_map = np.array(temp_data['niiIC_img'], dtype=np.float32)
        IC = IC_map[dictionary['MASK_ix'], dictionary['MASK_iy'], dictionary['MASK_iz']]
        ISO_map = np.array(temp_data['niiISO_img'], dtype=np.float32)
        ISO = ISO_map[dictionary['MASK_ix'], dictionary['MASK_iy'], dictionary['MASK_iz']]
        ISO_scaled = np.zeros_like(ISO)
        ISO_scaled[ISO>0] = (IC[ISO>0] - ISO[ISO>0]) / IC[ISO>0]
        ISO_scaled_save = np.zeros_like(ISO_map)
        ISO_scaled_save[dictionary['MASK_ix'], dictionary['MASK_iy'], dictionary['MASK_iz']] = ISO_scaled
        nib.save(nib.Nifti1Image(ISO_scaled_save, temp_data['affine']), pjoin(temp_data['RESULTS_path'],'compartment_IC_lesion_scaled.nii.gz'))
        del IC, IC_map, ISO, ISO_map, ISO_scaled_save
        if np.count_nonzero(ISO_scaled>0) == 0:
            logger.error('No lesions found in the input image.')
            return

        # rescale each streamline weight
        kept = dictionary['TRK']['kept']
        cdef double [::1] xic_view = temp_data['streamline_weights'][kept==1]
        cdef double [::1] xic_scaled_view = temp_data['streamline_weights'][kept==1].copy()
        cdef float [::1] ISO_scaled_view = ISO_scaled
        cdef unsigned int [::1] idx_v_view = dictionary['IC']['v']
        cdef unsigned int [::1] idx_f_view = dictionary['IC']['fiber']
        cdef size_t i, idx_v, idx_f
        cdef double val

        log_list = []
        ret_subinfo = logger.subinfo('Rescaling streamline weights accounting for lesions', indent_lvl=2, indent_char='-', with_progress=True)
        with ProgressBar(disable=False, hide_on_exit=True, subinfo=ret_subinfo, log_list=log_list):
            for i in range(dictionary['IC']['v'].shape[0]):
                idx_v = idx_v_view[i]
                val = ISO_scaled_view[idx_v]
                if val > 0:
                    idx_f = idx_f_view[i]
                    #TODO: allow considering other than the min value
                    if xic_view[idx_f] * val < xic_scaled_view[idx_f]:
                        xic_scaled_view[idx_f] = xic_view[idx_f] * val

        # return rescaled streamline weights
        xic_scaled = np.zeros_like(kept, dtype=np.float32)
        xic_scaled[kept==1] = xic_scaled_view
        return xic_scaled