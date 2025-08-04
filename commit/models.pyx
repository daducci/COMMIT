#!python
#cython: language_level=3, boundscheck=False, wraparound=False, profile=False
from os.path import join as pjoin
import numpy as np
import nibabel
from amico.models import BaseModel, StickZeppelinBall as _StickZeppelinBall, CylinderZeppelinBall as _CylinderZeppelinBall
import amico.util as util
from dicelib.ui import setup_logger, set_verbose, ProgressBar
logger = setup_logger('models')


class StickZeppelinBall(_StickZeppelinBall):
    """Simulate the response functions according to the Stick-Zeppelin-Ball model.
    See the AMICO.model module for details.
    """
    def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
        #FIXME: use logger
        util.set_verbose(2)
        return super().resample(in_path, idx_out, Ylm_out, doMergeB0, ndirs)


class CylinderZeppelinBall(_CylinderZeppelinBall):
    """Simulate the response functions according to the Cylinder-Zeppelin-Ball model.
    See the AMICO.model module for details.
    """
    def resample(self, in_path, idx_out, Ylm_out, doMergeB0, ndirs):
        #FIXME: use logger
        util.set_verbose(2)
        return super().resample(in_path, idx_out, Ylm_out, doMergeB0, ndirs)


class ScalarMap( BaseModel ) :
    """Implements a simple model where each compartment contributes to a scalar map,
    e.g. intra-axonsl signal fraction or myelin water fraction, proportionally to
    its local length inside each voxel.
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


    def _postprocess(self, evaluation, xic):
        """Rescale the streamline weights using the local tissue damage estimated
        in all imaging voxels with the COMMIT_lesion new model.

        Parameters
        ----------
        evaluation : object
            Evaluation object, to enable accessing the whole content of the
            whole evaluation object

        xic : np.array
            Streamline weights

        Returns
        -------
        np.array
            The rescaled streamline weights accounting for lesions
        """
        if not self.lesion_mask:
            # nothing to do if lesion mask is not given
            return xic

        RESULTS_path = evaluation.get_config('RESULTS_path')
        niiISO_img = np.asanyarray( nibabel.load( pjoin(RESULTS_path,'compartment_ISO.nii.gz') ).dataobj ).astype(np.float32)
        ISO = niiISO_img[evaluation.DICTIONARY['MASK_ix'], evaluation.DICTIONARY['MASK_iy'], evaluation.DICTIONARY['MASK_iz']]
        if np.count_nonzero(ISO>0) == 0:
            logger.warning('No lesions found')
            return xic

        # rescale the input scalar map in each voxel according to estimated lesion contributions
        niiIC_img = np.asanyarray( nibabel.load( pjoin(RESULTS_path,'compartment_IC.nii.gz') ).dataobj ).astype(np.float32)
        IC = niiIC_img[evaluation.DICTIONARY['MASK_ix'], evaluation.DICTIONARY['MASK_iy'], evaluation.DICTIONARY['MASK_iz']]
        ISO_scaled = np.zeros_like(ISO, dtype=np.float32)
        ISO_scaled[ISO>0] = (IC[ISO>0] - ISO[ISO>0]) / IC[ISO>0]
        ISO_scaled_save = np.zeros_like(niiISO_img, dtype=np.float32)
        ISO_scaled_save[evaluation.DICTIONARY['MASK_ix'], evaluation.DICTIONARY['MASK_iy'], evaluation.DICTIONARY['MASK_iz']] = ISO_scaled
        affine = evaluation.niiDWI.affine if nibabel.__version__ >= '2.0.0' else evaluation.niiDWI.get_affine()
        nibabel.save(nibabel.Nifti1Image(ISO_scaled_save, affine), pjoin(RESULTS_path,'compartment_IC_lesion_scaled.nii.gz'))

        # save the map of local tissue damage estimated in each voxel
        nibabel.save( nibabel.Nifti1Image( niiISO_img, affine ), pjoin(RESULTS_path,'compartment_lesion.nii.gz') )

        # override ISO map and set it to 0
        nibabel.save( nibabel.Nifti1Image( 0*niiISO_img, affine), pjoin(RESULTS_path,'compartment_ISO.nii.gz') )

        # rescale each streamline weight
        kept = evaluation.DICTIONARY['TRK']['kept']
        cdef double [::1] xic_view = xic[kept==1]
        cdef double [::1] xic_scaled_view = xic[kept==1].copy()
        cdef float [::1] ISO_scaled_view = ISO_scaled
        cdef unsigned int [::1] idx_v_view = evaluation.DICTIONARY['IC']['v']
        cdef unsigned int [::1] idx_f_view = evaluation.DICTIONARY['IC']['fiber']
        cdef size_t i, idx_v, idx_f
        cdef double val

        # Rescaling streamline weights accounting for lesions
        for i in range(evaluation.DICTIONARY['IC']['v'].shape[0]):
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