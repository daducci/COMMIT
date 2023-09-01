import torch


def COMMIT_A__block(x, Y, t_f, t_v, t_l,xPtr):

    # intra-cellular compartments
    for i in range(len(t_v)):
        x0 = x[t_f[i]]
        if x0 != 0:
            Y[t_v[i]] += t_l[i] * x0

    # isotropic compartments
    for i in range(len(t_v)):
        x0 = xPtr[i]
        if x0 != 0:
            Y[t_v[i]] += x0
    return Y


class LinearOperatorGPU:

    def __init__( self, DICTIONARY, KERNELS ) :
        # unpack the dictionary and kernels
        self.ICf        = torch.from_numpy(DICTIONARY['IC']['fiber']).cuda()
        self.ICv        = torch.from_numpy(DICTIONARY['IC']['v']).cuda()
        self.ICl        = torch.from_numpy(DICTIONARY['IC']['len']).cuda()
        self.ICo        = torch.from_numpy(DICTIONARY['IC']['o']).cuda()
        self.ECv        = torch.from_numpy(DICTIONARY['EC']['v']).cuda()
        self.ECo        = torch.from_numpy(DICTIONARY['EC']['o']).cuda()
        self.ISOv       = torch.from_numpy(DICTIONARY['ISO']['v']).cuda()

        self.wmrSFP     = torch.from_numpy(KERNELS['wmr']).cuda()
        self.wmhSFP     = torch.from_numpy(KERNELS['wmh']).cuda()
        self.isoSFP     = torch.from_numpy(KERNELS['iso']).cuda()

        self.nF         = DICTIONARY['IC']['nF']    # number of FIBERS
        self.nR         = KERNELS['wmr'].shape[0]   # number of FIBER RADII
        self.nE         = DICTIONARY['EC']['nE']    # number of EC segments
        self.nT         = KERNELS['wmh'].shape[0]   # number of EC TORTUOSITY values
        self.nV         = DICTIONARY['nV']          # number of VOXELS
        self.nI         = KERNELS['iso'].shape[0]   # number of ISO contributions
        self.n          = DICTIONARY['IC']['n']     # numbner of IC segments
        self.ndirs      = KERNELS['wmr'].shape[1]   # number of directions

        if KERNELS['wmr'].size > 0 :
            self.nS = KERNELS['wmr'].shape[2]       # number of SAMPLES
        elif KERNELS['wmh'].size > 0 :
            self.nS = KERNELS['wmh'].shape[2]
        else :
            self.nS = KERNELS['wmr'].shape[1]

        self.adjoint    = 0                         # direct of inverse product
        self.n1 = self.nV*self.nS
        self.n2 = self.nR*self.nF + self.nT*self.nE + self.nI*self.nV

    @property
    def T( self ) :
        """Transpose of the explicit matrix."""
        C = LinearOperator( self.DICTIONARY, self.KERNELS )
        C.adjoint = 1 - C.adjoint
        return C

    @property
    def shape( self ) :
        """Size of the explicit matrix."""
        if not self.adjoint :
            return ( self.n1, self.n2 )
        else :
            return ( self.n2, self.n1 )
        

    def dot( self,  v_in):
        """
        Parameters
        ----------
        v_in : 1D numpy.array of double
            Input vector for the matrix-vector multiplication

        Returns
        -------
        v_out : 1D numpy.array of double
            Results of the multiplication
        """
        # Permit only matrix-vector multiplications
        if v_in.size != self.shape[1] :
            raise RuntimeError( "A.dot(): dimensions do not match" )

        if not self.adjoint :
            # compute A*x dividing A into sub matrices with size the maximum GPU memory availabe
            # get gpu memory available
            mem = torch.cuda.mem_get_info()
            # get the maximum size of the sub matrices
            max_size = mem[0]//8
            # get the number of sub matrices
            n_sub = self.n2//max_size

            # create output array
            v_out = torch.zeros(self.shape[0], dtype=torch.float64)

            

            # compute A*x dividing A into sub matrices with size the maximum GPU memory availabe
            for i in range(n_sub):
                # compute the start and end indices for the sub matrix
                start = i*max_size
                end = (i+1)*max_size
                # compute the sub matrices
                x_sub = x[start:end]
                Y_sub = Y[start:end]
                # compute the sub COMMIT_A__block
                Y_sub = COMMIT_A__block(x_sub, Y_sub, t_f, t_v, t_l, xPtr)
                # update the output array
                v_out[start:end] = Y_sub

            # compute the last sub matrix
            if self.n2%max_size != 0:
                # compute the start and end indices for the sub matrix
                start = n_sub*max_size
                end = self.n2
                # compute the sub matrices
                x_sub = x[start:end]
                Y_sub = Y[start:end]
                # compute the sub COMMIT_A__block
                Y_sub = COMMIT_A__block(x_sub, Y_sub, t_f, t_v, t_l, xPtr)
                # update the output array
                v_out[start:end] = Y_sub
        else :
            # compute A^T*x dividing A into sub matrices with size the maximum GPU memory availabe
            # get gpu memory available
            mem = torch.cuda.mem_get_info()
            # get the maximum size of the sub matrices
            max_size = mem[0]//8
            # get the number of sub matrices
            n_sub = self.n1//max_size

            # create output array
            v_out = torch.zeros(self.shape[1], dtype=torch.float64)


            # compute A*x dividing A into sub matrices with size the maximum GPU memory availabe
            for i in range(n_sub):
                # compute the start and end indices for the sub matrix
                start = i*max_size
                end = (i+1)*max_size
                # compute the sub matrices
                x_sub = x[start:end]
                Y_sub = Y[start:end]
                # compute the sub COMMIT_A__block
                Y_sub = COMMIT_A__block(x_sub, Y_sub, t_f, t_v, t_l, xPtr)
                # update the output array
                v_out[start:end] = Y_sub

            # compute the last sub matrix
            if self.n1%max_size != 0:
                # compute the start and end indices for the sub matrix
                start = n_sub*max_size
                end = self.n1
                # compute the sub matrices
                x_sub = x[start:end]
                Y_sub = Y[start:end]
                # compute the sub COMMIT_A__block
                Y_sub = COMMIT_A__block(x_sub, Y_sub, t_f, t_v, t_l, xPtr)
                # update the output array
                v_out[start:end] = Y_sub

        return v_out


