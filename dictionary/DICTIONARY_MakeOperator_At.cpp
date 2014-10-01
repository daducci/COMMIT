// the following 2 lines are to fix a bug with XCode 5.1 in OSX
#include <stdint.h>
typedef uint16_t char16_t;

#include <stdint.h>
typedef uint16_t char16_t;

#include <cmath>
#include <matrix.h>
#include <mex.h>

// avoid to perform parameter checking each time
#define DO_CHECK 0

// number of intra-axonal compartments (different RADII)
#ifdef nIC
	#if (nIC<0 || nIC>4)
	#error nIC must be >= 0 and <= 4
	#endif	
#else
	#error "nIC" parameter must be passed to the compiler as "-DnIC=<value>"
#endif

// number of extra-axonal compartments (different TORTUOSITY)
#ifdef nEC
	#if (nEC<0 || nEC>4)
	#error nEC must be >= 0 and <= 4
	#endif	
#else
	#error "nEC" parameter must be passed to the compiler as "-DnEC=<value>"
#endif

// number of isotropic compartments (CSF and others)
#ifdef nISO
	#if (nISO<0 || nISO>4)
	#error nISO must be >= 0 and <= 4
	#endif	
#else
	#error "nISO" parameter must be passed to the compiler as "-DnISO=<value>"
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray* tmp;

	/* ================ */
	/* Check the INPUTS */
	/* ================ */
	#if DO_CHECK > 0
 	if( nrhs != 3 )
		mexErrMsgIdAndTxt("InvalidInput:nrhs", "Require 3 inputs.");
	if ( !mxIsStruct(prhs[0]) )
		mexErrMsgIdAndTxt("InvalidInput:DICTIONARY", "'DICTIONARY' must be a struct");
	#endif

	mxArray* IC  = mxGetField( prhs[0], 0, "IC" );
	mxArray* EC  = mxGetField( prhs[0], 0, "EC" );
	mxArray* ISO = mxGetField( prhs[0], 0, "ISO" );
	#if DO_CHECK > 0
	if ( !mxIsStruct(IC) || !mxIsStruct(EC) || !mxIsStruct(ISO) )
		mexErrMsgIdAndTxt("InvalidInput:DICTIONARY", "'DICTIONARY' format is wrong");
	#endif

	// Parse "n"
	tmp = mxGetField( IC, 0, "n" );
	#if DO_CHECK > 0
	if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
		mexErrMsgIdAndTxt("InvalidInput:n","'n' must be a real scalar");
	#endif
	const int n = mxGetScalar( tmp );

	// Parse "nF"
	tmp = mxGetField( IC, 0, "nF" );
	#if DO_CHECK > 0
	if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
		mexErrMsgIdAndTxt("InvalidInput:IC.nF","'IC.nF' must be a real scalar");
	#endif
	const int nF = mxGetScalar( tmp );

	// Parse "nS"
	tmp = mxGetField( prhs[1], 0, "nS" );
	#if DO_CHECK > 0
	if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
		mexErrMsgIdAndTxt("InvalidInput:nS","'nS' must be a real scalar");
	#endif
	const int nS = mxGetScalar( tmp );

	// Parse "nV"
	tmp = mxGetField( prhs[0], 0, "nV" );
	#if DO_CHECK > 0
	if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
		mexErrMsgIdAndTxt("InvalidInput:nV","'nV' must be a real scalar");
	#endif
	const int nV = mxGetScalar( tmp );

	// Parse "nE"
	tmp = mxGetField( EC, 0, "nE" );
	#if DO_CHECK > 0
	if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
		mexErrMsgIdAndTxt("InvalidInput:EC.nE","'EC.nE' must be a real scalar");
	#endif
	const int nE = mxGetScalar( tmp );

	// Parse "fiber"
	tmp = mxGetField( IC, 0, "fiber" );
	#if DO_CHECK > 0
	if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 )
		mexErrMsgIdAndTxt("InvalidInput:fiber","'fiber' must be a n*1 vector");
	#endif
 	UINT32_T* fiber = (UINT32_T*) mxGetData( tmp );

	// Parse "ICv", "ICo"
	tmp = mxGetField( IC, 0, "v" );
	#if DO_CHECK > 0
	if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 )
		mexErrMsgIdAndTxt("InvalidInput:IC.v","'IC.v' must be a n*1 vector");
	#endif
	UINT32_T* ICv = (UINT32_T*) mxGetData(tmp);
	tmp = mxGetField( IC, 0, "o" );
	#if DO_CHECK > 0
	if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 )
		mexErrMsgIdAndTxt("InvalidInput:IC.o","'IC.o' must be a n*1 vector");
	#endif
	UINT16_T* ICo = (UINT16_T*) mxGetData(tmp);

	tmp = mxGetField( IC, 0, "len" );
	#if DO_CHECK > 0
	if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 )
		mexErrMsgIdAndTxt("InvalidInput:len","'len' must be a n*1 vector");
	#endif
 	float* len = (float*) mxGetData(tmp);

	// Parse "ECv","ECo"
	tmp = mxGetField( EC, 0, "v" );
	#if DO_CHECK > 0
	if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 )
		mexErrMsgIdAndTxt("InvalidInput:EC.v","'EC.v' must be a n*1 vector");
	#endif
	UINT32_T* ECv = (UINT32_T*) mxGetData(tmp);
	tmp = mxGetField( EC, 0, "o" );
	#if DO_CHECK > 0
	if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 )
		mexErrMsgIdAndTxt("InvalidInput:EC.o","'EC.x' must be a n*1 vector");
	#endif
	UINT16_T* ECo = (UINT16_T*) mxGetData(tmp);

	// Parse "ISOv"
	tmp = mxGetField( ISO, 0, "v" );
	#if DO_CHECK > 0
	if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 )
		mexErrMsgIdAndTxt("InvalidInput:ISO.v","'ISO.v' must be a n*1 vector");
	#endif
	UINT32_T* ISOv = (UINT32_T*) mxGetData(tmp);
	
	// Parse "KERNELS.wmr", "KERNELS.wmh" and "KERNELS.iso"
	mxArray* wmr = mxGetField( prhs[1], 0, "wmr" );
	#if DO_CHECK > 0
	if ( !mxIsCell(wmr) )
		mexErrMsgIdAndTxt("InvalidInput:wmr","'wmr' must be a cell array");
	#endif
	#if nIC>=1
	double* wmrSFP0 = (double*) mxGetData( mxGetCell(wmr,0) );
	#endif
	#if nIC>=2
	double* wmrSFP1 = (double*) mxGetData( mxGetCell(wmr,1) );
	#endif
	#if nIC>=3
	double* wmrSFP2 = (double*) mxGetData( mxGetCell(wmr,2) );
	#endif
	#if nIC>=4
	double* wmrSFP3 = (double*) mxGetData( mxGetCell(wmr,3) );
	#endif
	
	#if nEC>=1
	mxArray* wmh = mxGetField( prhs[1], 0, "wmh" );
	#if DO_CHECK > 0
	if ( !mxIsCell(wmh) )
		mexErrMsgIdAndTxt("InvalidInput:wmh","'wmh' must be a cell array");
	#endif
	double* wmhSFP0 = (double*) mxGetData( mxGetCell(wmh,0) );
	#if nEC>=2
	double* wmhSFP1 = (double*) mxGetData( mxGetCell(wmh,1) );
	#endif
	#if nEC>=3
	double* wmhSFP2 = (double*) mxGetData( mxGetCell(wmh,2) );
	#endif
	#if nEC>=4
	double* wmhSFP3 = (double*) mxGetData( mxGetCell(wmh,3) );
	#endif
	#endif

	#if nISO>=1
	mxArray* iso = mxGetField( prhs[1], 0, "iso" );
	#if DO_CHECK > 0
	if ( !mxIsCell(iso) )
		mexErrMsgIdAndTxt("InvalidInput:iso","'iso' must be a cell array");
	#endif
	double* isoSFP0 = (double*) mxGetData( mxGetCell(iso,0) );
	#if nISO>=2
	double* isoSFP1 = (double*) mxGetData( mxGetCell(iso,1) );
	#endif
	#if nISO>=3
	double* isoSFP2 = (double*) mxGetData( mxGetCell(iso,2) );
	#endif
	#if nISO>=4
	double* isoSFP3 = (double*) mxGetData( mxGetCell(iso,3) );
	#endif
	#endif

	// Parse "Y"
	#if DO_CHECK > 0
	if ( mxGetNumberOfDimensions( prhs[2] ) != 4 )
		mexErrMsgIdAndTxt("InvalidInput:Y","'Y' must be a 4D matrix");
	#endif
 	double* Y = (double*) mxGetData( prhs[2] );

	const mwSize* Y_dims = mxGetDimensions( prhs[2] );
	const int Y_dimx = (int)Y_dims[1];
	const int Y_dimy = (int)Y_dims[2];
	const int Y_dimz = (int)Y_dims[3];


	/* =============== */
	/* Set the OUTPUTS */
	/* =============== */
	#if DO_CHECK > 0
	if( nlhs != 1 )
		mexErrMsgIdAndTxt("InvalidOutput:nlhs", "Required 1 output.");
	#endif

	const int outDims[4] = { nIC*nF + nEC*nE + nISO*nV, 1 };
    plhs[0] = mxCreateNumericArray(2, outDims, mxDOUBLE_CLASS, mxREAL);
    double *x = (double*)mxGetData( plhs[0] );


	/* ========= */
	/* Main loop */
	/* ========= */
	double x0_tmp, x1_tmp, x2_tmp, x3_tmp, Y_tmp;
	double *Yptr, *YptrEnd;
	int offset;
	

	/* intra-cellular compartments */
	#if nIC>=1

	double *wmrSFP0ptr;
	#if nIC>=2
	double *wmrSFP1ptr;
	#endif
	#if nIC>=3
	double *wmrSFP2ptr;
	#endif
	#if nIC>=4
	double *wmrSFP3ptr;
	#endif

	UINT32_T* ICvEnd = ICv+n;
	while( ICv != ICvEnd )
	{
		Yptr         = Y         + nS * (*ICv++);
		YptrEnd      = Yptr      + nS;

		offset = nS * (*ICo++);
		Y_tmp = *Yptr;
		wmrSFP0ptr   = wmrSFP0   + offset;
		x0_tmp = (*wmrSFP0ptr++) * Y_tmp;
		#if nIC>=2
		wmrSFP1ptr   = wmrSFP1   + offset;
		x1_tmp = (*wmrSFP1ptr++) * Y_tmp;
		#endif
		#if nIC>=3
		wmrSFP2ptr   = wmrSFP2   + offset;
		x2_tmp = (*wmrSFP2ptr++) * Y_tmp;
		#endif
		#if nIC>=4
		wmrSFP3ptr   = wmrSFP3   + offset;
		x3_tmp = (*wmrSFP3ptr++) * Y_tmp;
		#endif
		
		double w;
		while( ++Yptr != YptrEnd )
		{
			Y_tmp = *Yptr;
			x0_tmp += (*wmrSFP0ptr++) * Y_tmp;
			#if nIC>=2
			x1_tmp += (*wmrSFP1ptr++) * Y_tmp;
			#endif
			#if nIC>=3
			x2_tmp += (*wmrSFP2ptr++) * Y_tmp;
			#endif
			#if nIC>=4
			x3_tmp += (*wmrSFP3ptr++) * Y_tmp;
			#endif
		}

		w = (double)(*len++);
		x[*fiber]      += w * x0_tmp;
		#if nIC>=2
		x[*fiber+nF]   += w * x1_tmp;
		#endif
		#if nIC>=3
		x[*fiber+2*nF] += w * x2_tmp;
		#endif
		#if nIC>=4
		x[*fiber+3*nF] += w * x3_tmp;
		#endif

		fiber++;
	}
	#endif

	
	/* extra-cellular compartments */
	#if nEC>=1

	double *wmhSFP0ptr, *x_wmhPtr0 = x + nIC*nF, *x_wmhPtr0End = x_wmhPtr0 + nE;
	#if nEC>=2
	double *wmhSFP1ptr, *x_wmhPtr1 = x_wmhPtr0 + nE;
	#endif
	#if nEC>=3
	double *wmhSFP2ptr, *x_wmhPtr2 = x_wmhPtr1 + nE;
	#endif
	#if nEC>=4
	double *wmhSFP3ptr, *x_wmhPtr3 = x_wmhPtr2 + nE;
	#endif

	while( x_wmhPtr0 != x_wmhPtr0End )
	{
		Yptr         = Y         + nS * (*ECv++);
		YptrEnd      = Yptr      + nS;

		offset = nS * (*ECo++);
		
		Y_tmp = *Yptr;
		wmhSFP0ptr    = wmhSFP0  + offset;
		x0_tmp = (*wmhSFP0ptr++) * Y_tmp;
		#if nEC>=2
		wmhSFP1ptr    = wmhSFP1  + offset;
		x1_tmp = (*wmhSFP1ptr++) * Y_tmp;
		#endif
		#if nEC>=3
		wmhSFP2ptr    = wmhSFP2  + offset;
		x2_tmp = (*wmhSFP2ptr++) * Y_tmp;
		#endif
		#if nEC>=4
		wmhSFP3ptr    = wmhSFP3  + offset;
		x3_tmp = (*wmhSFP3ptr++) * Y_tmp;
		#endif

		while( ++Yptr != YptrEnd )
		{
			Y_tmp = *Yptr;
			x0_tmp += (*wmhSFP0ptr++) * Y_tmp;
			#if nEC>=2
			x1_tmp += (*wmhSFP1ptr++) * Y_tmp;
			#endif
			#if nEC>=3
			x2_tmp += (*wmhSFP2ptr++) * Y_tmp;
			#endif
			#if nEC>=4
			x3_tmp += (*wmhSFP3ptr++) * Y_tmp;
			#endif
		}
		(*x_wmhPtr0++) += x0_tmp;
		#if nEC>=2
		(*x_wmhPtr1++) += x1_tmp;
		#endif
		#if nEC>=3
		(*x_wmhPtr2++) += x2_tmp;
		#endif
		#if nEC>=4
		(*x_wmhPtr3++) += x3_tmp;
		#endif
	}
	#endif


	/* isotropic compartments */
	#if nISO>=1

	double *isoSFP0ptr, *x_isoPtr0 = x + nIC*nF + nEC*nE, *x_isoPtr0End = x_isoPtr0 + nV;
	#if nISO>=2
	double *isoSFP1ptr, *x_isoPtr1 = x_isoPtr0 + nV;
	#endif
	#if nISO>=3
	double *isoSFP2ptr, *x_isoPtr2 = x_isoPtr1 + nV;
	#endif
	#if nISO>=4
	double *isoSFP3ptr, *x_isoPtr3 = x_isoPtr2 + nV;
	#endif

	while( x_isoPtr0 != x_isoPtr0End )
	{
		Yptr         = Y         + nS * (*ISOv++);
		YptrEnd      = Yptr      + nS;

		isoSFP0ptr   = isoSFP0;
		#if nISO>=2
		isoSFP1ptr   = isoSFP1;
		#endif
		#if nISO>=3
		isoSFP2ptr   = isoSFP2;
		#endif
		#if nISO>=4
		isoSFP3ptr   = isoSFP3;
		#endif
		
		Y_tmp = *Yptr;
		x0_tmp = (*isoSFP0ptr++) * Y_tmp;
		#if nISO>=2
		x1_tmp = (*isoSFP1ptr++) * Y_tmp;
		#endif
		#if nISO>=3
		x2_tmp = (*isoSFP2ptr++) * Y_tmp;
		#endif
		#if nISO>=4
		x3_tmp = (*isoSFP3ptr++) * Y_tmp;
		#endif
		while( ++Yptr != YptrEnd )
		{
			Y_tmp = *Yptr;
			x0_tmp  += (*isoSFP0ptr++) * Y_tmp;
			#if nISO>=2
			x1_tmp  += (*isoSFP1ptr++) * Y_tmp;
			#endif
			#if nISO>=3
			x2_tmp  += (*isoSFP2ptr++) * Y_tmp;
			#endif
			#if nISO>=4
			x3_tmp  += (*isoSFP3ptr++) * Y_tmp;
			#endif
		}
		(*x_isoPtr0++) += x0_tmp;
		#if nISO>=2
		(*x_isoPtr1++) += x1_tmp;
		#endif
		#if nISO>=3
		(*x_isoPtr2++) += x2_tmp;
		#endif
		#if nISO>=4
		(*x_isoPtr3++) += x3_tmp;
		#endif
	}
	#endif

    return;
}
