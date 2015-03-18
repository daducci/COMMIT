// NOTE: uncomment the following 2 lines are to fix a bug with XCode 5.1 in OSX
//#include <stdint.h>
//typedef uint16_t char16_t;

#include <pthread.h>
#include <cmath>
#include <matrix.h>
#include <mex.h>

// avoid to perform parameter checking each time
#define DO_CHECK 1

// number of THREADS
#ifdef nTHREADS
    #if (nTHREADS<1 || nTHREADS>32)
     #error "nTHREADS" must be >= 1 and <= 32
    #endif
#else
    #error "nTHREADS" parameter must be passed to the compiler as "-DnTHREADS=<value>"
#endif

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


/* global variables */
int				n, nV, nE, nF, nS;
double			*x, *Y;
UINT8_T		    *ICthreads;
UINT32_T        *ECthreads, *ISOthreads;
UINT32_T		*ICf, *ICv, *ECv, *ISOv;
UINT16_T		*ICo, *ECo;
float			*ICl;

#if nIC>=1
float *wmrSFP0;
#endif
#if nIC>=2
float *wmrSFP1;
#endif
#if nIC>=3
float *wmrSFP2;
#endif
#if nIC>=4
float *wmrSFP3;
#endif

#if nEC>=1
float *wmhSFP0;
#endif
#if nEC>=2
float *wmhSFP1;
#endif
#if nEC>=3
float *wmhSFP2;
#endif
#if nEC>=4
float *wmhSFP3;
#endif

#if nISO>=1
float *isoSFP0;
#endif
#if nISO>=2
float *isoSFP1;
#endif
#if nISO>=3
float *isoSFP2;
#endif
#if nISO>=4
float *isoSFP3;
#endif


/* ================================================ */
/* Compute a sub-block of the MATRIX-VECTOR product */
/* ================================================ */
void* computeProductBlock( void *ptr )
{
    int      id = (long)ptr;
    int      offset;
     double   x0, x1, x2, x3, w, Y_tmp;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3;
    double   *Yptr, *YptrEnd;
    float   *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr;
    UINT32_T *t_v, *t_vEnd, *t_f;
    UINT16_T *t_o;
    float    *t_l;
    UINT8_T  *t_t;

    /* intra-cellular compartments */
#if nIC>=1
    t_v    = ICv;
    t_vEnd = ICv + n;
    t_o    = ICo;
    t_l    = ICl;
    t_f    = ICf;
    t_t    = ICthreads;

    while( t_v != t_vEnd )
    {
        // in this case, I need to walk throug because the segments are ordered in "voxel order"
        if ( *t_t == id )
        {
            Yptr    = Y    + (*t_v);
            YptrEnd = Yptr + nS;
            offset  = nS * (*t_o);

            Y_tmp = *Yptr;
            SFP0ptr   = wmrSFP0 + offset;
            x0 = (*SFP0ptr++) * Y_tmp;
            #if nIC>=2
            SFP1ptr   = wmrSFP1 + offset;
            x1 = (*SFP1ptr++) * Y_tmp;
            #endif
            #if nIC>=3
            SFP2ptr   = wmrSFP2 + offset;
            x2 = (*SFP2ptr++) * Y_tmp;
            #endif
            #if nIC>=4
            SFP3ptr   = wmrSFP3 + offset;
            x3 = (*SFP3ptr++) * Y_tmp;
            #endif

            while( ++Yptr != YptrEnd )
            {
                Y_tmp = *Yptr;
                x0 += (*SFP0ptr++) * Y_tmp;
                #if nIC>=2
                x1 += (*SFP1ptr++) * Y_tmp;
                #endif
                #if nIC>=3
                x2 += (*SFP2ptr++) * Y_tmp;
                #endif
                #if nIC>=4
                x3 += (*SFP3ptr++) * Y_tmp;
                #endif
            }

            w = (double)(*t_l);
            x[*t_f]      += w * x0;
            #if nIC>=2
            x[*t_f+nF]   += w * x1;
            #endif
            #if nIC>=3
            x[*t_f+2*nF] += w * x2;
            #endif
            #if nIC>=4
            x[*t_f+3*nF] += w * x3;
            #endif
        }

        t_f++;
        t_v++;
        t_o++;
        t_l++;
        t_t++;
    }
#endif

    /* extra-cellular compartments */
#if nEC>=1
    t_v    = ECv + ECthreads[id];
    t_vEnd = ECv + ECthreads[id+1];
    t_o    = ECo + ECthreads[id];

    x_Ptr0 = x + nIC*nF + ECthreads[id];
    #if nEC>=2
    x_Ptr1 = x_Ptr0 + nE;
    #endif
    #if nEC>=3
    x_Ptr2 = x_Ptr1 + nE;
    #endif
    #if nEC>=4
    x_Ptr3 = x_Ptr2 + nE;
    #endif

    while( t_v != t_vEnd )
    {
        Yptr    = Y    + (*t_v++);
        YptrEnd = Yptr + nS;
        offset  = nS * (*t_o++);

        Y_tmp = *Yptr;
        SFP0ptr = wmhSFP0 + offset;
        x0 = (*SFP0ptr++) * Y_tmp;
        #if nEC>=2
        SFP1ptr = wmhSFP1 + offset;
        x1 = (*SFP1ptr++) * Y_tmp;
        #endif
        #if nEC>=3
        SFP2ptr = wmhSFP2 + offset;
        x2 = (*SFP2ptr++) * Y_tmp;
        #endif
        #if nEC>=4
        SFP3ptr = wmhSFP3 + offset;
        x3 = (*SFP3ptr++) * Y_tmp;
        #endif

        while( ++Yptr != YptrEnd )
        {
            Y_tmp = *Yptr;
            x0 += (*SFP0ptr++) * Y_tmp;
            #if nEC>=2
            x1 += (*SFP1ptr++) * Y_tmp;
            #endif
            #if nEC>=3
            x2 += (*SFP2ptr++) * Y_tmp;
            #endif
            #if nEC>=4
            x3 += (*SFP3ptr++) * Y_tmp;
            #endif
        }
        (*x_Ptr0++) += x0;
        #if nEC>=2
        (*x_Ptr1++) += x1;
        #endif
        #if nEC>=3
        (*x_Ptr2++) += x2;
        #endif
        #if nEC>=4
        (*x_Ptr3++) += x3;
        #endif
    }
#endif

    /* isotropic compartments */
#if nISO>=1
    t_v    = ISOv + ISOthreads[id];
    t_vEnd = ISOv + ISOthreads[id+1];

    x_Ptr0 = x + nIC*nF + nEC*nE + ISOthreads[id];
    #if nISO>=2
    x_Ptr1 = x_Ptr0 + nV;
    #endif
    #if nISO>=3
    x_Ptr2 = x_Ptr1 + nV;
    #endif
    #if nISO>=4
    x_Ptr3 = x_Ptr2 + nV;
    #endif

    while( t_v != t_vEnd )
    {
        Yptr    = Y    + (*t_v++);
        YptrEnd = Yptr + nS;

        SFP0ptr = isoSFP0;
        #if nISO>=2
        SFP1ptr = isoSFP1;
        #endif
        #if nISO>=3
        SFP2ptr = isoSFP2;
        #endif
        #if nISO>=4
        SFP3ptr = isoSFP3;
        #endif

        Y_tmp = *Yptr;
        x0 = (*SFP0ptr++) * Y_tmp;
        #if nISO>=2
        x1 = (*SFP1ptr++) * Y_tmp;
        #endif
        #if nISO>=3
        x2 = (*SFP2ptr++) * Y_tmp;
        #endif
        #if nISO>=4
        x3 = (*SFP3ptr++) * Y_tmp;
        #endif

        while( ++Yptr != YptrEnd )
        {
            Y_tmp = *Yptr;
            x0  += (*SFP0ptr++) * Y_tmp;
            #if nISO>=2
            x1  += (*SFP1ptr++) * Y_tmp;
            #endif
            #if nISO>=3
            x2  += (*SFP2ptr++) * Y_tmp;
            #endif
            #if nISO>=4
            x3  += (*SFP3ptr++) * Y_tmp;
            #endif
        }

        (*x_Ptr0++) += x0;
        #if nISO>=2
        (*x_Ptr1++) += x1;
        #endif
        #if nISO>=3
        (*x_Ptr2++) += x2;
        #endif
        #if nISO>=4
        (*x_Ptr3++) += x3;
        #endif
    }
#endif

    pthread_exit( 0 );
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray* tmp;

    /* ================ */
    /* Check the INPUTS */
    /* ================ */
    #if DO_CHECK > 0
     if( nrhs != 4 )
        mexErrMsgIdAndTxt("InvalidInput:nrhs", "Require 4 inputs.");
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

    // Parse "DICTIONARY.IC.n"
    tmp = mxGetField( IC, 0, "n" );
    #if DO_CHECK > 0
    if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.IC.n","'DICTIONARY.IC.n' must be a real scalar");
    #endif
    n = mxGetScalar( tmp );

    // Parse "DICTIONARY.IC.nF"
    tmp = mxGetField( IC, 0, "nF" );
    #if DO_CHECK > 0
    if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.IC.nF","'DICTIONARY.IC.nF' must be a real scalar");
    #endif
    nF = mxGetScalar( tmp );

    // Parse "KERNELS.nS"
    tmp = mxGetField( prhs[1], 0, "nS" );
    #if DO_CHECK > 0
    if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.nS","'KERNELS.nS' must be a real scalar");
    #endif
    nS = mxGetScalar( tmp );

    // Parse "DICTIONARY.nV"
    tmp = mxGetField( prhs[0], 0, "nV" );
    #if DO_CHECK > 0
    if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.nV","'DICTIONARY.nV' must be a real scalar");
    #endif
    nV = mxGetScalar( tmp );

    // Parse "DICTIONARY.EC.nE"
    tmp = mxGetField( EC, 0, "nE" );
    #if DO_CHECK > 0
    if( !mxIsDouble(tmp) || mxIsComplex(tmp) ||  mxGetNumberOfElements(tmp)!=1 )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.EC.nE","'DICTIONARY.EC.nE' must be a real scalar");
    #endif
    nE = mxGetScalar( tmp );

    // Parse "DICTIONARY.IC.fiber"
    tmp = mxGetField( IC, 0, "fiber" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint32") )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.IC.fiber","'DICTIONARY.IC.fiber' must be a n*1 vector (uint32)");
    #endif
     ICf = (UINT32_T*) mxGetData( tmp );

    // Parse "DICTIONARY.IC.v", "DICTIONARY.IC.o"
    tmp = mxGetField( IC, 0, "v" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint32") )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.IC.v","'DICTIONARY.IC.v' must be a n*1 vector (uint32)");
    #endif
    ICv = (UINT32_T*) mxGetData(tmp);
    tmp = mxGetField( IC, 0, "o" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint16") )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.IC.o","'DICTIONARY.IC.o' must be a n*1 vector (uint16)");
    #endif
    ICo = (UINT16_T*) mxGetData(tmp);

    tmp = mxGetField( IC, 0, "len" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"single") )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.IC.len","'DICTIONARY.IC.len' must be a n*1 vector (single)");
    #endif
     ICl = (float*) mxGetData(tmp);

    // Parse "DICTIONARY.EC.v","DICTIONARY.EC.o"
    tmp = mxGetField( EC, 0, "v" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint32") )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.EC.v","'DICTIONARY.EC.v' must be a n*1 vector (uint32)");
    #endif
    ECv = (UINT32_T*) mxGetData(tmp);
    tmp = mxGetField( EC, 0, "o" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint16") )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.EC.o","'DICTIONARY.EC.o' must be a n*1 vector (uint16)");
    #endif
    ECo = (UINT16_T*) mxGetData(tmp);

    // Parse "DICTIONARY.ISO.v"
    tmp = mxGetField( ISO, 0, "v" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint32") )
        mexErrMsgIdAndTxt("InvalidInput:DICTIONARY.ISO.v","'DICTIONARY.ISO.v' must be a n*1 vector (uint32)");
    #endif
    ISOv = (UINT32_T*) mxGetData(tmp);

    // Parse "KERNELS.wmr", "KERNELS.wmh" and "KERNELS.iso"
    mxArray* wmr = mxGetField( prhs[1], 0, "wmr" );
    #if DO_CHECK > 0
    if ( !mxIsCell(wmr) )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmr","'KERNELS.wmr' must be a cell array");
    #endif
    #if nIC>=1
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(wmr,0),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmr","'KERNELS.wmr' must contain single");
    #endif
    wmrSFP0 = (float*) mxGetData( mxGetCell(wmr,0) );
    #endif
    #if nIC>=2
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(wmr,1),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmr","'KERNELS.wmr' must contain single");
    #endif
    wmrSFP1 = (float*) mxGetData( mxGetCell(wmr,1) );
    #endif
    #if nIC>=3
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(wmr,2),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmr","'KERNELS.wmr' must contain single");
    #endif
    wmrSFP2 = (float*) mxGetData( mxGetCell(wmr,2) );
    #endif
    #if nIC>=4
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(wmr,3),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmr","'KERNELS.wmr' must contain single");
    #endif
    wmrSFP3 = (float*) mxGetData( mxGetCell(wmr,3) );
    #endif

    #if nEC>=1
    mxArray* wmh = mxGetField( prhs[1], 0, "wmh" );
    #if DO_CHECK > 0
    if ( !mxIsCell(wmh) )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmh","'KERNELS.wmh' must be a cell array");
    if( !mxIsClass(mxGetCell(wmh,0),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmh","'KERNELS.wmh' must contain single");
    #endif
    wmhSFP0 = (float*) mxGetData( mxGetCell(wmh,0) );
    #if nEC>=2
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(wmh,1),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmh","'KERNELS.wmh' must contain single");
    #endif
    wmhSFP1 = (float*) mxGetData( mxGetCell(wmh,1) );
    #endif
    #if nEC>=3
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(wmh,2),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmh","'KERNELS.wmh' must contain single");
    #endif
    wmhSFP2 = (float*) mxGetData( mxGetCell(wmh,2) );
    #endif
    #if nEC>=4
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(wmh,3),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.wmh","'KERNELS.wmh' must contain single");
    #endif
    wmhSFP3 = (float*) mxGetData( mxGetCell(wmh,3) );
    #endif
    #endif

    #if nISO>=1
    mxArray* iso = mxGetField( prhs[1], 0, "iso" );
    #if DO_CHECK > 0
    if ( !mxIsCell(iso) )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.iso","'KERNELS.iso' must be a cell array");
    if( !mxIsClass(mxGetCell(iso,0),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.iso","'KERNELS.iso' must contain single");
    #endif
    isoSFP0 = (float*) mxGetData( mxGetCell(iso,0) );
    #if nISO>=2
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(iso,1),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.iso","'KERNELS.iso' must contain single");
    #endif
    isoSFP1 = (float*) mxGetData( mxGetCell(iso,1) );
    #endif
    #if nISO>=3
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(iso,2),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.iso","'KERNELS.iso' must contain single");
    #endif
    isoSFP2 = (float*) mxGetData( mxGetCell(iso,2) );
    #endif
    #if nISO>=4
    #if DO_CHECK > 0
    if( !mxIsClass(mxGetCell(iso,3),"single") )
        mexErrMsgIdAndTxt("InvalidInput:KERNELS.iso","'KERNELS.iso' must contain single");
    #endif
    isoSFP3 = (float*) mxGetData( mxGetCell(iso,3) );
    #endif
    #endif

    // Parse "Y"
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(prhs[2]) != 2 || mxGetN(prhs[2]) != 1 || !mxIsClass(prhs[2],"double") )
        mexErrMsgIdAndTxt("InvalidInput:Y","'Y' must be a n*1 vector (double)");
    #endif
     Y = (double*) mxGetData( prhs[2] );

    // parse "THREADS"
    #if DO_CHECK > 0
    if ( !mxIsStruct(prhs[3]) )
        mexErrMsgIdAndTxt("InvalidInput:THREADS", "'THREADS' must be a struct");
    #endif

    tmp  = mxGetField( prhs[3], 0, "ICt" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint8") )
        mexErrMsgIdAndTxt("InvalidInput:THREADS.ICt","'THREADS.ICt' must be a n*1 vector (uint8)");
    #endif
     ICthreads = (UINT8_T*) mxGetData( tmp );

    tmp  = mxGetField( prhs[3], 0, "ECt" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint32") )
        mexErrMsgIdAndTxt("InvalidInput:THREADS.ECt","'THREADS.ECt' must be a n*1 vector (uint32)");
    #endif
    ECthreads = (UINT32_T*) mxGetData( tmp );

    tmp  = mxGetField( prhs[3], 0, "ISOt" );
    #if DO_CHECK > 0
    if ( mxGetNumberOfDimensions(tmp) != 2 || mxGetN(tmp) != 1 || !mxIsClass(tmp,"uint32") )
        mexErrMsgIdAndTxt("InvalidInput:THREADS.ISOt","'THREADS.ISOt' must be a n*1 vector (uint32)");
    #endif
    ISOthreads = (UINT32_T*) mxGetData( tmp );


    /* =============== */
    /* Set the OUTPUTS */
    /* =============== */
    #if DO_CHECK > 0
    if( nlhs != 1 )
        mexErrMsgIdAndTxt("InvalidOutput:nlhs", "Required 1 output.");
    #endif

    const int outDims[4] = { nIC*nF + nEC*nE + nISO*nV, 1 };
    plhs[0] = mxCreateNumericArray(2, outDims, mxDOUBLE_CLASS, mxREAL);
    x = (double*)mxGetData( plhs[0] );


    /* ================================================== */
    /* Run SEPARATE THREADS to perform the multiplication */
    /* ================================================== */
    pthread_t threads[nTHREADS];
    for(int t=0; t<nTHREADS ; t++)
        pthread_create( &threads[t], NULL, computeProductBlock, (void *) (long int)t );
    for(int t=0; t<nTHREADS ; t++)
        pthread_join( threads[t], NULL );

    return;
}
