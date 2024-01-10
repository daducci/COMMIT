#include <pthread.h>
#include <stdint.h> // uint32_t etc

// number of THREADS
#ifdef nTHREADS
    #if (nTHREADS<1 || nTHREADS>255)
    #error "nTHREADS" must be in the range 1..255
    #endif
#else
    #error "nTHREADS" parameter must be passed to the compiler as "-DnTHREADS=<value>"
#endif


/* global variables */
int         nF, n;
double      *x, *Y;
uint32_t    *ICthreads, *ISOthreads;
uint8_t     *ICthreadsT;
uint32_t    *ISOthreadsT;
uint32_t    *ICf, *ICv, *ISOv;
float       *ICl;


// ====================================================
// Compute a sub-block of the A*x MAtRIX-VECTOR product
// ====================================================
void* COMMIT_A__block( void *ptr )
{
    int      id = (long)ptr;
    double   x0;
    double   *xPtr;
    uint32_t *t_v, *t_vEnd, *t_f;
    float    *t_l;

    // intra-cellular compartments
    t_v    = ICv + ICthreads[id];
    t_vEnd = ICv + ICthreads[id+1];
    t_l    = ICl + ICthreads[id];
    t_f    = ICf + ICthreads[id];

    while( t_v != t_vEnd )
    {
        x0 = x[*t_f];
        if ( x0 != 0 )
            Y[*t_v] += (double)(*t_l) * x0;
        t_f++;
        t_v++;
        t_l++;
    }
    #if nIC>=1
        // DCT basis functions 
        t_v    = ICv + ICthreads[id];
        t_vEnd = ICv + ICthreads[id+1];
        t_l    = ICl + ICthreads[id];
        t_f    = ICf + ICthreads[id];

        while( t_v != t_vEnd )
        {
            x_Ptr0 = x + *t_f;
            x0 = *x_Ptr0;
            #if nIC>=2
            x_Ptr1 = x_Ptr0 + nF;
            x1 = *x_Ptr1;
            #endif
            #if nIC>=3
            x_Ptr2 = x_Ptr1 + nF;
            x2 = *x_Ptr2;
            #endif
            #if nIC>=4
            x_Ptr3 = x_Ptr2 + nF;
            x3 = *x_Ptr3;
            #endif
            #if nIC>=5
            x_Ptr4 = x_Ptr3 + nF;
            x4 = *x_Ptr4;
            #endif
            #if nIC>=6
            x_Ptr5 = x_Ptr4 + nF;
            x5 = *x_Ptr5;
            #endif
            #if nIC>=7
            x_Ptr6 = x_Ptr5 + nF;
            x6 = *x_Ptr6;
            #endif
            #if nIC>=8
            x_Ptr7 = x_Ptr6 + nF;
            x7 = *x_Ptr7;
            #endif
            #if nIC>=9
            x_Ptr8 = x_Ptr7 + nF;
            x8 = *x_Ptr8;
            #endif
            #if nIC>=10
            x_Ptr9 = x_Ptr8 + nF;
            x9 = *x_Ptr9;
            #endif
            #if nIC>=11
            x_Ptr10 = x_Ptr9 + nF;
            x10 = *x_Ptr10;
            #endif
            #if nIC>=12
            x_Ptr11 = x_Ptr10 + nF;
            x11 = *x_Ptr11;
            #endif
            #if nIC>=13
            x_Ptr12 = x_Ptr11 + nF;
            x12 = *x_Ptr12;
            #endif
            #if nIC>=14
            x_Ptr13 = x_Ptr12 + nF;
            x13 = *x_Ptr13;
            #endif
            #if nIC>=15
            x_Ptr14 = x_Ptr13 + nF;
            x14 = *x_Ptr14;
            #endif
            #if nIC>=16
            x_Ptr15 = x_Ptr14 + nF;
            x15 = *x_Ptr15;
            #endif
            #if nIC>=17
            x_Ptr16 = x_Ptr15 + nF;
            x16 = *x_Ptr16;
            #endif
            #if nIC>=18
            x_Ptr17 = x_Ptr16 + nF;
            x17 = *x_Ptr17;
            #endif
            #if nIC>=19
            x_Ptr18 = x_Ptr17 + nF;
            x18 = *x_Ptr18;
            #endif
            #if nIC>=20
            x_Ptr19 = x_Ptr18 + nF;
            x19 = *x_Ptr19;
            #endif

            if ( x0 != 0
            #if nIC>=2
                || x1 != 0
            #endif
            #if nIC>=3
                || x2 != 0
            #endif
            #if nIC>=4
                || x3 != 0
            #endif
            )
            {
                //  Y[*t_v] += (double)(*t_l) * x0;
                Yptr    = Y    + nS * (*t_v);
                YptrEnd = Yptr + nS;
                w       = (double)(*t_l);
                offset  = nS * (*t_o);
                SFP0ptr = wmrSFP0 + offset;
                #if nIC>=2
                SFP1ptr = wmrSFP1 + offset;
                #endif
                #if nIC>=3
                SFP2ptr = wmrSFP2 + offset;
                #endif
                #if nIC>=4
                SFP3ptr = wmrSFP3 + offset;
                #endif
                
                while( Yptr != YptrEnd )
                    (*Yptr++) += w * (
                            x0 * (*SFP0ptr++)
                            #if nIC>=2
                            + x1 * (*SFP1ptr++)
                            #endif
                            #if nIC>=3
                            + x2 * (*SFP2ptr++)
                            #endif
                            #if nIC>=4
                            + x3 * (*SFP3ptr++)
                            #endif
                    );
            }

            t_f++;
            t_v++;
            t_l++;
        }
    #endif
    #if nISO>=1
        // isotropic compartments
        t_v    = ISOv + ISOthreads[id];
        t_vEnd = ISOv + ISOthreads[id+1];
        xPtr   = x + nF + ISOthreads[id];

        while( t_v != t_vEnd )
        {
            x0 = *xPtr++;
            if ( x0 != 0 )
                Y[*t_v] += x0;
            t_v++;
        }
    #endif

    pthread_exit( 0 );
}


// =========================
// Function called by CYTHON
// =========================
void COMMIT_A(
    int _nF, int _n, int _nE, int _nV, int _nS, int _ndirs,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICv, uint16_t *_ICo, float *_ICl,
    uint32_t *_ECv, uint16_t *_ECo,
    uint32_t *_ISOv,
    float *_wmrSFP, float *_ICmod, float *_wmhSFP, float *_isoSFP,
    uint32_t* _ICthreads, uint32_t* _ECthreads, uint32_t* _ISOthreads
)
{
    nF = _nF;
    n  = _n;

    x = _vIN;
    Y = _vOUT;

    ICf  = _ICf;
    ICv  = _ICv;
    ICl  = _ICl;
    ISOv = _ISOv;

    ICthreads  = _ICthreads;
    ISOthreads = _ISOthreads;

    // Run SEPARATE THREADS to perform the multiplication
    pthread_t threads[nTHREADS];
    int t;
    for(t=0; t<nTHREADS ; t++)
        pthread_create( &threads[t], NULL, COMMIT_A__block, (void *) (long int)t );
    for(t=0; t<nTHREADS ; t++)
        pthread_join( threads[t], NULL );
    return;
}


/* ===================================================== */
/* Compute a sub-block of the A'*y MAtRIX-VECTOR product */
/* ===================================================== */
void* COMMIT_At__block( void *ptr )
{
    int      id = (long)ptr;
    double   *xPtr;
    uint32_t *t_v, *t_vEnd, *t_f;
    float    *t_l;
    uint8_t  *t_t;

    // intra-cellular compartments
    t_v    = ICv;
    t_vEnd = ICv + n;
    t_l    = ICl;
    t_f    = ICf;
    t_t    = ICthreadsT;

    while( t_v != t_vEnd )
    {
        // in this case, I need to walk throug because the segments are ordered in "voxel order"
        if ( *t_t == id )
            x[*t_f] += (double)(*t_l) * Y[*t_v];
        t_t++;
        t_f++;
        t_v++;
        t_l++;
    }

#if nISO>=1
    // isotropic compartments
    t_v    = ISOv + ISOthreadsT[id];
    t_vEnd = ISOv + ISOthreadsT[id+1];
    xPtr   = x + nF + ISOthreadsT[id];

    while( t_v != t_vEnd )
        (*xPtr++) += Y[*t_v++];
#endif

    pthread_exit( 0 );
}


// =========================
// Function called by CYTHON
// =========================
void COMMIT_At(
    int _nF, int _n, int _nE, int _nV, int _nS, int _ndirs,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICv, uint16_t *_ICo, float *_ICl,
    uint32_t *_ECv, uint16_t *_ECo,
    uint32_t *_ISOv,
    float *_wmrSFP, float *_ICmod, float *_wmhSFP, float *_isoSFP,
    uint8_t* _ICthreadsT, uint32_t* _ECthreadsT, uint32_t* _ISOthreadsT
)
{
    nF = _nF;
    n  = _n;

    x = _vOUT;
    Y = _vIN;

    ICf  = _ICf;
    ICv  = _ICv;
    ICl  = _ICl;
    ISOv = _ISOv;

    ICthreadsT  = _ICthreadsT;
    ISOthreadsT = _ISOthreadsT;

    // Run SEPARATE THREADS to perform the multiplication
    pthread_t threads[nTHREADS];
    int t;
    for(t=0; t<nTHREADS ; t++)
        pthread_create( &threads[t], NULL, COMMIT_At__block, (void *) (long int)t );
    for(t=0; t<nTHREADS ; t++)
        pthread_join( threads[t], NULL );
    return;
}
