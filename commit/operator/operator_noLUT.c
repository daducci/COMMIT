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
    float *_wmrSFP, float *_wmhSFP, float *_isoSFP,
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
    float *_wmrSFP, float *_wmhSFP, float *_isoSFP,
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
