#include <pthread.h>
#include <stdint.h> // uint32_t etc
#include <stdio.h>  // printf
#include <inttypes.h>

// number of THREADS
#ifdef nTHREADS
    #if (nTHREADS<1 || nTHREADS>255)
    #error "nTHREADS" must be in the range 1..255
    #endif
#else
    #error "nTHREADS" parameter must be passed to the compiler as "-DnTHREADS=<value>"
#endif


/* global variables */
int         nF, n, nSf;
double      *x, *Y;
uint32_t    *ICthreads, *ISOthreads;
uint8_t     *ICthreadsT;
uint32_t    *ISOthreadsT;
uint32_t    *ICf, *ICv, *ISOv;
float       *ICl, *ICp;
float       *icSFB0, *icSFB1, *icSFB2, *icSFB3;


// ====================================================
// Compute a sub-block of the A*x MAtRIX-VECTOR product
// ====================================================
void* COMMIT_A__block( void *ptr )
{
    int      id = (long)ptr;
    double   x0, x1, x2, x3, w;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3;
    double   *Yptr, *YptrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr;
    uint32_t *t_v, *t_vEnd, *t_f;
    float    *t_l, *t_p;
    int      offset;
    // intra-cellular compartments
    #if nICs>=1
        // DCT basis functions 
        t_v    = ICv + ICthreads[id];
        t_vEnd = ICv + ICthreads[id+1];
        t_l    = ICl + ICthreads[id];
        t_f    = ICf + ICthreads[id];
        t_p    = ICp + ICthreads[id];      

        while( t_v != t_vEnd )
        {
            x_Ptr0 = x + *t_f;
            x0 = *x_Ptr0;
            #if nICs>=2
            x_Ptr1 = x_Ptr0 + nF;
            x1 = *x_Ptr1;
            #endif
            #if nICs>=3
            x_Ptr2 = x_Ptr1 + nF;
            x2 = *x_Ptr2;
            #endif
            #if nICs>=4
            x_Ptr3 = x_Ptr2 + nF;
            x3 = *x_Ptr3;
            #endif
            if ( x0 != 0
            #if nICs>=2
                || x1 != 0
            #endif
            #if nICs>=3
                || x2 != 0
            #endif
            #if nICs>=4
                || x3 != 0
            #endif
            )
            {
                Yptr    = Y + (*t_v);
                YptrEnd = Yptr + nICs;
                w       = (double)(*t_l);
                offset  = (*t_p);
                SFP0ptr = icSFB0;
                #if nICs>=2
                SFP1ptr = icSFB1 + offset;
                #endif
                #if nICs>=3
                SFP2ptr = icSFB2 + offset;
                #endif
                #if nICs>=4
                SFP3ptr = icSFB3 + offset;
                #endif
                
                // while( Yptr != YptrEnd )
                (*Yptr++) = w * (
                        x0 * (*SFP0ptr)
                        #if nICs>=2
                        + x1 * (*SFP1ptr)
                        #endif
                        #if nICs>=3
                        + x2 * (*SFP2ptr)
                        #endif
                        #if nICs>=4
                        + x3 * (*SFP3ptr)
                        #endif
                    );

            }
            // printf("x0 = %f, x1 = %f\n", x0, x1 );
            t_f++;
            t_v++;
            t_l++;
            t_p++;
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
    int _nF, int _n, int _nE, int _nV, int _nS, int _nSf, int _ndirs,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICv, uint16_t *_ICo, float *_ICl, float *_ICp,
    uint32_t *_ECv, uint16_t *_ECo,
    uint32_t *_ISOv,
    float *_wmrSFP, float *_ICmod, float *_wmhSFP, float *_isoSFP,
    uint32_t* _ICthreads, uint32_t* _ECthreads, uint32_t* _ISOthreads
)
{
    
    nF = _nF;
    n  = _n;
    nSf = _nSf;

    x = _vIN;
    Y = _vOUT;

    ICf  = _ICf;
    ICv  = _ICv;
    ICl  = _ICl;
    ICp  = _ICp;
    ISOv = _ISOv;
    
    #if nICs>=1
    icSFB0 = _ICmod;
    #if nICs>=2
    icSFB1 = icSFB0 + nSf-1;
    #if nICs>=3
    icSFB2 = icSFB1 + nSf-1;
    #if nICs>=4
    icSFB3 = icSFB2 + nSf-1;
    #endif
    #endif
    #endif
    #endif

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
    double   x0, x1, x2, x3, w, Y_tmp;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr;
    double   *Yptr, *YptrEnd;
    uint32_t *t_v, *t_vEnd, *t_f;
    float    *t_l, *t_p;
    uint8_t  *t_t;
    int      offset;
    double   *xPtr;

    #if ( nICs >= 1)
    // intra-cellular compartments
        t_v    = ICv;
        t_vEnd = ICv + n;
        t_l    = ICl;
        t_f    = ICf;
        t_p    = ICp;
        t_t    = ICthreadsT;
        

        while( t_v != t_vEnd )
        {
            // in this case, I need to walk throug because the segments are ordered in "voxel order"
            if ( *t_t == id )
            {
                Yptr    = Y + (*t_v);
                YptrEnd = Yptr + nICs;
                Y_tmp = *Yptr;

                offset = (*t_p);
                SFP0ptr   = icSFB0;
                x0 = (*SFP0ptr) * Y_tmp;
                #if nICs>=2
                SFP1ptr   = icSFB1 + offset;
                x1 = (*SFP1ptr) * Y_tmp;
                #endif
                #if nICs>=3
                SFP2ptr   = icSFB2 + offset;
                x2 = (*SFP2ptr) * Y_tmp;
                #endif
                #if nICs>=4
                SFP3ptr   = icSFB3 + offset;
                x3 = (*SFP3ptr) * Y_tmp;
                #endif

                // while( ++Yptr != YptrEnd )
                // {
                //     Y_tmp = *Yptr;
                x0 = (*SFP0ptr) * Y_tmp;
                #if nICs>=2
                x1 = (*SFP1ptr) * Y_tmp;
                #endif
                #if nICs>=3
                x2 = (*SFP2ptr) * Y_tmp;
                #endif
                #if nICs>=4
                x3 = (*SFP3ptr) * Y_tmp;
                #endif
                //     Yptr++;
                // }

                w = (double)(*t_l);
                x[*t_f]      += w * x0;
                #if nICs>=2
                x[*t_f+nF]   += w * x1;
                #endif
                #if nICs>=3
                x[*t_f+2*nF] += w * x2;
                #endif
                #if nICs>=4
                x[*t_f+3*nF] += w * x3;
                #endif
            }

            t_f++;
            t_v++;
            t_l++;
            t_p++;
            t_t++;
        }
    #endif
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
    int _nF, int _n, int _nE, int _nV, int _nS, int _nSf, int _ndirs,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICv, uint16_t *_ICo, float *_ICl, float *_ICp,
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
    ICp  = _ICp;
    ISOv = _ISOv;

    #if nICs>=1
    icSFB0 = _ICmod;
    #if nICs>=2
    icSFB1 = icSFB0 + nSf;
    #if nICs>=3
    icSFB2 = icSFB1 + nSf;
    #if nICs>=4
    icSFB3 = icSFB2 + nSf;
    #endif
    #endif
    #endif
    #endif


    ICthreadsT  = _ICthreadsT;
    ISOthreadsT = _ISOthreadsT;
    // printf("ICthreadsT = %d, ISOthreadsT = %d\n", ICthreadsT, ISOthreadsT);

    // Run SEPARATE THREADS to perform the multiplication
    pthread_t threads[nTHREADS];
    int t;
    for(t=0; t<nTHREADS ; t++)
        pthread_create( &threads[t], NULL, COMMIT_At__block, (void *) (long int)t );
    for(t=0; t<nTHREADS ; t++)
        pthread_join( threads[t], NULL );
    return;
}
