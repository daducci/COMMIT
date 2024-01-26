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
double      *x, *Y, *icSFB0, *icSFB1, *icSFB2, *icSFB3, *icSFB4, *icSFB5, *icSFB6, *icSFB7, *icSFB8, *icSFB9;
uint32_t    *ICthreads, *ISOthreads;
uint8_t     *ICthreadsT;
uint32_t    *ISOthreadsT;
uint32_t    *ICf, *ICv, *ISOv;
float       *ICl, *ICp;


// ====================================================
// Compute a sub-block of the A*x MAtRIX-VECTOR product
// ====================================================
void* COMMIT_A__block( void *ptr )
{
    int      id = (long)ptr;
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, w;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3, *x_Ptr4, *x_Ptr5, *x_Ptr6, *x_Ptr7, *x_Ptr8, *x_Ptr9;
    double   *Yptr, *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr;
    uint32_t *t_v, *t_vEnd, *t_f, *t_p;
    float    *t_l;
    int      start, offset;
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
            #if nICs>=5
            x_Ptr4 = x_Ptr3 + nF;
            x4 = *x_Ptr4;
            #endif
            #if nICs>=6
            x_Ptr5 = x_Ptr4 + nF;
            x5 = *x_Ptr5;
            #endif
            #if nICs>=7
            x_Ptr6 = x_Ptr5 + nF;
            x6 = *x_Ptr6;
            #endif
            #if nICs>=8
            x_Ptr7 = x_Ptr6 + nF;
            x7 = *x_Ptr7;
            #endif
            #if nICs>=9
            x_Ptr8 = x_Ptr7 + nF;
            x8 = *x_Ptr8;
            #endif
            #if nICs>=10
            x_Ptr9 = x_Ptr8 + nF;
            x9 = *x_Ptr9;
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
            #if nICs>=5
                || x4 != 0
            #endif
            #if nICs>=6
                || x5 != 0
            #endif
            #if nICs>=7
                || x6 != 0
            #endif
            #if nICs>=8
                || x7 != 0
            #endif
            #if nICs>=9
                || x8 != 0
            #endif
            #if nICs>=10
                || x9 != 0
            #endif

            )
            {
                Yptr    = Y + (*t_v);
                // YptrEnd = Yptr + nICs;
                w       = (double)(*t_l);
                offset  = (*t_p) - 1;
                // printf("offset = %d\n", offset);
                SFP0ptr = icSFB0 + offset; 
                #if nICs>=2
                SFP1ptr = icSFB1 + offset;
                #endif
                #if nICs>=3
                SFP2ptr = icSFB2 + offset;
                #endif
                #if nICs>=4
                SFP3ptr = icSFB3 + offset;
                #endif
                #if nICs>=5
                SFP4ptr = icSFB4 + offset;
                #endif
                #if nICs>=6
                SFP5ptr = icSFB5 + offset;
                #endif
                #if nICs>=7
                SFP6ptr = icSFB6 + offset;
                #endif
                #if nICs>=8
                SFP7ptr = icSFB7 + offset;
                #endif
                #if nICs>=9
                SFP8ptr = icSFB8 + offset;
                #endif
                #if nICs>=10
                SFP9ptr = icSFB9 + offset;
                #endif

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
                        #if nICs>=5
                        + x4 * (*SFP4ptr)
                        #endif
                        #if nICs>=6
                        + x5 * (*SFP5ptr)
                        #endif
                        #if nICs>=7
                        + x6 * (*SFP6ptr)
                        #endif
                        #if nICs>=8
                        + x7 * (*SFP7ptr)
                        #endif
                        #if nICs>=9
                        + x8 * (*SFP8ptr)
                        #endif
                        #if nICs>=10
                        + x9 * (*SFP9ptr)
                        #endif
                    );
                    printf("w = %f\n", w);
                    // printf("SFP0ptr = %f, SFP1ptr = %f, SFP2ptr = %f, SFP3ptr = %f, SFP4ptr = %f, SFP5ptr = %f, SFP6ptr = %f, SFP7ptr = %f, SFP8ptr = %f, SFP9ptr = %f\n", *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr);
                    // printf("SFP1ptr = %f, SFP2ptr = %f, SFP7ptr = %f, SFP9ptr = %f\n", *SFP1ptr, *SFP2ptr, *SFP7ptr, *SFP9ptr);

            }
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
    float *_wmrSFP, double *_ICmod, float *_wmhSFP, float *_isoSFP,
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
    icSFB1 = icSFB0 + nSf;
    #if nICs>=3
    icSFB2 = icSFB1 + nSf;
    #if nICs>=4
    icSFB3 = icSFB2 + nSf;
    #if nICs>=5
    icSFB4 = icSFB3 + nSf;
    #if nICs>=6
    icSFB5 = icSFB4 + nSf;
    #if nICs>=7
    icSFB6 = icSFB5 + nSf;
    #if nICs>=8
    icSFB7 = icSFB6 + nSf;
    #if nICs>=9
    icSFB8 = icSFB7 + nSf;
    #if nICs>=10
    icSFB9 = icSFB8 + nSf;
    #endif
    #endif
    #endif
    #endif
    #endif
    #endif
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
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, w, Y_tmp, *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr;
    double   *Yptr;
    uint32_t *t_v, *t_vEnd, *t_f, *t_p;
    float    *t_l;
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
                Y_tmp = *Yptr;

                offset = (*t_p);
                SFP0ptr   = icSFB0 + offset;
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
                #if nICs>=5
                SFP4ptr   = icSFB4 + offset;
                x4 = (*SFP4ptr) * Y_tmp;
                #endif
                #if nICs>=6
                SFP5ptr   = icSFB5 + offset;
                x5 = (*SFP5ptr) * Y_tmp;
                #endif
                #if nICs>=7
                SFP6ptr   = icSFB6 + offset;
                x6 = (*SFP6ptr) * Y_tmp;
                #endif
                #if nICs>=8
                SFP7ptr   = icSFB7 + offset;
                x7 = (*SFP7ptr) * Y_tmp;
                #endif
                #if nICs>=9
                SFP8ptr   = icSFB8 + offset;
                x8 = (*SFP8ptr) * Y_tmp;
                #endif
                #if nICs>=10
                SFP9ptr   = icSFB9 + offset;
                x9 = (*SFP9ptr) * Y_tmp;
                #endif

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
                #if nICs>=5
                x[*t_f+4*nF] += w * x4;
                #endif
                #if nICs>=6
                x[*t_f+5*nF] += w * x5;
                #endif
                #if nICs>=7
                x[*t_f+6*nF] += w * x6;
                #endif
                #if nICs>=8
                x[*t_f+7*nF] += w * x7;
                #endif
                #if nICs>=9
                x[*t_f+8*nF] += w * x8;
                #endif
                #if nICs>=10
                x[*t_f+9*nF] += w * x9;
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
    float *_wmrSFP, double *_ICmod, float *_wmhSFP, float *_isoSFP,
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
    #if nICs>=5
    icSFB4 = icSFB3 + nSf;
    #if nICs>=6
    icSFB5 = icSFB4 + nSf;
    #if nICs>=7
    icSFB6 = icSFB5 + nSf;
    #if nICs>=8
    icSFB7 = icSFB6 + nSf;
    #if nICs>=9
    icSFB8 = icSFB7 + nSf;
    #if nICs>=10
    icSFB9 = icSFB8 + nSf;
    #endif
    #endif
    #endif
    #endif
    #endif
    #endif
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