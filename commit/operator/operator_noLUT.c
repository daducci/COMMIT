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
    double   x0_tmp, x1_tmp, x2_tmp, x3_tmp, x4_tmp, x5_tmp, x6_tmp, x7_tmp, x8_tmp, x9_tmp;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3, *x_Ptr4, *x_Ptr5, *x_Ptr6, *x_Ptr7, *x_Ptr8, *x_Ptr9;
    double   *Yptr, *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr;
    uint32_t *t_v, *t_vEnd, *t_f, *t_p;
    float    *t_l;
    int      start_offset, offset;
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
            x0_tmp = 0;
            #if nICs>=2
            x_Ptr1 = x_Ptr0 + nF;
            x1 = *x_Ptr1;
            x1_tmp = 0;
            #endif
            #if nICs>=3
            x_Ptr2 = x_Ptr1 + nF;
            x2 = *x_Ptr2;
            x2_tmp = 0;
            #endif
            #if nICs>=4
            x_Ptr3 = x_Ptr2 + nF;
            x3 = *x_Ptr3;
            x3_tmp = 0;
            #endif
            #if nICs>=5
            x_Ptr4 = x_Ptr3 + nF;
            x4 = *x_Ptr4;
            x4_tmp = 0;
            #endif
            #if nICs>=6
            x_Ptr5 = x_Ptr4 + nF;
            x5 = *x_Ptr5;
            x5_tmp = 0;
            #endif
            #if nICs>=7
            x_Ptr6 = x_Ptr5 + nF;
            x6 = *x_Ptr6;
            x6_tmp = 0;
            #endif
            #if nICs>=8
            x_Ptr7 = x_Ptr6 + nF;
            x7 = *x_Ptr7;
            x7_tmp = 0;
            #endif
            #if nICs>=9
            x_Ptr8 = x_Ptr7 + nF;
            x8 = *x_Ptr8;
            x8_tmp = 0;
            #endif
            #if nICs>=10
            x_Ptr9 = x_Ptr8 + nF;
            x9 = *x_Ptr9;
            x9_tmp = 0;
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
            w       = (double)(*t_l);
            start_offset = 0;
            offset  = (*t_p);
            while (start_offset != offset)
            {
                // printf("start_offset = %d, offset = %d\n", start_offset, offset);
                SFP0ptr = icSFB0 + start_offset;
                x0_tmp += x0 * (*SFP0ptr);
                // printf("x0_tmp = %f\n", x0_tmp);
                // printf("SFP0ptr = %f\n", *SFP0ptr);
                #if nICs>=2
                SFP1ptr = icSFB1 + start_offset;
                x1_tmp += x1 * (*SFP1ptr);
                #endif
                #if nICs>=3
                SFP2ptr = icSFB2 + start_offset;
                x2_tmp += x2 * (*SFP2ptr);
                #endif
                #if nICs>=4
                SFP3ptr = icSFB3 + start_offset;
                x3_tmp += x3 * (*SFP3ptr);
                #endif
                #if nICs>=5
                SFP4ptr = icSFB4 + start_offset;
                x4_tmp += x4 * (*SFP4ptr);
                #endif
                #if nICs>=6
                SFP5ptr = icSFB5 + start_offset;
                x5_tmp += x5 * (*SFP5ptr);
                #endif
                #if nICs>=7
                SFP6ptr = icSFB6 + start_offset;
                x6_tmp += x6 * (*SFP6ptr);
                #endif
                #if nICs>=8
                SFP7ptr = icSFB7 + start_offset;
                x7_tmp += x7 * (*SFP7ptr);
                #endif
                #if nICs>=9
                SFP8ptr = icSFB8 + start_offset;
                x8_tmp += x8 * (*SFP8ptr);
                #endif
                #if nICs>=10
                SFP9ptr = icSFB9 + start_offset;
                x9_tmp += x9 * (*SFP9ptr);
                #endif                
                start_offset++;
                }

                (*Yptr++) = w * (
                        x0_tmp
                        #if nICs>=2
                        + x1_tmp
                        #endif
                        #if nICs>=3
                        + x2_tmp
                        #endif
                        #if nICs>=4
                        + x3_tmp
                        #endif
                        #if nICs>=5
                        + x4_tmp
                        #endif
                        #if nICs>=6
                        + x5_tmp
                        #endif
                        #if nICs>=7
                        + x6_tmp
                        #endif
                        #if nICs>=8
                        + x7_tmp
                        #endif
                        #if nICs>=9
                        + x8_tmp
                        #endif
                        #if nICs>=10
                        + x9_tmp
                        #endif
                    );
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
    int      start_offset, offset;
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
                start_offset = 0;
                offset = (*t_p);
                while (start_offset != offset)
                {
                    SFP0ptr = icSFB0 + start_offset;
                    x0 = (*SFP0ptr) * Y_tmp;
                    #if nICs>=2
                    SFP1ptr = icSFB1 + start_offset;
                    x1 = (*SFP1ptr) * Y_tmp;
                    #endif
                    #if nICs>=3
                    SFP2ptr = icSFB2 + start_offset;
                    x2 = (*SFP2ptr) * Y_tmp;
                    #endif
                    #if nICs>=4
                    SFP3ptr = icSFB3 + start_offset;
                    x3 = (*SFP3ptr) * Y_tmp;
                    #endif
                    #if nICs>=5
                    SFP4ptr = icSFB4 + start_offset;
                    x4 = (*SFP4ptr) * Y_tmp;
                    #endif
                    #if nICs>=6
                    SFP5ptr = icSFB5 + start_offset;
                    x5 = (*SFP5ptr) * Y_tmp;
                    #endif
                    #if nICs>=7
                    SFP6ptr = icSFB6 + start_offset;
                    x6 = (*SFP6ptr) * Y_tmp;
                    #endif
                    #if nICs>=8
                    SFP7ptr = icSFB7 + start_offset;
                    x7 = (*SFP7ptr) * Y_tmp;
                    #endif
                    #if nICs>=9
                    SFP8ptr = icSFB8 + start_offset;
                    x8 = (*SFP8ptr) * Y_tmp;
                    #endif
                    #if nICs>=10
                    SFP9ptr = icSFB9 + start_offset;
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
                    x[*t_f] += w * x3;
                    #endif
                    #if nICs>=5
                    x[*t_f] += w * x4;
                    #endif
                    #if nICs>=6
                    x[*t_f] += w * x5;
                    #endif
                    #if nICs>=7
                    x[*t_f] += w * x6;
                    #endif
                    #if nICs>=8
                    x[*t_f] += w * x7;
                    #endif
                    #if nICs>=9
                    x[*t_f] += w * x8;
                    #endif
                    #if nICs>=10
                    x[*t_f] += w * x9;
                    #endif
                    start_offset++;
                }
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
