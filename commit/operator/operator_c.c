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
int         nF, n, nE, nV, nS;
double      *x, *Y;
uint32_t    *ICthreads, *ECthreads, *ISOthreads;
uint8_t     *ICthreadsT;
uint32_t    *ECthreadsT, *ISOthreadsT;
uint32_t    *ICf, *ICv, *ECv, *ISOv;
uint16_t    *ICo, *ECo;
float       *ICl;
float       *wmrSFP0, *wmrSFP1, *wmrSFP2, *wmrSFP3;
float       *wmhSFP0, *wmhSFP1, *wmhSFP2, *wmhSFP3;
float       *isoSFP0, *isoSFP1, *isoSFP2, *isoSFP3;



// ====================================================
// Compute a sub-block of the A*x MAtRIX-VECTOR product
// ====================================================
void* COMMIT_A__block( void *ptr )
{
    int      id = (long)ptr;
    int      offset;
    double   x0, x1, x2, x3, w;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3;
    double   *Yptr, *YptrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr;
    uint32_t *t_v, *t_vEnd, *t_f;
    uint16_t *t_o;
    float    *t_l;

#if nIC>=1
    // intra-cellular compartments
    t_v    = ICv + ICthreads[id];
    t_vEnd = ICv + ICthreads[id+1];
    t_o    = ICo + ICthreads[id];
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
        t_o++;
        t_l++;
    }
#endif

#if nEC>=1
    // extra-cellular compartments
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
        x0 = *x_Ptr0++;
        #if nEC>=2
        x1 = *x_Ptr1++;
        #endif
        #if nEC>=3
        x2 = *x_Ptr2++;
        #endif
        #if nEC>=4
        x3 = *x_Ptr3++;
        #endif
        if (
               x0 != 0
            #if nEC>=2
            || x1 != 0
            #endif
            #if nEC>=3
            || x2 != 0
            #endif
            #if nEC>=4
            || x3 != 0
            #endif
          )
        {
            Yptr    = Y    + nS * (*t_v);
            YptrEnd = Yptr + nS;
            offset  = nS * (*t_o);
            SFP0ptr = wmhSFP0 + offset;
            #if nEC>=2
            SFP1ptr = wmhSFP1 + offset;
            #endif
            #if nEC>=3
            SFP2ptr = wmhSFP2 + offset;
            #endif
            #if nEC>=4
            SFP3ptr = wmhSFP3 + offset;
            #endif

            while( Yptr != YptrEnd )
                (*Yptr++) += (
                      x0 * (*SFP0ptr++)
                    #if nEC>=2
                    + x1 * (*SFP1ptr++)
                    #endif
                    #if nEC>=3
                    + x2 * (*SFP2ptr++)
                    #endif
                    #if nEC>=4
                    + x3 * (*SFP3ptr++)
                    #endif
                );
        }
        t_v++;
        t_o++;
    }
#endif

#if nISO>=1
    // isotropic compartments
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
        x0 = *x_Ptr0++;
        #if nISO>=2
        x1 = *x_Ptr1++;
        #endif
        #if nISO>=3
        x2 = *x_Ptr2++;
        #endif
        #if nISO>=4
        x3 = *x_Ptr3++;
        #endif
        if (
               x0 != 0
            #if nISO>=2
            || x1 != 0
            #endif
            #if nISO>=3
            || x2 != 0
            #endif
            #if nISO>=4
            || x3 != 0
            #endif
          )
        {
            Yptr    = Y    + nS * (*t_v);
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

            while( Yptr != YptrEnd )
                (*Yptr++) += (
                      x0 * (*SFP0ptr++)
                    #if nISO>=2
                    + x1 * (*SFP1ptr++)
                    #endif
                    #if nISO>=3
                    + x2 * (*SFP2ptr++)
                    #endif
                    #if nISO>=4
                    + x3 * (*SFP3ptr++)
                    #endif
                );
        }
        t_v++;
    }
#endif

    pthread_exit( 0 );
}


// =========================
// Function called by CYTHON
// =========================
void COMMIT_A(
    int _nF, int _n, int _nE, int _nV, int _nS,
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
    nE = _nE;
    nV = _nV;
    nS = _nS;

    x = _vIN;
    Y = _vOUT;

    ICf  = _ICf;
    ICv  = _ICv;
    ICo  = _ICo;
    ICl  = _ICl;
    ECv  = _ECv;
    ECo  = _ECo;
    ISOv = _ISOv;

    #if nIC>=1
    wmrSFP0 = _wmrSFP;
    #if nIC>=2
    wmrSFP1 = wmrSFP0 + 181*181*_nS;
    #if nIC>=3
    wmrSFP2 = wmrSFP1 + 181*181*_nS;
    #if nIC>=4
    wmrSFP3 = wmrSFP2 + 181*181*_nS;
    #endif
    #endif
    #endif
    #endif
    #if nEC>=1
    wmhSFP0 = _wmhSFP;
    #if nEC>=2
    wmhSFP1 = wmhSFP0 + 181*181*_nS;
    #if nEC>=3
    wmhSFP2 = wmhSFP1 + 181*181*_nS;
    #if nEC>=4
    wmhSFP3 = wmhSFP2 + 181*181*_nS;
    #endif
    #endif
    #endif
    #endif
    #if nISO>=1
    isoSFP0 = _isoSFP;
    #if nISO>=2
    isoSFP1 = isoSFP0 + _nS;
    #if nISO>=3
    isoSFP2 = isoSFP1 + _nS;
    #if nISO>=4
    isoSFP3 = isoSFP2 + _nS;
    #endif
    #endif
    #endif
    #endif

    ICthreads  = _ICthreads;
    ECthreads  = _ECthreads;
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
    int      offset;
    double   x0, x1, x2, x3, w, Y_tmp;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3;
    double   *Yptr, *YptrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr;
    uint32_t *t_v, *t_vEnd, *t_f;
    uint16_t *t_o;
    float    *t_l;
    uint8_t  *t_t;

#if nIC>=1
    // intra-cellular compartments
    t_v    = ICv;
    t_vEnd = ICv + n;
    t_o    = ICo;
    t_l    = ICl;
    t_f    = ICf;
    t_t    = ICthreadsT;

    while( t_v != t_vEnd )
    {
        // in this case, I need to walk throug because the segments are ordered in "voxel order"
        if ( *t_t == id )
        {
            Yptr    = Y    + nS * (*t_v);
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

#if nEC>=1
    // extra-cellular compartments
    t_v    = ECv + ECthreadsT[id];
    t_vEnd = ECv + ECthreadsT[id+1];
    t_o    = ECo + ECthreadsT[id];

    x_Ptr0 = x + nIC*nF + ECthreadsT[id];
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
        Yptr    = Y    + nS * (*t_v++);
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

#if nISO>=1
    // isotropic compartments
    t_v    = ISOv + ISOthreadsT[id];
    t_vEnd = ISOv + ISOthreadsT[id+1];

    x_Ptr0 = x + nIC*nF + nEC*nE + ISOthreadsT[id];
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
        Yptr    = Y    + nS * (*t_v++);
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


// =========================
// Function called by CYTHON
// =========================
void COMMIT_At(
    int _nF, int _n, int _nE, int _nV, int _nS,
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
    nE = _nE;
    nV = _nV;
    nS = _nS;

    x = _vOUT;
    Y = _vIN;

    ICf  = _ICf;
    ICv  = _ICv;
    ICo  = _ICo;
    ICl  = _ICl;
    ECv  = _ECv;
    ECo  = _ECo;
    ISOv = _ISOv;

    #if nIC>=1
    wmrSFP0 = _wmrSFP;
    #if nIC>=2
    wmrSFP1 = wmrSFP0 + 181*181*_nS;
    #if nIC>=3
    wmrSFP2 = wmrSFP1 + 181*181*_nS;
    #if nIC>=4
    wmrSFP3 = wmrSFP2 + 181*181*_nS;
    #endif
    #endif
    #endif
    #endif
    #if nEC>=1
    wmhSFP0 = _wmhSFP;
    #if nEC>=2
    wmhSFP1 = wmhSFP0 + 181*181*_nS;
    #if nEC>=3
    wmhSFP2 = wmhSFP1 + 181*181*_nS;
    #if nEC>=4
    wmhSFP3 = wmhSFP2 + 181*181*_nS;
    #endif
    #endif
    #endif
    #endif
    #if nISO>=1
    isoSFP0 = _isoSFP;
    #if nISO>=2
    isoSFP1 = isoSFP0 + _nS;
    #if nISO>=3
    isoSFP2 = isoSFP1 + _nS;
    #if nISO>=4
    isoSFP3 = isoSFP2 + _nS;
    #endif
    #endif
    #endif
    #endif

    ICthreadsT  = _ICthreadsT;
    ECthreadsT  = _ECthreadsT;
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
