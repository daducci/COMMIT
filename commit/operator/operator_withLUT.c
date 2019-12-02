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
int         nF, n, nE, nV, nS, ndirs;
double      *x, *Y;
uint32_t    *ICthreads, *ECthreads, *ISOthreads;
uint8_t     *ICthreadsT;
uint32_t    *ECthreadsT, *ISOthreadsT;
uint32_t    *ICf, *ICv, *ECv, *ISOv;
uint16_t    *ICo, *ECo;
float       *ICl;
float       *wmrSFP0, *wmrSFP1, *wmrSFP2, *wmrSFP3, *wmrSFP4, *wmrSFP5, *wmrSFP6, *wmrSFP7, *wmrSFP8, *wmrSFP9, *wmrSFP10, *wmrSFP11, *wmrSFP12, *wmrSFP13, *wmrSFP14, *wmrSFP15, *wmrSFP16, *wmrSFP17, *wmrSFP18, *wmrSFP19;
float       *wmhSFP0, *wmhSFP1, *wmhSFP2, *wmhSFP3, *wmhSFP4, *wmhSFP5, *wmhSFP6, *wmhSFP7, *wmhSFP8, *wmhSFP9, *wmhSFP10, *wmhSFP11, *wmhSFP12, *wmhSFP13, *wmhSFP14, *wmhSFP15, *wmhSFP16, *wmhSFP17, *wmhSFP18, *wmhSFP19;
float       *isoSFP0, *isoSFP1, *isoSFP2, *isoSFP3, *isoSFP4, *isoSFP5, *isoSFP6, *isoSFP7, *isoSFP8, *isoSFP9, *isoSFP10, *isoSFP11, *isoSFP12, *isoSFP13, *isoSFP14, *isoSFP15, *isoSFP16, *isoSFP17, *isoSFP18, *isoSFP19;



// ====================================================
// Compute a sub-block of the A*x MAtRIX-VECTOR product
// ====================================================
void* COMMIT_A__block( void *ptr )
{
    int      id = (long)ptr;
    int      offset;
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, w;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3, *x_Ptr4, *x_Ptr5, *x_Ptr6, *x_Ptr7, *x_Ptr8, *x_Ptr9, *x_Ptr10, *x_Ptr11, *x_Ptr12, *x_Ptr13, *x_Ptr14, *x_Ptr15, *x_Ptr16, *x_Ptr17, *x_Ptr18, *x_Ptr19;
    double   *Yptr, *YptrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr, *SFP10ptr, *SFP11ptr, *SFP12ptr, *SFP13ptr, *SFP14ptr, *SFP15ptr, *SFP16ptr, *SFP17ptr, *SFP18ptr, *SFP19ptr;
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
        #if nIC>=5
            || x4 != 0
        #endif
        #if nIC>=6
            || x5 != 0
        #endif
        #if nIC>=7
            || x6 != 0
        #endif
        #if nIC>=8
            || x7 != 0
        #endif
        #if nIC>=9
            || x8 != 0
        #endif
        #if nIC>=10
            || x9 != 0
        #endif
        #if nIC>=11
            || x10 != 0
        #endif
        #if nIC>=12
            || x11 != 0
        #endif
        #if nIC>=13
            || x12 != 0
        #endif
        #if nIC>=14
            || x13 != 0
        #endif
        #if nIC>=15
            || x14 != 0
        #endif
        #if nIC>=16
            || x15 != 0
        #endif
        #if nIC>=17
            || x16 != 0
        #endif
        #if nIC>=18
            || x17 != 0
        #endif
        #if nIC>=19
            || x18 != 0
        #endif
        #if nIC>=20
            || x19 != 0
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
            #if nIC>=5
            SFP4ptr = wmrSFP4 + offset;
            #endif
            #if nIC>=6
            SFP5ptr = wmrSFP5 + offset;
            #endif
            #if nIC>=7
            SFP6ptr = wmrSFP6 + offset;
            #endif
            #if nIC>=8
            SFP7ptr = wmrSFP7 + offset;
            #endif
            #if nIC>=9
            SFP8ptr = wmrSFP8 + offset;
            #endif
            #if nIC>=10
            SFP9ptr = wmrSFP9 + offset;
            #endif
            #if nIC>=11
            SFP10ptr = wmrSFP10 + offset;
            #endif
            #if nIC>=12
            SFP11ptr = wmrSFP11 + offset;
            #endif
            #if nIC>=13
            SFP12ptr = wmrSFP12 + offset;
            #endif
            #if nIC>=14
            SFP13ptr = wmrSFP13 + offset;
            #endif
            #if nIC>=15
            SFP14ptr = wmrSFP14 + offset;
            #endif
            #if nIC>=16
            SFP15ptr = wmrSFP15 + offset;
            #endif
            #if nIC>=17
            SFP16ptr = wmrSFP16 + offset;
            #endif
            #if nIC>=18
            SFP17ptr = wmrSFP17 + offset;
            #endif
            #if nIC>=19
            SFP18ptr = wmrSFP18 + offset;
            #endif
            #if nIC>=20
            SFP19ptr = wmrSFP19 + offset;
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
                        #if nIC>=5
                        + x4 * (*SFP4ptr++)
                        #endif
                        #if nIC>=6
                        + x5 * (*SFP5ptr++)
                        #endif
                        #if nIC>=7
                        + x6 * (*SFP6ptr++)
                        #endif
                        #if nIC>=8
                        + x7 * (*SFP7ptr++)
                        #endif
                        #if nIC>=9
                        + x8 * (*SFP8ptr++)
                        #endif
                        #if nIC>=10
                        + x9 * (*SFP9ptr++)
                        #endif
                        #if nIC>=11
                        + x10 * (*SFP10ptr++)
                        #endif
                        #if nIC>=12
                        + x11 * (*SFP11ptr++)
                        #endif
                        #if nIC>=13
                        + x12 * (*SFP12ptr++)
                        #endif
                        #if nIC>=14
                        + x13 * (*SFP13ptr++)
                        #endif
                        #if nIC>=15
                        + x14 * (*SFP14ptr++)
                        #endif
                        #if nIC>=16
                        + x15 * (*SFP15ptr++)
                        #endif
                        #if nIC>=17
                        + x16 * (*SFP16ptr++)
                        #endif
                        #if nIC>=18
                        + x17 * (*SFP17ptr++)
                        #endif
                        #if nIC>=19
                        + x18 * (*SFP18ptr++)
                        #endif
                        #if nIC>=20
                        + x19 * (*SFP19ptr++)
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
    #if nEC>=5
    x_Ptr4 = x_Ptr3 + nE;
    #endif
    #if nEC>=6
    x_Ptr5 = x_Ptr4 + nE;
    #endif
    #if nEC>=7
    x_Ptr6 = x_Ptr5 + nE;
    #endif
    #if nEC>=8
    x_Ptr7 = x_Ptr6 + nE;
    #endif
    #if nEC>=9
    x_Ptr8 = x_Ptr7 + nE;
    #endif
    #if nEC>=10
    x_Ptr9 = x_Ptr8 + nE;
    #endif
    #if nEC>=11
    x_Ptr10 = x_Ptr9 + nE;
    #endif
    #if nEC>=12
    x_Ptr11 = x_Ptr10 + nE;
    #endif
    #if nEC>=13
    x_Ptr12 = x_Ptr11 + nE;
    #endif
    #if nEC>=14
    x_Ptr13 = x_Ptr12 + nE;
    #endif
    #if nEC>=15
    x_Ptr14 = x_Ptr13 + nE;
    #endif
    #if nEC>=16
    x_Ptr15 = x_Ptr14 + nE;
    #endif
    #if nEC>=17
    x_Ptr16 = x_Ptr15 + nE;
    #endif
    #if nEC>=18
    x_Ptr17 = x_Ptr16 + nE;
    #endif
    #if nEC>=19
    x_Ptr18 = x_Ptr17 + nE;
    #endif
    #if nEC>=20
    x_Ptr19 = x_Ptr18 + nE;
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
        #if nEC>=5
        x4 = *x_Ptr4++;
        #endif
        #if nEC>=6
        x5 = *x_Ptr5++;
        #endif
        #if nEC>=7
        x6 = *x_Ptr6++;
        #endif
        #if nEC>=8
        x7 = *x_Ptr7++;
        #endif
        #if nEC>=9
        x8 = *x_Ptr8++;
        #endif
        #if nEC>=10
        x9 = *x_Ptr9++;
        #endif
        #if nEC>=11
        x10 = *x_Ptr10++;
        #endif
        #if nEC>=12
        x11 = *x_Ptr11++;
        #endif
        #if nEC>=13
        x12 = *x_Ptr12++;
        #endif
        #if nEC>=14
        x13 = *x_Ptr13++;
        #endif
        #if nEC>=15
        x14 = *x_Ptr14++;
        #endif
        #if nEC>=16
        x15 = *x_Ptr15++;
        #endif
        #if nEC>=17
        x16 = *x_Ptr16++;
        #endif
        #if nEC>=18
        x17 = *x_Ptr17++;
        #endif
        #if nEC>=19
        x18 = *x_Ptr18++;
        #endif
        #if nEC>=20
        x19 = *x_Ptr19++;
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
            #if nEC>=5
            || x4 != 0
            #endif
            #if nEC>=6
            || x5 != 0
            #endif
            #if nEC>=7
            || x6 != 0
            #endif
            #if nEC>=8
            || x7 != 0
            #endif
            #if nEC>=9
            || x8 != 0
            #endif
            #if nEC>=10
            || x9 != 0
            #endif
            #if nEC>=11
            || x10 != 0
            #endif
            #if nEC>=12
            || x11 != 0
            #endif
            #if nEC>=13
            || x12 != 0
            #endif
            #if nEC>=14
            || x13 != 0
            #endif
            #if nEC>=15
            || x14 != 0
            #endif
            #if nEC>=16
            || x15 != 0
            #endif
            #if nEC>=17
            || x16 != 0
            #endif
            #if nEC>=18
            || x17 != 0
            #endif
            #if nEC>=19
            || x18 != 0
            #endif
            #if nEC>=20
            || x19 != 0
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
            #if nEC>=5
            SFP4ptr = wmhSFP4 + offset;
            #endif
            #if nEC>=6
            SFP5ptr = wmhSFP5 + offset;
            #endif
            #if nEC>=7
            SFP6ptr = wmhSFP6 + offset;
            #endif
            #if nEC>=8
            SFP7ptr = wmhSFP7 + offset;
            #endif
            #if nEC>=9
            SFP8ptr = wmhSFP8 + offset;
            #endif
            #if nEC>=10
            SFP9ptr = wmhSFP9 + offset;
            #endif
            #if nEC>=11
            SFP10ptr = wmhSFP10 + offset;
            #endif
            #if nEC>=12
            SFP11ptr = wmhSFP11 + offset;
            #endif
            #if nEC>=13
            SFP12ptr = wmhSFP12 + offset;
            #endif
            #if nEC>=14
            SFP13ptr = wmhSFP13 + offset;
            #endif
            #if nEC>=15
            SFP14ptr = wmhSFP14 + offset;
            #endif
            #if nEC>=16
            SFP15ptr = wmhSFP15 + offset;
            #endif
            #if nEC>=17
            SFP16ptr = wmhSFP16 + offset;
            #endif
            #if nEC>=18
            SFP17ptr = wmhSFP17 + offset;
            #endif
            #if nEC>=19
            SFP18ptr = wmhSFP18 + offset;
            #endif
            #if nEC>=20
            SFP19ptr = wmhSFP19 + offset;
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
                    #if nEC>=5
                    + x4 * (*SFP4ptr++)
                    #endif
                    #if nEC>=6
                    + x5 * (*SFP5ptr++)
                    #endif
                    #if nEC>=7
                    + x6 * (*SFP6ptr++)
                    #endif
                    #if nEC>=8
                    + x7 * (*SFP7ptr++)
                    #endif
                    #if nEC>=9
                    + x8 * (*SFP8ptr++)
                    #endif
                    #if nEC>=10
                    + x9 * (*SFP9ptr++)
                    #endif
                    #if nEC>=11
                    + x10 * (*SFP10ptr++)
                    #endif
                    #if nEC>=12
                    + x11 * (*SFP11ptr++)
                    #endif
                    #if nEC>=13
                    + x12 * (*SFP12ptr++)
                    #endif
                    #if nEC>=14
                    + x13 * (*SFP13ptr++)
                    #endif
                    #if nEC>=15
                    + x14 * (*SFP14ptr++)
                    #endif
                    #if nEC>=16
                    + x15 * (*SFP15ptr++)
                    #endif
                    #if nEC>=17
                    + x16 * (*SFP16ptr++)
                    #endif
                    #if nEC>=18
                    + x17 * (*SFP17ptr++)
                    #endif
                    #if nEC>=19
                    + x18 * (*SFP18ptr++)
                    #endif
                    #if nEC>=20
                    + x19 * (*SFP19ptr++)
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
    #if nISO>=5
    x_Ptr4 = x_Ptr3 + nV;
    #endif
    #if nISO>=6
    x_Ptr5 = x_Ptr4 + nV;
    #endif
    #if nISO>=7
    x_Ptr6 = x_Ptr5 + nV;
    #endif
    #if nISO>=8
    x_Ptr7 = x_Ptr6 + nV;
    #endif
    #if nISO>=9
    x_Ptr8 = x_Ptr7 + nV;
    #endif
    #if nISO>=10
    x_Ptr9 = x_Ptr8 + nV;
    #endif
    #if nISO>=11
    x_Ptr10 = x_Ptr9 + nV;
    #endif
    #if nISO>=12
    x_Ptr11 = x_Ptr10 + nV;
    #endif
    #if nISO>=13
    x_Ptr12 = x_Ptr11 + nV;
    #endif
    #if nISO>=14
    x_Ptr13 = x_Ptr12 + nV;
    #endif
    #if nISO>=15
    x_Ptr14 = x_Ptr13 + nV;
    #endif
    #if nISO>=16
    x_Ptr15 = x_Ptr14 + nV;
    #endif
    #if nISO>=17
    x_Ptr16 = x_Ptr15 + nV;
    #endif
    #if nISO>=18
    x_Ptr17 = x_Ptr16 + nV;
    #endif
    #if nISO>=19
    x_Ptr18 = x_Ptr17 + nV;
    #endif
    #if nISO>=20
    x_Ptr19 = x_Ptr18 + nV;
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
        #if nISO>=5
        x4 = *x_Ptr4++;
        #endif
        #if nISO>=6
        x5 = *x_Ptr5++;
        #endif
        #if nISO>=7
        x6 = *x_Ptr6++;
        #endif
        #if nISO>=8
        x7 = *x_Ptr7++;
        #endif
        #if nISO>=9
        x8 = *x_Ptr8++;
        #endif
        #if nISO>=10
        x9 = *x_Ptr9++;
        #endif
        #if nISO>=11
        x10 = *x_Ptr10++;
        #endif
        #if nISO>=12
        x11 = *x_Ptr11++;
        #endif
        #if nISO>=13
        x12 = *x_Ptr12++;
        #endif
        #if nISO>=14
        x13 = *x_Ptr13++;
        #endif
        #if nISO>=15
        x14 = *x_Ptr14++;
        #endif
        #if nISO>=16
        x15 = *x_Ptr15++;
        #endif
        #if nISO>=17
        x16 = *x_Ptr16++;
        #endif
        #if nISO>=18
        x17 = *x_Ptr17++;
        #endif
        #if nISO>=19
        x18 = *x_Ptr18++;
        #endif
        #if nISO>=20
        x19 = *x_Ptr19++;
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
            #if nISO>=5
            || x4 != 0
            #endif
            #if nISO>=6
            || x5 != 0
            #endif
            #if nISO>=7
            || x6 != 0
            #endif
            #if nISO>=8
            || x7 != 0
            #endif
            #if nISO>=9
            || x8 != 0
            #endif
            #if nISO>=10
            || x9 != 0
            #endif
            #if nISO>=11
            || x10 != 0
            #endif
            #if nISO>=12
            || x11 != 0
            #endif
            #if nISO>=13
            || x12 != 0
            #endif
            #if nISO>=14
            || x13 != 0
            #endif
            #if nISO>=15
            || x14 != 0
            #endif
            #if nISO>=16
            || x15 != 0
            #endif
            #if nISO>=17
            || x16 != 0
            #endif
            #if nISO>=18
            || x17 != 0
            #endif
            #if nISO>=19
            || x18 != 0
            #endif
            #if nISO>=20
            || x19 != 0
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
            #if nISO>=5
            SFP4ptr = isoSFP4;
            #endif
            #if nISO>=6
            SFP5ptr = isoSFP5;
            #endif
            #if nISO>=7
            SFP6ptr = isoSFP6;
            #endif
            #if nISO>=8
            SFP7ptr = isoSFP7;
            #endif
            #if nISO>=9
            SFP8ptr = isoSFP8;
            #endif
            #if nISO>=10
            SFP9ptr = isoSFP9;
            #endif
            #if nISO>=11
            SFP10ptr = isoSFP10;
            #endif
            #if nISO>=12
            SFP11ptr = isoSFP11;
            #endif
            #if nISO>=13
            SFP12ptr = isoSFP12;
            #endif
            #if nISO>=14
            SFP13ptr = isoSFP13;
            #endif
            #if nISO>=15
            SFP14ptr = isoSFP14;
            #endif
            #if nISO>=16
            SFP15ptr = isoSFP15;
            #endif
            #if nISO>=17
            SFP16ptr = isoSFP16;
            #endif
            #if nISO>=18
            SFP17ptr = isoSFP17;
            #endif
            #if nISO>=19
            SFP18ptr = isoSFP18;
            #endif
            #if nISO>=20
            SFP19ptr = isoSFP19;
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
                    #if nISO>=5
                    + x4 * (*SFP4ptr++)
                    #endif
                    #if nISO>=6
                    + x5 * (*SFP5ptr++)
                    #endif
                    #if nISO>=7
                    + x6 * (*SFP6ptr++)
                    #endif
                    #if nISO>=8
                    + x7 * (*SFP7ptr++)
                    #endif
                    #if nISO>=9
                    + x8 * (*SFP8ptr++)
                    #endif
                    #if nISO>=10
                    + x9 * (*SFP9ptr++)
                    #endif
                    #if nISO>=11
                    + x10 * (*SFP10ptr++)
                    #endif
                    #if nISO>=12
                    + x11 * (*SFP11ptr++)
                    #endif
                    #if nISO>=13
                    + x12 * (*SFP12ptr++)
                    #endif
                    #if nISO>=14
                    + x13 * (*SFP13ptr++)
                    #endif
                    #if nISO>=15
                    + x14 * (*SFP14ptr++)
                    #endif
                    #if nISO>=16
                    + x15 * (*SFP15ptr++)
                    #endif
                    #if nISO>=17
                    + x16 * (*SFP16ptr++)
                    #endif
                    #if nISO>=18
                    + x17 * (*SFP17ptr++)
                    #endif
                    #if nISO>=19
                    + x18 * (*SFP18ptr++)
                    #endif
                    #if nISO>=20
                    + x19 * (*SFP19ptr++)
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
    nE = _nE;
    nV = _nV;
    nS = _nS;
    ndirs = _ndirs;

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
    wmrSFP1 = wmrSFP0 + _ndirs*_nS;
    #if nIC>=3
    wmrSFP2 = wmrSFP1 + _ndirs*_nS;
    #if nIC>=4
    wmrSFP3 = wmrSFP2 + _ndirs*_nS;
    #if nIC>=5
    wmrSFP4 = wmrSFP3 + _ndirs*_nS;
    #if nIC>=6
    wmrSFP5 = wmrSFP4 + _ndirs*_nS;
    #if nIC>=7
    wmrSFP6 = wmrSFP5 + _ndirs*_nS;
    #if nIC>=8
    wmrSFP7 = wmrSFP6 + _ndirs*_nS;
    #if nIC>=9
    wmrSFP8 = wmrSFP7 + _ndirs*_nS;
    #if nIC>=10
    wmrSFP9 = wmrSFP8 + _ndirs*_nS;
    #if nIC>=11
    wmrSFP10 = wmrSFP9 + _ndirs*_nS;
    #if nIC>=12
    wmrSFP11 = wmrSFP10 + _ndirs*_nS;
    #if nIC>=13
    wmrSFP12 = wmrSFP11 + _ndirs*_nS;
    #if nIC>=14
    wmrSFP13 = wmrSFP12 + _ndirs*_nS;
    #if nIC>=15
    wmrSFP14 = wmrSFP13 + _ndirs*_nS;
    #if nIC>=16
    wmrSFP15 = wmrSFP14 + _ndirs*_nS;
    #if nIC>=17
    wmrSFP16 = wmrSFP15 + _ndirs*_nS;
    #if nIC>=18
    wmrSFP17 = wmrSFP16 + _ndirs*_nS;
    #if nIC>=19
    wmrSFP18 = wmrSFP17 + _ndirs*_nS;
    #if nIC>=20
    wmrSFP19 = wmrSFP18 + _ndirs*_nS;
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
    #if nEC>=1
    wmhSFP0 = _wmhSFP;
    #if nEC>=2
    wmhSFP1 = wmhSFP0 + _ndirs*_nS;
    #if nEC>=3
    wmhSFP2 = wmhSFP1 + _ndirs*_nS;
    #if nEC>=4
    wmhSFP3 = wmhSFP2 + _ndirs*_nS;
    #if nEC>=5
    wmhSFP4 = wmhSFP3 + _ndirs*_nS;
    #if nEC>=6
    wmhSFP5 = wmhSFP4 + _ndirs*_nS;
    #if nEC>=7
    wmhSFP6 = wmhSFP5 + _ndirs*_nS;
    #if nEC>=8
    wmhSFP7 = wmhSFP6 + _ndirs*_nS;
    #if nEC>=9
    wmhSFP8 = wmhSFP7 + _ndirs*_nS;
    #if nEC>=10
    wmhSFP9 = wmhSFP8 + _ndirs*_nS;
    #if nEC>=11
    wmhSFP10 = wmhSFP9 + _ndirs*_nS;
    #if nEC>=12
    wmhSFP11 = wmhSFP10 + _ndirs*_nS;
    #if nEC>=13
    wmhSFP12 = wmhSFP11 + _ndirs*_nS;
    #if nEC>=14
    wmhSFP13 = wmhSFP12 + _ndirs*_nS;
    #if nEC>=15
    wmhSFP14 = wmhSFP13 + _ndirs*_nS;
    #if nEC>=16
    wmhSFP15 = wmhSFP14 + _ndirs*_nS;
    #if nEC>=17
    wmhSFP16 = wmhSFP15 + _ndirs*_nS;
    #if nEC>=18
    wmhSFP17 = wmhSFP16 + _ndirs*_nS;
    #if nEC>=19
    wmhSFP18 = wmhSFP17 + _ndirs*_nS;
    #if nEC>=20
    wmhSFP19 = wmhSFP18 + _ndirs*_nS;
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
    #if nISO>=1
    isoSFP0 = _isoSFP;
    #if nISO>=2
    isoSFP1 = isoSFP0 + _nS;
    #if nISO>=3
    isoSFP2 = isoSFP1 + _nS;
    #if nISO>=4
    isoSFP3 = isoSFP2 + _nS;
    #if nISO>=5
    isoSFP4 = isoSFP3 + _nS;
    #if nISO>=6
    isoSFP5 = isoSFP4 + _nS;
    #if nISO>=7
    isoSFP6 = isoSFP5 + _nS;
    #if nISO>=8
    isoSFP7 = isoSFP6 + _nS;
    #if nISO>=9
    isoSFP8 = isoSFP7 + _nS;
    #if nISO>=10
    isoSFP9 = isoSFP8 + _nS;
    #if nISO>=11
    isoSFP10 = isoSFP9 + _nS;
    #if nISO>=12
    isoSFP11 = isoSFP10 + _nS;
    #if nISO>=13
    isoSFP12 = isoSFP11 + _nS;
    #if nISO>=14
    isoSFP13 = isoSFP12 + _nS;
    #if nISO>=15
    isoSFP14 = isoSFP13 + _nS;
    #if nISO>=16
    isoSFP15 = isoSFP14 + _nS;
    #if nISO>=17
    isoSFP16 = isoSFP15 + _nS;
    #if nISO>=18
    isoSFP17 = isoSFP16 + _nS;
    #if nISO>=19
    isoSFP18 = isoSFP17 + _nS;
    #if nISO>=20
    isoSFP19 = isoSFP18 + _nS;
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
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, w, Y_tmp;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3, *x_Ptr4, *x_Ptr5, *x_Ptr6, *x_Ptr7, *x_Ptr8, *x_Ptr9, *x_Ptr10, *x_Ptr11, *x_Ptr12, *x_Ptr13, *x_Ptr14, *x_Ptr15, *x_Ptr16, *x_Ptr17, *x_Ptr18, *x_Ptr19;
    double   *Yptr, *YptrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr, *SFP10ptr, *SFP11ptr, *SFP12ptr, *SFP13ptr, *SFP14ptr, *SFP15ptr, *SFP16ptr, *SFP17ptr, *SFP18ptr, *SFP19ptr;
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
            #if nIC>=5
            SFP4ptr   = wmrSFP4 + offset;
            x4 = (*SFP4ptr++) * Y_tmp;
            #endif
            #if nIC>=6
            SFP5ptr   = wmrSFP5 + offset;
            x5 = (*SFP5ptr++) * Y_tmp;
            #endif
            #if nIC>=7
            SFP6ptr   = wmrSFP6 + offset;
            x6 = (*SFP6ptr++) * Y_tmp;
            #endif
            #if nIC>=8
            SFP7ptr   = wmrSFP7 + offset;
            x7 = (*SFP7ptr++) * Y_tmp;
            #endif
            #if nIC>=9
            SFP8ptr   = wmrSFP8 + offset;
            x8 = (*SFP8ptr++) * Y_tmp;
            #endif
            #if nIC>=10
            SFP9ptr   = wmrSFP9 + offset;
            x9 = (*SFP9ptr++) * Y_tmp;
            #endif
            #if nIC>=11
            SFP10ptr   = wmrSFP10 + offset;
            x10 = (*SFP10ptr++) * Y_tmp;
            #endif
            #if nIC>=12
            SFP11ptr   = wmrSFP11 + offset;
            x11 = (*SFP11ptr++) * Y_tmp;
            #endif
            #if nIC>=13
            SFP12ptr   = wmrSFP12 + offset;
            x12 = (*SFP12ptr++) * Y_tmp;
            #endif
            #if nIC>=14
            SFP13ptr   = wmrSFP13 + offset;
            x13 = (*SFP13ptr++) * Y_tmp;
            #endif
            #if nIC>=15
            SFP14ptr   = wmrSFP14 + offset;
            x14 = (*SFP14ptr++) * Y_tmp;
            #endif
            #if nIC>=16
            SFP15ptr   = wmrSFP15 + offset;
            x15 = (*SFP15ptr++) * Y_tmp;
            #endif
            #if nIC>=17
            SFP16ptr   = wmrSFP16 + offset;
            x16 = (*SFP16ptr++) * Y_tmp;
            #endif
            #if nIC>=18
            SFP17ptr   = wmrSFP17 + offset;
            x17 = (*SFP17ptr++) * Y_tmp;
            #endif
            #if nIC>=19
            SFP18ptr   = wmrSFP18 + offset;
            x18 = (*SFP18ptr++) * Y_tmp;
            #endif
            #if nIC>=20
            SFP19ptr   = wmrSFP19 + offset;
            x19 = (*SFP19ptr++) * Y_tmp;
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
                #if nIC>=5
                x4 += (*SFP4ptr++) * Y_tmp;
                #endif
                #if nIC>=6
                x5 += (*SFP5ptr++) * Y_tmp;
                #endif
                #if nIC>=7
                x6 += (*SFP6ptr++) * Y_tmp;
                #endif
                #if nIC>=8
                x7 += (*SFP7ptr++) * Y_tmp;
                #endif
                #if nIC>=9
                x8 += (*SFP8ptr++) * Y_tmp;
                #endif
                #if nIC>=10
                x9 += (*SFP9ptr++) * Y_tmp;
                #endif
                #if nIC>=11
                x10 += (*SFP10ptr++) * Y_tmp;
                #endif
                #if nIC>=12
                x11 += (*SFP11ptr++) * Y_tmp;
                #endif
                #if nIC>=13
                x12 += (*SFP12ptr++) * Y_tmp;
                #endif
                #if nIC>=14
                x13 += (*SFP13ptr++) * Y_tmp;
                #endif
                #if nIC>=15
                x14 += (*SFP14ptr++) * Y_tmp;
                #endif
                #if nIC>=16
                x15 += (*SFP15ptr++) * Y_tmp;
                #endif
                #if nIC>=17
                x16 += (*SFP16ptr++) * Y_tmp;
                #endif
                #if nIC>=18
                x17 += (*SFP17ptr++) * Y_tmp;
                #endif
                #if nIC>=19
                x18 += (*SFP18ptr++) * Y_tmp;
                #endif
                #if nIC>=20
                x19 += (*SFP19ptr++) * Y_tmp;
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
            #if nIC>=5
            x[*t_f+4*nF] += w * x4;
            #endif
            #if nIC>=6
            x[*t_f+5*nF] += w * x5;
            #endif
            #if nIC>=7
            x[*t_f+6*nF] += w * x6;
            #endif
            #if nIC>=8
            x[*t_f+7*nF] += w * x7;
            #endif
            #if nIC>=9
            x[*t_f+8*nF] += w * x8;
            #endif
            #if nIC>=10
            x[*t_f+9*nF] += w * x9;
            #endif
            #if nIC>=11
            x[*t_f+10*nF] += w * x10;
            #endif
            #if nIC>=12
            x[*t_f+11*nF] += w * x11;
            #endif
            #if nIC>=13
            x[*t_f+12*nF] += w * x12;
            #endif
            #if nIC>=14
            x[*t_f+13*nF] += w * x13;
            #endif
            #if nIC>=15
            x[*t_f+14*nF] += w * x14;
            #endif
            #if nIC>=16
            x[*t_f+15*nF] += w * x15;
            #endif
            #if nIC>=17
            x[*t_f+16*nF] += w * x16;
            #endif
            #if nIC>=18
            x[*t_f+17*nF] += w * x17;
            #endif
            #if nIC>=19
            x[*t_f+18*nF] += w * x18;
            #endif
            #if nIC>=20
            x[*t_f+19*nF] += w * x19;
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
    #if nEC>=5
    x_Ptr4 = x_Ptr3 + nE;
    #endif
    #if nEC>=6
    x_Ptr5 = x_Ptr4 + nE;
    #endif
    #if nEC>=7
    x_Ptr6 = x_Ptr5 + nE;
    #endif
    #if nEC>=8
    x_Ptr7 = x_Ptr6 + nE;
    #endif
    #if nEC>=9
    x_Ptr8 = x_Ptr7 + nE;
    #endif
    #if nEC>=10
    x_Ptr9 = x_Ptr8 + nE;
    #endif
    #if nEC>=11
    x_Ptr10 = x_Ptr9 + nE;
    #endif
    #if nEC>=12
    x_Ptr11 = x_Ptr10 + nE;
    #endif
    #if nEC>=13
    x_Ptr12 = x_Ptr11 + nE;
    #endif
    #if nEC>=14
    x_Ptr13 = x_Ptr12 + nE;
    #endif
    #if nEC>=15
    x_Ptr14 = x_Ptr13 + nE;
    #endif
    #if nEC>=16
    x_Ptr15 = x_Ptr14 + nE;
    #endif
    #if nEC>=17
    x_Ptr16 = x_Ptr15 + nE;
    #endif
    #if nEC>=18
    x_Ptr17 = x_Ptr16 + nE;
    #endif
    #if nEC>=19
    x_Ptr18 = x_Ptr17 + nE;
    #endif
    #if nEC>=20
    x_Ptr19 = x_Ptr18 + nE;
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
        #if nEC>=5
        SFP4ptr = wmhSFP4 + offset;
        x4 = (*SFP4ptr++) * Y_tmp;
        #endif
        #if nEC>=6
        SFP5ptr = wmhSFP5 + offset;
        x5 = (*SFP5ptr++) * Y_tmp;
        #endif
        #if nEC>=7
        SFP6ptr = wmhSFP6 + offset;
        x6 = (*SFP6ptr++) * Y_tmp;
        #endif
        #if nEC>=8
        SFP7ptr = wmhSFP7 + offset;
        x7 = (*SFP7ptr++) * Y_tmp;
        #endif
        #if nEC>=9
        SFP8ptr = wmhSFP8 + offset;
        x8 = (*SFP8ptr++) * Y_tmp;
        #endif
        #if nEC>=10
        SFP9ptr = wmhSFP9 + offset;
        x9 = (*SFP9ptr++) * Y_tmp;
        #endif
        #if nEC>=11
        SFP10ptr = wmhSFP10 + offset;
        x10 = (*SFP10ptr++) * Y_tmp;
        #endif
        #if nEC>=12
        SFP11ptr = wmhSFP11 + offset;
        x11 = (*SFP11ptr++) * Y_tmp;
        #endif
        #if nEC>=13
        SFP12ptr = wmhSFP12 + offset;
        x12 = (*SFP12ptr++) * Y_tmp;
        #endif
        #if nEC>=14
        SFP13ptr = wmhSFP13 + offset;
        x13 = (*SFP13ptr++) * Y_tmp;
        #endif
        #if nEC>=15
        SFP14ptr = wmhSFP14 + offset;
        x14 = (*SFP14ptr++) * Y_tmp;
        #endif
        #if nEC>=16
        SFP15ptr = wmhSFP15 + offset;
        x15 = (*SFP15ptr++) * Y_tmp;
        #endif
        #if nEC>=17
        SFP16ptr = wmhSFP16 + offset;
        x16 = (*SFP16ptr++) * Y_tmp;
        #endif
        #if nEC>=18
        SFP17ptr = wmhSFP17 + offset;
        x17 = (*SFP17ptr++) * Y_tmp;
        #endif
        #if nEC>=19
        SFP18ptr = wmhSFP18 + offset;
        x18 = (*SFP18ptr++) * Y_tmp;
        #endif
        #if nEC>=20
        SFP19ptr = wmhSFP19 + offset;
        x19 = (*SFP19ptr++) * Y_tmp;
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
            #if nEC>=5
            x4 += (*SFP4ptr++) * Y_tmp;
            #endif
            #if nEC>=6
            x5 += (*SFP5ptr++) * Y_tmp;
            #endif
            #if nEC>=7
            x6 += (*SFP6ptr++) * Y_tmp;
            #endif
            #if nEC>=8
            x7 += (*SFP7ptr++) * Y_tmp;
            #endif
            #if nEC>=9
            x8 += (*SFP8ptr++) * Y_tmp;
            #endif
            #if nEC>=10
            x9 += (*SFP9ptr++) * Y_tmp;
            #endif
            #if nEC>=11
            x10 += (*SFP10ptr++) * Y_tmp;
            #endif
            #if nEC>=12
            x11 += (*SFP11ptr++) * Y_tmp;
            #endif
            #if nEC>=13
            x12 += (*SFP12ptr++) * Y_tmp;
            #endif
            #if nEC>=14
            x13 += (*SFP13ptr++) * Y_tmp;
            #endif
            #if nEC>=15
            x14 += (*SFP14ptr++) * Y_tmp;
            #endif
            #if nEC>=16
            x15 += (*SFP15ptr++) * Y_tmp;
            #endif
            #if nEC>=17
            x16 += (*SFP16ptr++) * Y_tmp;
            #endif
            #if nEC>=18
            x17 += (*SFP17ptr++) * Y_tmp;
            #endif
            #if nEC>=19
            x18 += (*SFP18ptr++) * Y_tmp;
            #endif
            #if nEC>=20
            x19 += (*SFP19ptr++) * Y_tmp;
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
        #if nEC>=5
        (*x_Ptr4++) += x4;
        #endif
        #if nEC>=6
        (*x_Ptr5++) += x5;
        #endif
        #if nEC>=7
        (*x_Ptr6++) += x6;
        #endif
        #if nEC>=8
        (*x_Ptr7++) += x7;
        #endif
        #if nEC>=9
        (*x_Ptr8++) += x8;
        #endif
        #if nEC>=10
        (*x_Ptr9++) += x9;
        #endif
        #if nEC>=11
        (*x_Ptr10++) += x10;
        #endif
        #if nEC>=12
        (*x_Ptr11++) += x11;
        #endif
        #if nEC>=13
        (*x_Ptr12++) += x12;
        #endif
        #if nEC>=14
        (*x_Ptr13++) += x13;
        #endif
        #if nEC>=15
        (*x_Ptr14++) += x14;
        #endif
        #if nEC>=16
        (*x_Ptr15++) += x15;
        #endif
        #if nEC>=17
        (*x_Ptr16++) += x16;
        #endif
        #if nEC>=18
        (*x_Ptr17++) += x17;
        #endif
        #if nEC>=19
        (*x_Ptr18++) += x18;
        #endif
        #if nEC>=20
        (*x_Ptr19++) += x19;
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
    #if nISO>=5
    x_Ptr4 = x_Ptr3 + nV;
    #endif
    #if nISO>=6
    x_Ptr5 = x_Ptr4 + nV;
    #endif
    #if nISO>=7
    x_Ptr6 = x_Ptr5 + nV;
    #endif
    #if nISO>=8
    x_Ptr7 = x_Ptr6 + nV;
    #endif
    #if nISO>=9
    x_Ptr8 = x_Ptr7 + nV;
    #endif
    #if nISO>=10
    x_Ptr9 = x_Ptr8 + nV;
    #endif
    #if nISO>=11
    x_Ptr10 = x_Ptr9 + nV;
    #endif
    #if nISO>=12
    x_Ptr11 = x_Ptr10 + nV;
    #endif
    #if nISO>=13
    x_Ptr12 = x_Ptr11 + nV;
    #endif
    #if nISO>=14
    x_Ptr13 = x_Ptr12 + nV;
    #endif
    #if nISO>=15
    x_Ptr14 = x_Ptr13 + nV;
    #endif
    #if nISO>=16
    x_Ptr15 = x_Ptr14 + nV;
    #endif
    #if nISO>=17
    x_Ptr16 = x_Ptr15 + nV;
    #endif
    #if nISO>=18
    x_Ptr17 = x_Ptr16 + nV;
    #endif
    #if nISO>=19
    x_Ptr18 = x_Ptr17 + nV;
    #endif
    #if nISO>=20
    x_Ptr19 = x_Ptr18 + nV;
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
        #if nISO>=5
        SFP4ptr = isoSFP4;
        #endif
        #if nISO>=6
        SFP5ptr = isoSFP5;
        #endif
        #if nISO>=7
        SFP6ptr = isoSFP6;
        #endif
        #if nISO>=8
        SFP7ptr = isoSFP7;
        #endif
        #if nISO>=9
        SFP8ptr = isoSFP8;
        #endif
        #if nISO>=10
        SFP9ptr = isoSFP9;
        #endif
        #if nISO>=11
        SFP10ptr = isoSFP10;
        #endif
        #if nISO>=12
        SFP11ptr = isoSFP11;
        #endif
        #if nISO>=13
        SFP12ptr = isoSFP12;
        #endif
        #if nISO>=14
        SFP13ptr = isoSFP13;
        #endif
        #if nISO>=15
        SFP14ptr = isoSFP14;
        #endif
        #if nISO>=16
        SFP15ptr = isoSFP15;
        #endif
        #if nISO>=17
        SFP16ptr = isoSFP16;
        #endif
        #if nISO>=18
        SFP17ptr = isoSFP17;
        #endif
        #if nISO>=19
        SFP18ptr = isoSFP18;
        #endif
        #if nISO>=20
        SFP19ptr = isoSFP19;
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
        #if nISO>=5
        x4 = (*SFP4ptr++) * Y_tmp;
        #endif
        #if nISO>=6
        x5 = (*SFP5ptr++) * Y_tmp;
        #endif
        #if nISO>=7
        x6 = (*SFP6ptr++) * Y_tmp;
        #endif
        #if nISO>=8
        x7 = (*SFP7ptr++) * Y_tmp;
        #endif
        #if nISO>=9
        x8 = (*SFP8ptr++) * Y_tmp;
        #endif
        #if nISO>=10
        x9 = (*SFP9ptr++) * Y_tmp;
        #endif
        #if nISO>=11
        x10 = (*SFP10ptr++) * Y_tmp;
        #endif
        #if nISO>=12
        x11 = (*SFP11ptr++) * Y_tmp;
        #endif
        #if nISO>=13
        x12 = (*SFP12ptr++) * Y_tmp;
        #endif
        #if nISO>=14
        x13 = (*SFP13ptr++) * Y_tmp;
        #endif
        #if nISO>=15
        x14 = (*SFP14ptr++) * Y_tmp;
        #endif
        #if nISO>=16
        x15 = (*SFP15ptr++) * Y_tmp;
        #endif
        #if nISO>=17
        x16 = (*SFP16ptr++) * Y_tmp;
        #endif
        #if nISO>=18
        x17 = (*SFP17ptr++) * Y_tmp;
        #endif
        #if nISO>=19
        x18 = (*SFP18ptr++) * Y_tmp;
        #endif
        #if nISO>=20
        x19 = (*SFP19ptr++) * Y_tmp;
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
            #if nISO>=5
            x4  += (*SFP4ptr++) * Y_tmp;
            #endif
            #if nISO>=6
            x5  += (*SFP5ptr++) * Y_tmp;
            #endif
            #if nISO>=7
            x6  += (*SFP6ptr++) * Y_tmp;
            #endif
            #if nISO>=8
            x7  += (*SFP7ptr++) * Y_tmp;
            #endif
            #if nISO>=9
            x8  += (*SFP8ptr++) * Y_tmp;
            #endif
            #if nISO>=10
            x9  += (*SFP9ptr++) * Y_tmp;
            #endif
            #if nISO>=11
            x10  += (*SFP10ptr++) * Y_tmp;
            #endif
            #if nISO>=12
            x11  += (*SFP11ptr++) * Y_tmp;
            #endif
            #if nISO>=13
            x12  += (*SFP12ptr++) * Y_tmp;
            #endif
            #if nISO>=14
            x13  += (*SFP13ptr++) * Y_tmp;
            #endif
            #if nISO>=15
            x14  += (*SFP14ptr++) * Y_tmp;
            #endif
            #if nISO>=16
            x15  += (*SFP15ptr++) * Y_tmp;
            #endif
            #if nISO>=17
            x16  += (*SFP16ptr++) * Y_tmp;
            #endif
            #if nISO>=18
            x17  += (*SFP17ptr++) * Y_tmp;
            #endif
            #if nISO>=19
            x18  += (*SFP18ptr++) * Y_tmp;
            #endif
            #if nISO>=20
            x19  += (*SFP19ptr++) * Y_tmp;
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
        #if nISO>=5
        (*x_Ptr4++) += x4;
        #endif
        #if nISO>=6
        (*x_Ptr5++) += x5;
        #endif
        #if nISO>=7
        (*x_Ptr6++) += x6;
        #endif
        #if nISO>=8
        (*x_Ptr7++) += x7;
        #endif
        #if nISO>=9
        (*x_Ptr8++) += x8;
        #endif
        #if nISO>=10
        (*x_Ptr9++) += x9;
        #endif
        #if nISO>=11
        (*x_Ptr10++) += x10;
        #endif
        #if nISO>=12
        (*x_Ptr11++) += x11;
        #endif
        #if nISO>=13
        (*x_Ptr12++) += x12;
        #endif
        #if nISO>=14
        (*x_Ptr13++) += x13;
        #endif
        #if nISO>=15
        (*x_Ptr14++) += x14;
        #endif
        #if nISO>=16
        (*x_Ptr15++) += x15;
        #endif
        #if nISO>=17
        (*x_Ptr16++) += x16;
        #endif
        #if nISO>=18
        (*x_Ptr17++) += x17;
        #endif
        #if nISO>=19
        (*x_Ptr18++) += x18;
        #endif
        #if nISO>=20
        (*x_Ptr19++) += x19;
        #endif
    }
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
    nE = _nE;
    nV = _nV;
    nS = _nS;
    ndirs = _ndirs;

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
    wmrSFP1 = wmrSFP0 + _ndirs*_nS;
    #if nIC>=3
    wmrSFP2 = wmrSFP1 + _ndirs*_nS;
    #if nIC>=4
    wmrSFP3 = wmrSFP2 + _ndirs*_nS;
    #if nIC>=5
    wmrSFP4 = wmrSFP3 + _ndirs*_nS;
    #if nIC>=6
    wmrSFP5 = wmrSFP4 + _ndirs*_nS;
    #if nIC>=7
    wmrSFP6 = wmrSFP5 + _ndirs*_nS;
    #if nIC>=8
    wmrSFP7 = wmrSFP6 + _ndirs*_nS;
    #if nIC>=9
    wmrSFP8 = wmrSFP7 + _ndirs*_nS;
    #if nIC>=10
    wmrSFP9 = wmrSFP8 + _ndirs*_nS;
    #if nIC>=11
    wmrSFP10 = wmrSFP9 + _ndirs*_nS;
    #if nIC>=12
    wmrSFP11 = wmrSFP10 + _ndirs*_nS;
    #if nIC>=13
    wmrSFP12 = wmrSFP11 + _ndirs*_nS;
    #if nIC>=14
    wmrSFP13 = wmrSFP12 + _ndirs*_nS;
    #if nIC>=15
    wmrSFP14 = wmrSFP13 + _ndirs*_nS;
    #if nIC>=16
    wmrSFP15 = wmrSFP14 + _ndirs*_nS;
    #if nIC>=17
    wmrSFP16 = wmrSFP15 + _ndirs*_nS;
    #if nIC>=18
    wmrSFP17 = wmrSFP16 + _ndirs*_nS;
    #if nIC>=19
    wmrSFP18 = wmrSFP17 + _ndirs*_nS;
    #if nIC>=20
    wmrSFP19 = wmrSFP18 + _ndirs*_nS;
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
    #if nEC>=1
    wmhSFP0 = _wmhSFP;
    #if nEC>=2
    wmhSFP1 = wmhSFP0 + _ndirs*_nS;
    #if nEC>=3
    wmhSFP2 = wmhSFP1 + _ndirs*_nS;
    #if nEC>=4
    wmhSFP3 = wmhSFP2 + _ndirs*_nS;
    #if nEC>=5
    wmhSFP4 = wmhSFP3 + _ndirs*_nS;
    #if nEC>=6
    wmhSFP5 = wmhSFP4 + _ndirs*_nS;
    #if nEC>=7
    wmhSFP6 = wmhSFP5 + _ndirs*_nS;
    #if nEC>=8
    wmhSFP7 = wmhSFP6 + _ndirs*_nS;
    #if nEC>=9
    wmhSFP8 = wmhSFP7 + _ndirs*_nS;
    #if nEC>=10
    wmhSFP9 = wmhSFP8 + _ndirs*_nS;
    #if nEC>=11
    wmhSFP10 = wmhSFP9 + _ndirs*_nS;
    #if nEC>=12
    wmhSFP11 = wmhSFP10 + _ndirs*_nS;
    #if nEC>=13
    wmhSFP12 = wmhSFP11 + _ndirs*_nS;
    #if nEC>=14
    wmhSFP13 = wmhSFP12 + _ndirs*_nS;
    #if nEC>=15
    wmhSFP14 = wmhSFP13 + _ndirs*_nS;
    #if nEC>=16
    wmhSFP15 = wmhSFP14 + _ndirs*_nS;
    #if nEC>=17
    wmhSFP16 = wmhSFP15 + _ndirs*_nS;
    #if nEC>=18
    wmhSFP17 = wmhSFP16 + _ndirs*_nS;
    #if nEC>=19
    wmhSFP18 = wmhSFP17 + _ndirs*_nS;
    #if nEC>=20
    wmhSFP19 = wmhSFP18 + _ndirs*_nS;
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
    #if nISO>=1
    isoSFP0 = _isoSFP;
    #if nISO>=2
    isoSFP1 = isoSFP0 + _nS;
    #if nISO>=3
    isoSFP2 = isoSFP1 + _nS;
    #if nISO>=4
    isoSFP3 = isoSFP2 + _nS;
    #if nISO>=5
    isoSFP4 = isoSFP3 + _nS;
    #if nISO>=6
    isoSFP5 = isoSFP4 + _nS;
    #if nISO>=7
    isoSFP6 = isoSFP5 + _nS;
    #if nISO>=8
    isoSFP7 = isoSFP6 + _nS;
    #if nISO>=9
    isoSFP8 = isoSFP7 + _nS;
    #if nISO>=10
    isoSFP9 = isoSFP8 + _nS;
    #if nISO>=11
    isoSFP10 = isoSFP9 + _nS;
    #if nISO>=12
    isoSFP11 = isoSFP10 + _nS;
    #if nISO>=13
    isoSFP12 = isoSFP11 + _nS;
    #if nISO>=14
    isoSFP13 = isoSFP12 + _nS;
    #if nISO>=15
    isoSFP14 = isoSFP13 + _nS;
    #if nISO>=16
    isoSFP15 = isoSFP14 + _nS;
    #if nISO>=17
    isoSFP16 = isoSFP15 + _nS;
    #if nISO>=18
    isoSFP17 = isoSFP16 + _nS;
    #if nISO>=19
    isoSFP18 = isoSFP17 + _nS;
    #if nISO>=20
    isoSFP19 = isoSFP18 + _nS;
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
