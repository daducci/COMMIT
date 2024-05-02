#include <pthread.h>
#include <stdint.h> // uint32_t etc

// max number of threads
#define MAX_THREADS 255

// global variables
int         nF, n, nE, nV, nS, ndirs;
double      *x, *Y;
uint32_t    *ICthreads, *ECthreads, *ISOthreads;
uint8_t     *ICthreadsT;
uint32_t    *ECthreadsT, *ISOthreadsT;
uint32_t    *ICf, *ICeval, *ICv, *ECv, *ISOv;
uint16_t    *ICo, *ECo;
float       *ICl;
float       *wmrSFP0, *wmrSFP1, *wmrSFP2, *wmrSFP3, *wmrSFP4, *wmrSFP5, *wmrSFP6, *wmrSFP7, *wmrSFP8, *wmrSFP9, *wmrSFP10, *wmrSFP11, *wmrSFP12, *wmrSFP13, *wmrSFP14, *wmrSFP15, *wmrSFP16, *wmrSFP17, *wmrSFP18, *wmrSFP19;
float       *wmhSFP0, *wmhSFP1, *wmhSFP2, *wmhSFP3, *wmhSFP4, *wmhSFP5, *wmhSFP6, *wmhSFP7, *wmhSFP8, *wmhSFP9, *wmhSFP10, *wmhSFP11, *wmhSFP12, *wmhSFP13, *wmhSFP14, *wmhSFP15, *wmhSFP16, *wmhSFP17, *wmhSFP18, *wmhSFP19;
float       *isoSFP0, *isoSFP1, *isoSFP2, *isoSFP3, *isoSFP4, *isoSFP5, *isoSFP6, *isoSFP7, *isoSFP8, *isoSFP9, *isoSFP10, *isoSFP11, *isoSFP12, *isoSFP13, *isoSFP14, *isoSFP15, *isoSFP16, *isoSFP17, *isoSFP18, *isoSFP19;
uint32_t    nIC, nEC, nISO;


//////////////////////////////////////////////////////////
// Compute a sub-block of the A*x MATRIX-VECTOR product //
//////////////////////////////////////////////////////////
void* COMMIT_A__block( void *ptr )
{
    int      id = (long)ptr;
    int      offset;
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, w;
    double   *xPtr0, *xPtr1, *xPtr2, *xPtr3, *xPtr4, *xPtr5, *xPtr6, *xPtr7, *xPtr8, *xPtr9, *xPtr10, *xPtr11, *xPtr12, *xPtr13, *xPtr14, *xPtr15, *xPtr16, *xPtr17, *xPtr18, *xPtr19;
    double   *YPtr, *YPtrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr, *SFP10ptr, *SFP11ptr, *SFP12ptr, *SFP13ptr, *SFP14ptr, *SFP15ptr, *SFP16ptr, *SFP17ptr, *SFP18ptr, *SFP19ptr;
    uint32_t *t_v, *t_vEnd, *t_f;
    uint16_t *t_o;
    float    *t_l;

    // intra-cellular compartments
    if (nIC > 0)
    {
        t_v = ICv + ICthreads[id];
        t_vEnd = ICv + ICthreads[id+1];
        t_o = ICo + ICthreads[id];
        t_l = ICl + ICthreads[id];
        t_f = ICf + ICthreads[id];
        switch (nIC)
        {
            case 1:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    if (x0 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 2:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    if (x0 != 0 || x1 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 3:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    if (x0 != 0 || x1 != 0 || x2 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 4:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 5:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 6:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 7:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 8:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 9:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 10:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 11:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 12:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 13:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    xPtr12 = xPtr11 + nF;
                    x12 = *xPtr12;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 14:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    xPtr12 = xPtr11 + nF;
                    x12 = *xPtr12;
                    xPtr13 = xPtr12 + nF;
                    x13 = *xPtr13;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 15:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    xPtr12 = xPtr11 + nF;
                    x12 = *xPtr12;
                    xPtr13 = xPtr12 + nF;
                    x13 = *xPtr13;
                    xPtr14 = xPtr13 + nF;
                    x14 = *xPtr14;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 16:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    xPtr12 = xPtr11 + nF;
                    x12 = *xPtr12;
                    xPtr13 = xPtr12 + nF;
                    x13 = *xPtr13;
                    xPtr14 = xPtr13 + nF;
                    x14 = *xPtr14;
                    xPtr15 = xPtr14 + nF;
                    x15 = *xPtr15;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 17:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    xPtr12 = xPtr11 + nF;
                    x12 = *xPtr12;
                    xPtr13 = xPtr12 + nF;
                    x13 = *xPtr13;
                    xPtr14 = xPtr13 + nF;
                    x14 = *xPtr14;
                    xPtr15 = xPtr14 + nF;
                    x15 = *xPtr15;
                    xPtr16 = xPtr15 + nF;
                    x16 = *xPtr16;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        SFP16ptr = wmrSFP16 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 18:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    xPtr12 = xPtr11 + nF;
                    x12 = *xPtr12;
                    xPtr13 = xPtr12 + nF;
                    x13 = *xPtr13;
                    xPtr14 = xPtr13 + nF;
                    x14 = *xPtr14;
                    xPtr15 = xPtr14 + nF;
                    x15 = *xPtr15;
                    xPtr16 = xPtr15 + nF;
                    x16 = *xPtr16;
                    xPtr17 = xPtr16 + nF;
                    x17 = *xPtr17;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        SFP16ptr = wmrSFP16 + offset;
                        SFP17ptr = wmrSFP17 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 19:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    xPtr12 = xPtr11 + nF;
                    x12 = *xPtr12;
                    xPtr13 = xPtr12 + nF;
                    x13 = *xPtr13;
                    xPtr14 = xPtr13 + nF;
                    x14 = *xPtr14;
                    xPtr15 = xPtr14 + nF;
                    x15 = *xPtr15;
                    xPtr16 = xPtr15 + nF;
                    x16 = *xPtr16;
                    xPtr17 = xPtr16 + nF;
                    x17 = *xPtr17;
                    xPtr18 = xPtr17 + nF;
                    x18 = *xPtr18;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0 || x18 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        SFP16ptr = wmrSFP16 + offset;
                        SFP17ptr = wmrSFP17 + offset;
                        SFP18ptr = wmrSFP18 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++) + x18 * (*SFP18ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
            case 20:
                while (t_v != t_vEnd)
                {
                    xPtr0 = x + (*t_f);
                    x0 = *xPtr0;
                    xPtr1 = xPtr0 + nF;
                    x1 = *xPtr1;
                    xPtr2 = xPtr1 + nF;
                    x2 = *xPtr2;
                    xPtr3 = xPtr2 + nF;
                    x3 = *xPtr3;
                    xPtr4 = xPtr3 + nF;
                    x4 = *xPtr4;
                    xPtr5 = xPtr4 + nF;
                    x5 = *xPtr5;
                    xPtr6 = xPtr5 + nF;
                    x6 = *xPtr6;
                    xPtr7 = xPtr6 + nF;
                    x7 = *xPtr7;
                    xPtr8 = xPtr7 + nF;
                    x8 = *xPtr8;
                    xPtr9 = xPtr8 + nF;
                    x9 = *xPtr9;
                    xPtr10 = xPtr9 + nF;
                    x10 = *xPtr10;
                    xPtr11 = xPtr10 + nF;
                    x11 = *xPtr11;
                    xPtr12 = xPtr11 + nF;
                    x12 = *xPtr12;
                    xPtr13 = xPtr12 + nF;
                    x13 = *xPtr13;
                    xPtr14 = xPtr13 + nF;
                    x14 = *xPtr14;
                    xPtr15 = xPtr14 + nF;
                    x15 = *xPtr15;
                    xPtr16 = xPtr15 + nF;
                    x16 = *xPtr16;
                    xPtr17 = xPtr16 + nF;
                    x17 = *xPtr17;
                    xPtr18 = xPtr17 + nF;
                    x18 = *xPtr18;
                    xPtr19 = xPtr18 + nF;
                    x19 = *xPtr19;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0 || x18 != 0 || x19 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        SFP0ptr = wmrSFP0 + offset;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        SFP8ptr = wmrSFP8 + offset;
                        SFP9ptr = wmrSFP9 + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        SFP16ptr = wmrSFP16 + offset;
                        SFP17ptr = wmrSFP17 + offset;
                        SFP18ptr = wmrSFP18 + offset;
                        SFP19ptr = wmrSFP19 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++) + x18 * (*SFP18ptr++) + x19 * (*SFP19ptr++));
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                }
                break;
        }
    }

    // extra-cellular compartments
    if (nEC > 0)
    {
        t_v = ECv + ECthreads[id];
        t_vEnd = ECv + ECthreads[id+1];
        t_o = ECo + ECthreads[id];
        xPtr0 = x + nIC*nF + ECthreads[id];
        switch (nEC)
        {
            case 1:
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    if (x0 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 2:
                xPtr1 = xPtr0 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    if (x0 != 0 || x1 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 3:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    if (x0 != 0 || x1 != 0 || x2 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 4:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 5:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 6:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 7:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 8:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 9:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 10:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 11:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 12:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 13:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        SFP12ptr = wmhSFP12 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 14:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        SFP12ptr = wmhSFP12 + offset;
                        SFP13ptr = wmhSFP13 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 15:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        SFP12ptr = wmhSFP12 + offset;
                        SFP13ptr = wmhSFP13 + offset;
                        SFP14ptr = wmhSFP14 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 16:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        SFP12ptr = wmhSFP12 + offset;
                        SFP13ptr = wmhSFP13 + offset;
                        SFP14ptr = wmhSFP14 + offset;
                        SFP15ptr = wmhSFP15 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 17:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                xPtr16 = xPtr15 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    x16 = *xPtr16++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        SFP12ptr = wmhSFP12 + offset;
                        SFP13ptr = wmhSFP13 + offset;
                        SFP14ptr = wmhSFP14 + offset;
                        SFP15ptr = wmhSFP15 + offset;
                        SFP16ptr = wmhSFP16 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 18:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                xPtr16 = xPtr15 + nE;
                xPtr17 = xPtr16 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    x16 = *xPtr16++;
                    x17 = *xPtr17++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        SFP12ptr = wmhSFP12 + offset;
                        SFP13ptr = wmhSFP13 + offset;
                        SFP14ptr = wmhSFP14 + offset;
                        SFP15ptr = wmhSFP15 + offset;
                        SFP16ptr = wmhSFP16 + offset;
                        SFP17ptr = wmhSFP17 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 19:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                xPtr16 = xPtr15 + nE;
                xPtr17 = xPtr16 + nE;
                xPtr18 = xPtr17 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    x16 = *xPtr16++;
                    x17 = *xPtr17++;
                    x18 = *xPtr18++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0 || x18 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        SFP12ptr = wmhSFP12 + offset;
                        SFP13ptr = wmhSFP13 + offset;
                        SFP14ptr = wmhSFP14 + offset;
                        SFP15ptr = wmhSFP15 + offset;
                        SFP16ptr = wmhSFP16 + offset;
                        SFP17ptr = wmhSFP17 + offset;
                        SFP18ptr = wmhSFP18 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++) + x18 * (*SFP18ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
            case 20:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                xPtr16 = xPtr15 + nE;
                xPtr17 = xPtr16 + nE;
                xPtr18 = xPtr17 + nE;
                xPtr19 = xPtr18 + nE;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    x16 = *xPtr16++;
                    x17 = *xPtr17++;
                    x18 = *xPtr18++;
                    x19 = *xPtr19++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0 || x18 != 0 || x19 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        SFP0ptr = wmhSFP0 + offset;
                        SFP1ptr = wmhSFP1 + offset;
                        SFP2ptr = wmhSFP2 + offset;
                        SFP3ptr = wmhSFP3 + offset;
                        SFP4ptr = wmhSFP4 + offset;
                        SFP5ptr = wmhSFP5 + offset;
                        SFP6ptr = wmhSFP6 + offset;
                        SFP7ptr = wmhSFP7 + offset;
                        SFP8ptr = wmhSFP8 + offset;
                        SFP9ptr = wmhSFP9 + offset;
                        SFP10ptr = wmhSFP10 + offset;
                        SFP11ptr = wmhSFP11 + offset;
                        SFP12ptr = wmhSFP12 + offset;
                        SFP13ptr = wmhSFP13 + offset;
                        SFP14ptr = wmhSFP14 + offset;
                        SFP15ptr = wmhSFP15 + offset;
                        SFP16ptr = wmhSFP16 + offset;
                        SFP17ptr = wmhSFP17 + offset;
                        SFP18ptr = wmhSFP18 + offset;
                        SFP19ptr = wmhSFP19 + offset;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++) + x18 * (*SFP18ptr++) + x19 * (*SFP19ptr++));
                    }
                    t_v++;
                    t_o++;
                }
                break;
        }
    }

    // isotropic compartments
    if (nISO > 0)
    {
        t_v = ISOv + ISOthreads[id];
        t_vEnd = ISOv + ISOthreads[id+1];
        xPtr0 = x + nIC*nF + nEC*nE + ISOthreads[id];
        switch (nISO)
        {
            case 1:
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    if (x0 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++));
                    }
                    t_v++;
                }
                break;
            case 2:
                xPtr1 = xPtr0 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    if (x0 != 0 || x1 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++));
                    }
                    t_v++;
                }
                break;
            case 3:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    if (x0 != 0 || x1 != 0 || x2 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++));
                    }
                    t_v++;
                }
                break;
            case 4:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++));
                    }
                    t_v++;
                }
                break;
            case 5:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++));
                    }
                    t_v++;
                }
                break;
            case 6:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++));
                    }
                    t_v++;
                }
                break;
            case 7:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++));
                    }
                    t_v++;
                }
                break;
            case 8:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++));
                    }
                    t_v++;
                }
                break;
            case 9:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++));
                    }
                    t_v++;
                }
                break;
            case 10:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++));
                    }
                    t_v++;
                }
                break;
            case 11:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++));
                    }
                    t_v++;
                }
                break;
            case 12:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++));
                    }
                    t_v++;
                }
                break;
            case 13:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        SFP12ptr = isoSFP12;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++));
                    }
                    t_v++;
                }
                break;
            case 14:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        SFP12ptr = isoSFP12;
                        SFP13ptr = isoSFP13;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++));
                    }
                    t_v++;
                }
                break;
            case 15:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        SFP12ptr = isoSFP12;
                        SFP13ptr = isoSFP13;
                        SFP14ptr = isoSFP14;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++));
                    }
                    t_v++;
                }
                break;
            case 16:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        SFP12ptr = isoSFP12;
                        SFP13ptr = isoSFP13;
                        SFP14ptr = isoSFP14;
                        SFP15ptr = isoSFP15;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++));
                    }
                    t_v++;
                }
                break;
            case 17:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                xPtr16 = xPtr15 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    x16 = *xPtr16++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        SFP12ptr = isoSFP12;
                        SFP13ptr = isoSFP13;
                        SFP14ptr = isoSFP14;
                        SFP15ptr = isoSFP15;
                        SFP16ptr = isoSFP16;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++));
                    }
                    t_v++;
                }
                break;
            case 18:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                xPtr16 = xPtr15 + nV;
                xPtr17 = xPtr16 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    x16 = *xPtr16++;
                    x17 = *xPtr17++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        SFP12ptr = isoSFP12;
                        SFP13ptr = isoSFP13;
                        SFP14ptr = isoSFP14;
                        SFP15ptr = isoSFP15;
                        SFP16ptr = isoSFP16;
                        SFP17ptr = isoSFP17;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++));
                    }
                    t_v++;
                }
                break;
            case 19:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                xPtr16 = xPtr15 + nV;
                xPtr17 = xPtr16 + nV;
                xPtr18 = xPtr17 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    x16 = *xPtr16++;
                    x17 = *xPtr17++;
                    x18 = *xPtr18++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0 || x18 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        SFP12ptr = isoSFP12;
                        SFP13ptr = isoSFP13;
                        SFP14ptr = isoSFP14;
                        SFP15ptr = isoSFP15;
                        SFP16ptr = isoSFP16;
                        SFP17ptr = isoSFP17;
                        SFP18ptr = isoSFP18;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++) + x18 * (*SFP18ptr++));
                    }
                    t_v++;
                }
                break;
            case 20:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                xPtr16 = xPtr15 + nV;
                xPtr17 = xPtr16 + nV;
                xPtr18 = xPtr17 + nV;
                xPtr19 = xPtr18 + nV;
                while (t_v != t_vEnd)
                {
                    x0 = *xPtr0++;
                    x1 = *xPtr1++;
                    x2 = *xPtr2++;
                    x3 = *xPtr3++;
                    x4 = *xPtr4++;
                    x5 = *xPtr5++;
                    x6 = *xPtr6++;
                    x7 = *xPtr7++;
                    x8 = *xPtr8++;
                    x9 = *xPtr9++;
                    x10 = *xPtr10++;
                    x11 = *xPtr11++;
                    x12 = *xPtr12++;
                    x13 = *xPtr13++;
                    x14 = *xPtr14++;
                    x15 = *xPtr15++;
                    x16 = *xPtr16++;
                    x17 = *xPtr17++;
                    x18 = *xPtr18++;
                    x19 = *xPtr19++;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0 || x10 != 0 || x11 != 0 || x12 != 0 || x13 != 0 || x14 != 0 || x15 != 0 || x16 != 0 || x17 != 0 || x18 != 0 || x19 != 0)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        SFP0ptr = isoSFP0;
                        SFP1ptr = isoSFP1;
                        SFP2ptr = isoSFP2;
                        SFP3ptr = isoSFP3;
                        SFP4ptr = isoSFP4;
                        SFP5ptr = isoSFP5;
                        SFP6ptr = isoSFP6;
                        SFP7ptr = isoSFP7;
                        SFP8ptr = isoSFP8;
                        SFP9ptr = isoSFP9;
                        SFP10ptr = isoSFP10;
                        SFP11ptr = isoSFP11;
                        SFP12ptr = isoSFP12;
                        SFP13ptr = isoSFP13;
                        SFP14ptr = isoSFP14;
                        SFP15ptr = isoSFP15;
                        SFP16ptr = isoSFP16;
                        SFP17ptr = isoSFP17;
                        SFP18ptr = isoSFP18;
                        SFP19ptr = isoSFP19;
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += (x0 * (*SFP0ptr++) + x1 * (*SFP1ptr++) + x2 * (*SFP2ptr++) + x3 * (*SFP3ptr++) + x4 * (*SFP4ptr++) + x5 * (*SFP5ptr++) + x6 * (*SFP6ptr++) + x7 * (*SFP7ptr++) + x8 * (*SFP8ptr++) + x9 * (*SFP9ptr++) + x10 * (*SFP10ptr++) + x11 * (*SFP11ptr++) + x12 * (*SFP12ptr++) + x13 * (*SFP13ptr++) + x14 * (*SFP14ptr++) + x15 * (*SFP15ptr++) + x16 * (*SFP16ptr++) + x17 * (*SFP17ptr++) + x18 * (*SFP18ptr++) + x19 * (*SFP19ptr++));
                    }
                    t_v++;
                }
                break;
        }
    }

    pthread_exit( 0 );
}

///////////////////////////////
// Function called by CYTHON //
///////////////////////////////
void COMMIT_A(
    int _nF, int _n, int _nE, int _nV, int _nS, int _ndirs,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICeval, uint32_t *_ICv, uint16_t *_ICo, float *_ICl,
    uint32_t *_ECv, uint16_t *_ECo,
    uint32_t *_ISOv,
    float *_wmrSFP, float *_wmhSFP, float *_isoSFP,
    uint32_t* _ICthreads, uint32_t* _ECthreads, uint32_t* _ISOthreads,
    uint32_t _nIC, uint32_t _nEC, uint32_t _nISO, uint32_t _nThreads
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
    ICeval = _ICeval;
    ICv  = _ICv;
    ICo  = _ICo;
    ICl  = _ICl;
    ECv  = _ECv;
    ECo  = _ECo;
    ISOv = _ISOv;

    nIC = _nIC;
    nEC = _nEC;
    nISO = _nISO;

    switch (nIC)
    {
        case 1:
            wmrSFP0 = _wmrSFP;
            break;
        case 2:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            break;
        case 3:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            break;
        case 4:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            break;
        case 5:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            break;
        case 6:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            wmrSFP5 = wmrSFP4 + _ndirs*_nS;
            break;
        case 7:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            wmrSFP5 = wmrSFP4 + _ndirs*_nS;
            wmrSFP6 = wmrSFP5 + _ndirs*_nS;
            break;
        case 8:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            wmrSFP5 = wmrSFP4 + _ndirs*_nS;
            wmrSFP6 = wmrSFP5 + _ndirs*_nS;
            wmrSFP7 = wmrSFP6 + _ndirs*_nS;
            break;
        case 9:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            wmrSFP5 = wmrSFP4 + _ndirs*_nS;
            wmrSFP6 = wmrSFP5 + _ndirs*_nS;
            wmrSFP7 = wmrSFP6 + _ndirs*_nS;
            wmrSFP8 = wmrSFP7 + _ndirs*_nS;
            break;
        case 10:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            break;
        case 11:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            break;
        case 12:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            break;
        case 13:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            break;
        case 14:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            break;
        case 15:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            break;
        case 16:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            break;
        case 17:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            wmrSFP16 = wmrSFP15 + _ndirs*_nS;
            break;
        case 18:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            wmrSFP16 = wmrSFP15 + _ndirs*_nS;
            wmrSFP17 = wmrSFP16 + _ndirs*_nS;
            break;
        case 19:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0  + _ndirs*_nS;
            wmrSFP2  = wmrSFP1  + _ndirs*_nS;
            wmrSFP3  = wmrSFP2  + _ndirs*_nS;
            wmrSFP4  = wmrSFP3  + _ndirs*_nS;
            wmrSFP5  = wmrSFP4  + _ndirs*_nS;
            wmrSFP6  = wmrSFP5  + _ndirs*_nS;
            wmrSFP7  = wmrSFP6  + _ndirs*_nS;
            wmrSFP8  = wmrSFP7  + _ndirs*_nS;
            wmrSFP9  = wmrSFP8  + _ndirs*_nS;
            wmrSFP10 = wmrSFP9  + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            wmrSFP16 = wmrSFP15 + _ndirs*_nS;
            wmrSFP17 = wmrSFP16 + _ndirs*_nS;
            wmrSFP18 = wmrSFP17 + _ndirs*_nS;
            break;
        case 20:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0  + _ndirs*_nS;
            wmrSFP2  = wmrSFP1  + _ndirs*_nS;
            wmrSFP3  = wmrSFP2  + _ndirs*_nS;
            wmrSFP4  = wmrSFP3  + _ndirs*_nS;
            wmrSFP5  = wmrSFP4  + _ndirs*_nS;
            wmrSFP6  = wmrSFP5  + _ndirs*_nS;
            wmrSFP7  = wmrSFP6  + _ndirs*_nS;
            wmrSFP8  = wmrSFP7  + _ndirs*_nS;
            wmrSFP9  = wmrSFP8  + _ndirs*_nS;
            wmrSFP10 = wmrSFP9  + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            wmrSFP16 = wmrSFP15 + _ndirs*_nS;
            wmrSFP17 = wmrSFP16 + _ndirs*_nS;
            wmrSFP18 = wmrSFP17 + _ndirs*_nS;
            wmrSFP19 = wmrSFP18 + _ndirs*_nS;
            break;
    }

    switch (nEC)
    {
        case 1:
            wmhSFP0 = _wmhSFP;
            break;
        case 2:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            break;
        case 3:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            break;
        case 4:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            break;
        case 5:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            break;
        case 6:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            wmhSFP5 = wmhSFP4 + _ndirs*_nS;
            break;
        case 7:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            wmhSFP5 = wmhSFP4 + _ndirs*_nS;
            wmhSFP6 = wmhSFP5 + _ndirs*_nS;
            break;
        case 8:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            wmhSFP5 = wmhSFP4 + _ndirs*_nS;
            wmhSFP6 = wmhSFP5 + _ndirs*_nS;
            wmhSFP7 = wmhSFP6 + _ndirs*_nS;
            break;
        case 9:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            wmhSFP5 = wmhSFP4 + _ndirs*_nS;
            wmhSFP6 = wmhSFP5 + _ndirs*_nS;
            wmhSFP7 = wmhSFP6 + _ndirs*_nS;
            wmhSFP8 = wmhSFP7 + _ndirs*_nS;
            break;
        case 10:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            break;
        case 11:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            break;
        case 12:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            break;
        case 13:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            break;
        case 14:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            break;
        case 15:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            break;
        case 16:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            break;
        case 17:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            wmhSFP16 = wmhSFP15 + _ndirs*_nS;
            break;
        case 18:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            wmhSFP16 = wmhSFP15 + _ndirs*_nS;
            wmhSFP17 = wmhSFP16 + _ndirs*_nS;
            break;
        case 19:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0  + _ndirs*_nS;
            wmhSFP2  = wmhSFP1  + _ndirs*_nS;
            wmhSFP3  = wmhSFP2  + _ndirs*_nS;
            wmhSFP4  = wmhSFP3  + _ndirs*_nS;
            wmhSFP5  = wmhSFP4  + _ndirs*_nS;
            wmhSFP6  = wmhSFP5  + _ndirs*_nS;
            wmhSFP7  = wmhSFP6  + _ndirs*_nS;
            wmhSFP8  = wmhSFP7  + _ndirs*_nS;
            wmhSFP9  = wmhSFP8  + _ndirs*_nS;
            wmhSFP10 = wmhSFP9  + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            wmhSFP16 = wmhSFP15 + _ndirs*_nS;
            wmhSFP17 = wmhSFP16 + _ndirs*_nS;
            wmhSFP18 = wmhSFP17 + _ndirs*_nS;
            break;
        case 20:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0  + _ndirs*_nS;
            wmhSFP2  = wmhSFP1  + _ndirs*_nS;
            wmhSFP3  = wmhSFP2  + _ndirs*_nS;
            wmhSFP4  = wmhSFP3  + _ndirs*_nS;
            wmhSFP5  = wmhSFP4  + _ndirs*_nS;
            wmhSFP6  = wmhSFP5  + _ndirs*_nS;
            wmhSFP7  = wmhSFP6  + _ndirs*_nS;
            wmhSFP8  = wmhSFP7  + _ndirs*_nS;
            wmhSFP9  = wmhSFP8  + _ndirs*_nS;
            wmhSFP10 = wmhSFP9  + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            wmhSFP16 = wmhSFP15 + _ndirs*_nS;
            wmhSFP17 = wmhSFP16 + _ndirs*_nS;
            wmhSFP18 = wmhSFP17 + _ndirs*_nS;
            wmhSFP19 = wmhSFP18 + _ndirs*_nS;
            break;
    }

    switch (nISO)
    {
        case 1:
            isoSFP0 = _isoSFP;
            break;
        case 2:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            break;
        case 3:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            break;
        case 4:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            break;
        case 5:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            break;
        case 6:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            isoSFP5 = isoSFP4 + _nS;
            break;
        case 7:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            isoSFP5 = isoSFP4 + _nS;
            isoSFP6 = isoSFP5 + _nS;
            break;
        case 8:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            isoSFP5 = isoSFP4 + _nS;
            isoSFP6 = isoSFP5 + _nS;
            isoSFP7 = isoSFP6 + _nS;
            break;
        case 9:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            isoSFP5 = isoSFP4 + _nS;
            isoSFP6 = isoSFP5 + _nS;
            isoSFP7 = isoSFP6 + _nS;
            isoSFP8 = isoSFP7 + _nS;
            break;
        case 10:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            break;
        case 11:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            break;
        case 12:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            break;
        case 13:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            break;
        case 14:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            break;
        case 15:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            break;
        case 16:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            break;
        case 17:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            isoSFP16 = isoSFP15 + _nS;
            break;
        case 18:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            isoSFP16 = isoSFP15 + _nS;
            isoSFP17 = isoSFP16 + _nS;
            break;
        case 19:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0  + _nS;
            isoSFP2  = isoSFP1  + _nS;
            isoSFP3  = isoSFP2  + _nS;
            isoSFP4  = isoSFP3  + _nS;
            isoSFP5  = isoSFP4  + _nS;
            isoSFP6  = isoSFP5  + _nS;
            isoSFP7  = isoSFP6  + _nS;
            isoSFP8  = isoSFP7  + _nS;
            isoSFP9  = isoSFP8  + _nS;
            isoSFP10 = isoSFP9  + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            isoSFP16 = isoSFP15 + _nS;
            isoSFP17 = isoSFP16 + _nS;
            isoSFP18 = isoSFP17 + _nS;
            break;
        case 20:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0  + _nS;
            isoSFP2  = isoSFP1  + _nS;
            isoSFP3  = isoSFP2  + _nS;
            isoSFP4  = isoSFP3  + _nS;
            isoSFP5  = isoSFP4  + _nS;
            isoSFP6  = isoSFP5  + _nS;
            isoSFP7  = isoSFP6  + _nS;
            isoSFP8  = isoSFP7  + _nS;
            isoSFP9  = isoSFP8  + _nS;
            isoSFP10 = isoSFP9  + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            isoSFP16 = isoSFP15 + _nS;
            isoSFP17 = isoSFP16 + _nS;
            isoSFP18 = isoSFP17 + _nS;
            isoSFP19 = isoSFP18 + _nS;
            break;
    }

    ICthreads  = _ICthreads;
    ECthreads  = _ECthreads;
    ISOthreads = _ISOthreads;

    // Run SEPARATE THREADS to perform the multiplication
    pthread_t threads[MAX_THREADS];
    int t;
    for(t=0; t<_nThreads ; t++)
        pthread_create( &threads[t], NULL, COMMIT_A__block, (void *) (long int)t );
    for(t=0; t<_nThreads ; t++)
        pthread_join( threads[t], NULL );
    return;
}


///////////////////////////////////////////////////////////
// Compute a sub-block of the A'*y MAtRIX-VECTOR product //
///////////////////////////////////////////////////////////
void* COMMIT_At__block( void *ptr )
{
    int      id = (long)ptr;
    int      offset;
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, w, YTmp;
    double   *xPtr0, *xPtr1, *xPtr2, *xPtr3, *xPtr4, *xPtr5, *xPtr6, *xPtr7, *xPtr8, *xPtr9, *xPtr10, *xPtr11, *xPtr12, *xPtr13, *xPtr14, *xPtr15, *xPtr16, *xPtr17, *xPtr18, *xPtr19;
    double   *YPtr, *YPtrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr, *SFP10ptr, *SFP11ptr, *SFP12ptr, *SFP13ptr, *SFP14ptr, *SFP15ptr, *SFP16ptr, *SFP17ptr, *SFP18ptr, *SFP19ptr;
    uint32_t *t_v, *t_vEnd, *t_f;
    uint16_t *t_o;
    float    *t_l;
    uint8_t  *t_t;

    // intra-cellular compartments
    if (nIC > 0)
    {
        t_v = ICv;
        t_vEnd = ICv + n;
        t_o = ICo;
        t_l = ICl;
        t_f = ICf;
        t_t = ICthreadsT;
        switch (nIC)
        {
            case 1:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr = wmrSFP0 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 2:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 3:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 4:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    t_t++;
                }
                break;

            case 5:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    t_t++;
                }
                break;

            case 6:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    t_t++;
                }
                break;

            case 7:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 8:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr = wmrSFP0 + offset;
                        SFP1ptr = wmrSFP1 + offset;
                        SFP2ptr = wmrSFP2 + offset;
                        SFP3ptr = wmrSFP3 + offset;
                        SFP4ptr = wmrSFP4 + offset;
                        SFP5ptr = wmrSFP5 + offset;
                        SFP6ptr = wmrSFP6 + offset;
                        SFP7ptr = wmrSFP7 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 9:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 10:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 11:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 12:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 13:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        x12 = (*SFP12ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                            x12 += (*SFP12ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                        x[*t_f+12*nF] += w * x12;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 14:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        x12 = (*SFP12ptr++) * YTmp;
                        x13 = (*SFP13ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                            x12 += (*SFP12ptr++) * YTmp;
                            x13 += (*SFP13ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                        x[*t_f+12*nF] += w * x12;
                        x[*t_f+13*nF] += w * x13;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 15:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        x12 = (*SFP12ptr++) * YTmp;
                        x13 = (*SFP13ptr++) * YTmp;
                        x14 = (*SFP14ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                            x12 += (*SFP12ptr++) * YTmp;
                            x13 += (*SFP13ptr++) * YTmp;
                            x14 += (*SFP14ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                        x[*t_f+12*nF] += w * x12;
                        x[*t_f+13*nF] += w * x13;
                        x[*t_f+14*nF] += w * x14;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 16:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        x12 = (*SFP12ptr++) * YTmp;
                        x13 = (*SFP13ptr++) * YTmp;
                        x14 = (*SFP14ptr++) * YTmp;
                        x15 = (*SFP15ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                            x12 += (*SFP12ptr++) * YTmp;
                            x13 += (*SFP13ptr++) * YTmp;
                            x14 += (*SFP14ptr++) * YTmp;
                            x15 += (*SFP15ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                        x[*t_f+12*nF] += w * x12;
                        x[*t_f+13*nF] += w * x13;
                        x[*t_f+14*nF] += w * x14;
                        x[*t_f+15*nF] += w * x15;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 17:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        SFP16ptr = wmrSFP16 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        x12 = (*SFP12ptr++) * YTmp;
                        x13 = (*SFP13ptr++) * YTmp;
                        x14 = (*SFP14ptr++) * YTmp;
                        x15 = (*SFP15ptr++) * YTmp;
                        x16 = (*SFP16ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                            x12 += (*SFP12ptr++) * YTmp;
                            x13 += (*SFP13ptr++) * YTmp;
                            x14 += (*SFP14ptr++) * YTmp;
                            x15 += (*SFP15ptr++) * YTmp;
                            x16 += (*SFP16ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                        x[*t_f+12*nF] += w * x12;
                        x[*t_f+13*nF] += w * x13;
                        x[*t_f+14*nF] += w * x14;
                        x[*t_f+15*nF] += w * x15;
                        x[*t_f+16*nF] += w * x16;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 18:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        SFP16ptr = wmrSFP16 + offset;
                        SFP17ptr = wmrSFP17 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        x12 = (*SFP12ptr++) * YTmp;
                        x13 = (*SFP13ptr++) * YTmp;
                        x14 = (*SFP14ptr++) * YTmp;
                        x15 = (*SFP15ptr++) * YTmp;
                        x16 = (*SFP16ptr++) * YTmp;
                        x17 = (*SFP17ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                            x12 += (*SFP12ptr++) * YTmp;
                            x13 += (*SFP13ptr++) * YTmp;
                            x14 += (*SFP14ptr++) * YTmp;
                            x15 += (*SFP15ptr++) * YTmp;
                            x16 += (*SFP16ptr++) * YTmp;
                            x17 += (*SFP17ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                        x[*t_f+12*nF] += w * x12;
                        x[*t_f+13*nF] += w * x13;
                        x[*t_f+14*nF] += w * x14;
                        x[*t_f+15*nF] += w * x15;
                        x[*t_f+16*nF] += w * x16;
                        x[*t_f+17*nF] += w * x17;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 19:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        SFP16ptr = wmrSFP16 + offset;
                        SFP17ptr = wmrSFP17 + offset;
                        SFP18ptr = wmrSFP18 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        x12 = (*SFP12ptr++) * YTmp;
                        x13 = (*SFP13ptr++) * YTmp;
                        x14 = (*SFP14ptr++) * YTmp;
                        x15 = (*SFP15ptr++) * YTmp;
                        x16 = (*SFP16ptr++) * YTmp;
                        x17 = (*SFP17ptr++) * YTmp;
                        x18 = (*SFP18ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                            x12 += (*SFP12ptr++) * YTmp;
                            x13 += (*SFP13ptr++) * YTmp;
                            x14 += (*SFP14ptr++) * YTmp;
                            x15 += (*SFP15ptr++) * YTmp;
                            x16 += (*SFP16ptr++) * YTmp;
                            x17 += (*SFP17ptr++) * YTmp;
                            x18 += (*SFP18ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                        x[*t_f+12*nF] += w * x12;
                        x[*t_f+13*nF] += w * x13;
                        x[*t_f+14*nF] += w * x14;
                        x[*t_f+15*nF] += w * x15;
                        x[*t_f+16*nF] += w * x16;
                        x[*t_f+17*nF] += w * x17;
                        x[*t_f+18*nF] += w * x18;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;

            case 20:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        SFP0ptr  = wmrSFP0  + offset;
                        SFP1ptr  = wmrSFP1  + offset;
                        SFP2ptr  = wmrSFP2  + offset;
                        SFP3ptr  = wmrSFP3  + offset;
                        SFP4ptr  = wmrSFP4  + offset;
                        SFP5ptr  = wmrSFP5  + offset;
                        SFP6ptr  = wmrSFP6  + offset;
                        SFP7ptr  = wmrSFP7  + offset;
                        SFP8ptr  = wmrSFP8  + offset;
                        SFP9ptr  = wmrSFP9  + offset;
                        SFP10ptr = wmrSFP10 + offset;
                        SFP11ptr = wmrSFP11 + offset;
                        SFP12ptr = wmrSFP12 + offset;
                        SFP13ptr = wmrSFP13 + offset;
                        SFP14ptr = wmrSFP14 + offset;
                        SFP15ptr = wmrSFP15 + offset;
                        SFP16ptr = wmrSFP16 + offset;
                        SFP17ptr = wmrSFP17 + offset;
                        SFP18ptr = wmrSFP18 + offset;
                        SFP19ptr = wmrSFP19 + offset;
                        x0 = (*SFP0ptr++) * YTmp;
                        x1 = (*SFP1ptr++) * YTmp;
                        x2 = (*SFP2ptr++) * YTmp;
                        x3 = (*SFP3ptr++) * YTmp;
                        x4 = (*SFP4ptr++) * YTmp;
                        x5 = (*SFP5ptr++) * YTmp;
                        x6 = (*SFP6ptr++) * YTmp;
                        x7 = (*SFP7ptr++) * YTmp;
                        x8 = (*SFP8ptr++) * YTmp;
                        x9 = (*SFP9ptr++) * YTmp;
                        x10 = (*SFP10ptr++) * YTmp;
                        x11 = (*SFP11ptr++) * YTmp;
                        x12 = (*SFP12ptr++) * YTmp;
                        x13 = (*SFP13ptr++) * YTmp;
                        x14 = (*SFP14ptr++) * YTmp;
                        x15 = (*SFP15ptr++) * YTmp;
                        x16 = (*SFP16ptr++) * YTmp;
                        x17 = (*SFP17ptr++) * YTmp;
                        x18 = (*SFP18ptr++) * YTmp;
                        x19 = (*SFP19ptr++) * YTmp;
                        while (++YPtr != YPtrEnd)
                        {
                            YTmp = *YPtr;
                            x0 += (*SFP0ptr++) * YTmp;
                            x1 += (*SFP1ptr++) * YTmp;
                            x2 += (*SFP2ptr++) * YTmp;
                            x3 += (*SFP3ptr++) * YTmp;
                            x4 += (*SFP4ptr++) * YTmp;
                            x5 += (*SFP5ptr++) * YTmp;
                            x6 += (*SFP6ptr++) * YTmp;
                            x7 += (*SFP7ptr++) * YTmp;
                            x8 += (*SFP8ptr++) * YTmp;
                            x9 += (*SFP9ptr++) * YTmp;
                            x10 += (*SFP10ptr++) * YTmp;
                            x11 += (*SFP11ptr++) * YTmp;
                            x12 += (*SFP12ptr++) * YTmp;
                            x13 += (*SFP13ptr++) * YTmp;
                            x14 += (*SFP14ptr++) * YTmp;
                            x15 += (*SFP15ptr++) * YTmp;
                            x16 += (*SFP16ptr++) * YTmp;
                            x17 += (*SFP17ptr++) * YTmp;
                            x18 += (*SFP18ptr++) * YTmp;
                            x19 += (*SFP19ptr++) * YTmp;
                        }
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                        x[*t_f+10*nF] += w * x10;
                        x[*t_f+11*nF] += w * x11;
                        x[*t_f+12*nF] += w * x12;
                        x[*t_f+13*nF] += w * x13;
                        x[*t_f+14*nF] += w * x14;
                        x[*t_f+15*nF] += w * x15;
                        x[*t_f+16*nF] += w * x16;
                        x[*t_f+17*nF] += w * x17;
                        x[*t_f+18*nF] += w * x18;
                        x[*t_f+19*nF] += w * x19;
                    }
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    
                    t_t++;
                }
                break;
        }
    }

    // extra-cellular compartments
    if (nEC > 0)
    {
        t_v = ECv + ECthreadsT[id];
        t_vEnd = ECv + ECthreadsT[id+1];
        t_o = ECo + ECthreadsT[id];
        xPtr0 = x + nIC*nF + ECthreadsT[id];
        switch (nEC)
        {
            case 1:
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    t_v++;
                    t_o++;
                }
                break;
            case 2:
                xPtr1 = xPtr0 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    t_v++;
                    t_o++;
                }
                break;
            case 3:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    t_v++;
                    t_o++;
                }
                break;

            case 4:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    t_v++;
                    t_o++;
                }
                break;

            case 5:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    t_v++;
                    t_o++;
                }
                break;

            case 6:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    t_v++;
                    t_o++;
                }
                break;

            case 7:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    t_v++;
                    t_o++;
                }
                break;

            case 8:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    t_v++;
                    t_o++;
                }
                break;

            case 9:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    t_v++;
                    t_o++;
                }
                break;

            case 10:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    t_v++;
                    t_o++;
                }
                break;

            case 11:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    t_v++;
                    t_o++;
                }
                break;

            case 12:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    t_v++;
                    t_o++;
                }
                break;

            case 13:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    SFP12ptr = wmhSFP12 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    t_v++;
                    t_o++;
                }
                break;

            case 14:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    SFP12ptr = wmhSFP12 + offset;
                    SFP13ptr = wmhSFP13 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    t_v++;
                    t_o++;
                }
                break;

            case 15:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    SFP12ptr = wmhSFP12 + offset;
                    SFP13ptr = wmhSFP13 + offset;
                    SFP14ptr = wmhSFP14 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    t_v++;
                    t_o++;
                }
                break;

            case 16:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    SFP12ptr = wmhSFP12 + offset;
                    SFP13ptr = wmhSFP13 + offset;
                    SFP14ptr = wmhSFP14 + offset;
                    SFP15ptr = wmhSFP15 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    t_v++;
                    t_o++;
                }
                break;

            case 17:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                xPtr16 = xPtr15 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    SFP12ptr = wmhSFP12 + offset;
                    SFP13ptr = wmhSFP13 + offset;
                    SFP14ptr = wmhSFP14 + offset;
                    SFP15ptr = wmhSFP15 + offset;
                    SFP16ptr = wmhSFP16 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    x16 = (*SFP16ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                        x16 += (*SFP16ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    (*xPtr16++) += x16;
                    t_v++;
                    t_o++;
                }
                break;

            case 18:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                xPtr16 = xPtr15 + nE;
                xPtr17 = xPtr16 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    SFP12ptr = wmhSFP12 + offset;
                    SFP13ptr = wmhSFP13 + offset;
                    SFP14ptr = wmhSFP14 + offset;
                    SFP15ptr = wmhSFP15 + offset;
                    SFP16ptr = wmhSFP16 + offset;
                    SFP17ptr = wmhSFP17 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    x16 = (*SFP16ptr++) * YTmp;
                    x17 = (*SFP17ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                        x16 += (*SFP16ptr++) * YTmp;
                        x17 += (*SFP17ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    (*xPtr16++) += x16;
                    (*xPtr17++) += x17;
                    t_v++;
                    t_o++;
                }
                break;

            case 19:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                xPtr16 = xPtr15 + nE;
                xPtr17 = xPtr16 + nE;
                xPtr18 = xPtr17 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    SFP12ptr = wmhSFP12 + offset;
                    SFP13ptr = wmhSFP13 + offset;
                    SFP14ptr = wmhSFP14 + offset;
                    SFP15ptr = wmhSFP15 + offset;
                    SFP16ptr = wmhSFP16 + offset;
                    SFP17ptr = wmhSFP17 + offset;
                    SFP18ptr = wmhSFP18 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    x16 = (*SFP16ptr++) * YTmp;
                    x17 = (*SFP17ptr++) * YTmp;
                    x18 = (*SFP18ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                        x16 += (*SFP16ptr++) * YTmp;
                        x17 += (*SFP17ptr++) * YTmp;
                        x18 += (*SFP18ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    (*xPtr16++) += x16;
                    (*xPtr17++) += x17;
                    (*xPtr18++) += x18;
                    t_v++;
                    t_o++;
                }
                break;

            case 20:
                xPtr1 = xPtr0 + nE;
                xPtr2 = xPtr1 + nE;
                xPtr3 = xPtr2 + nE;
                xPtr4 = xPtr3 + nE;
                xPtr5 = xPtr4 + nE;
                xPtr6 = xPtr5 + nE;
                xPtr7 = xPtr6 + nE;
                xPtr8 = xPtr7 + nE;
                xPtr9 = xPtr8 + nE;
                xPtr10 = xPtr9 + nE;
                xPtr11 = xPtr10 + nE;
                xPtr12 = xPtr11 + nE;
                xPtr13 = xPtr12 + nE;
                xPtr14 = xPtr13 + nE;
                xPtr15 = xPtr14 + nE;
                xPtr16 = xPtr15 + nE;
                xPtr17 = xPtr16 + nE;
                xPtr18 = xPtr17 + nE;
                xPtr19 = xPtr18 + nE;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    SFP0ptr = wmhSFP0 + offset;
                    SFP1ptr = wmhSFP1 + offset;
                    SFP2ptr = wmhSFP2 + offset;
                    SFP3ptr = wmhSFP3 + offset;
                    SFP4ptr = wmhSFP4 + offset;
                    SFP5ptr = wmhSFP5 + offset;
                    SFP6ptr = wmhSFP6 + offset;
                    SFP7ptr = wmhSFP7 + offset;
                    SFP8ptr = wmhSFP8 + offset;
                    SFP9ptr = wmhSFP9 + offset;
                    SFP10ptr = wmhSFP10 + offset;
                    SFP11ptr = wmhSFP11 + offset;
                    SFP12ptr = wmhSFP12 + offset;
                    SFP13ptr = wmhSFP13 + offset;
                    SFP14ptr = wmhSFP14 + offset;
                    SFP15ptr = wmhSFP15 + offset;
                    SFP16ptr = wmhSFP16 + offset;
                    SFP17ptr = wmhSFP17 + offset;
                    SFP18ptr = wmhSFP18 + offset;
                    SFP19ptr = wmhSFP19 + offset;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    x16 = (*SFP16ptr++) * YTmp;
                    x17 = (*SFP17ptr++) * YTmp;
                    x18 = (*SFP18ptr++) * YTmp;
                    x19 = (*SFP19ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                        x16 += (*SFP16ptr++) * YTmp;
                        x17 += (*SFP17ptr++) * YTmp;
                        x18 += (*SFP18ptr++) * YTmp;
                        x19 += (*SFP19ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    (*xPtr16++) += x16;
                    (*xPtr17++) += x17;
                    (*xPtr18++) += x18;
                    (*xPtr19++) += x19;
                    t_v++;
                    t_o++;
                }
                break;
        }
    }

    // isotropic compartments
    if (nISO > 0)
    {
        t_v = ISOv + ISOthreadsT[id];
        t_vEnd = ISOv + ISOthreadsT[id+1];
        xPtr0 = x + nIC*nF + nEC*nE + ISOthreadsT[id];
        switch (nISO)
        {
            case 1:
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    x0 = (*SFP0ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    t_v++;
                }
                break;
            
            case 2:
                xPtr1 = xPtr0 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    t_v++;
                }
                break;

            case 3:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    t_v++;
                }
                break;

            case 4:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    t_v++;
                }
                break;

            case 5:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    t_v++;
                }
                break;

            case 6:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    t_v++;
                }
                break;

            case 7:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    t_v++;
                }
                break;

            case 8:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    t_v++;
                }
                break;

            case 9:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    t_v++;
                }
                break;

            case 10:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    t_v++;
                }
                break;

            case 11:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    t_v++;
                }
                break;

            case 12:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    t_v++;
                }
                break;

            case 13:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    SFP12ptr = isoSFP12;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    t_v++;
                }
                break;

            case 14:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    SFP12ptr = isoSFP12;
                    SFP13ptr = isoSFP13;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    t_v++;
                }
                break;

            case 15:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    SFP12ptr = isoSFP12;
                    SFP13ptr = isoSFP13;
                    SFP14ptr = isoSFP14;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    t_v++;
                }
                break;

            case 16:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    SFP12ptr = isoSFP12;
                    SFP13ptr = isoSFP13;
                    SFP14ptr = isoSFP14;
                    SFP15ptr = isoSFP15;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    t_v++;
                }
                break;

            case 17:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                xPtr16 = xPtr15 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    SFP12ptr = isoSFP12;
                    SFP13ptr = isoSFP13;
                    SFP14ptr = isoSFP14;
                    SFP15ptr = isoSFP15;
                    SFP16ptr = isoSFP16;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    x16 = (*SFP16ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                        x16 += (*SFP16ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    (*xPtr16++) += x16;
                    t_v++;
                }
                break;

            case 18:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                xPtr16 = xPtr15 + nV;
                xPtr17 = xPtr16 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    SFP12ptr = isoSFP12;
                    SFP13ptr = isoSFP13;
                    SFP14ptr = isoSFP14;
                    SFP15ptr = isoSFP15;
                    SFP16ptr = isoSFP16;
                    SFP17ptr = isoSFP17;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    x16 = (*SFP16ptr++) * YTmp;
                    x17 = (*SFP17ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                        x16 += (*SFP16ptr++) * YTmp;
                        x17 += (*SFP17ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    (*xPtr16++) += x16;
                    (*xPtr17++) += x17;
                    t_v++;
                }
                break;

            case 19:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                xPtr16 = xPtr15 + nV;
                xPtr17 = xPtr16 + nV;
                xPtr18 = xPtr17 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    SFP12ptr = isoSFP12;
                    SFP13ptr = isoSFP13;
                    SFP14ptr = isoSFP14;
                    SFP15ptr = isoSFP15;
                    SFP16ptr = isoSFP16;
                    SFP17ptr = isoSFP17;
                    SFP18ptr = isoSFP18;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    x16 = (*SFP16ptr++) * YTmp;
                    x17 = (*SFP17ptr++) * YTmp;
                    x18 = (*SFP18ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                        x16 += (*SFP16ptr++) * YTmp;
                        x17 += (*SFP17ptr++) * YTmp;
                        x18 += (*SFP18ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    (*xPtr16++) += x16;
                    (*xPtr17++) += x17;
                    (*xPtr18++) += x18;
                    t_v++;
                }
                break;

            case 20:
                xPtr1 = xPtr0 + nV;
                xPtr2 = xPtr1 + nV;
                xPtr3 = xPtr2 + nV;
                xPtr4 = xPtr3 + nV;
                xPtr5 = xPtr4 + nV;
                xPtr6 = xPtr5 + nV;
                xPtr7 = xPtr6 + nV;
                xPtr8 = xPtr7 + nV;
                xPtr9 = xPtr8 + nV;
                xPtr10 = xPtr9 + nV;
                xPtr11 = xPtr10 + nV;
                xPtr12 = xPtr11 + nV;
                xPtr13 = xPtr12 + nV;
                xPtr14 = xPtr13 + nV;
                xPtr15 = xPtr14 + nV;
                xPtr16 = xPtr15 + nV;
                xPtr17 = xPtr16 + nV;
                xPtr18 = xPtr17 + nV;
                xPtr19 = xPtr18 + nV;
                while (t_v != t_vEnd)
                {
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    SFP0ptr = isoSFP0;
                    SFP1ptr = isoSFP1;
                    SFP2ptr = isoSFP2;
                    SFP3ptr = isoSFP3;
                    SFP4ptr = isoSFP4;
                    SFP5ptr = isoSFP5;
                    SFP6ptr = isoSFP6;
                    SFP7ptr = isoSFP7;
                    SFP8ptr = isoSFP8;
                    SFP9ptr = isoSFP9;
                    SFP10ptr = isoSFP10;
                    SFP11ptr = isoSFP11;
                    SFP12ptr = isoSFP12;
                    SFP13ptr = isoSFP13;
                    SFP14ptr = isoSFP14;
                    SFP15ptr = isoSFP15;
                    SFP16ptr = isoSFP16;
                    SFP17ptr = isoSFP17;
                    SFP18ptr = isoSFP18;
                    SFP19ptr = isoSFP19;
                    x0 = (*SFP0ptr++) * YTmp;
                    x1 = (*SFP1ptr++) * YTmp;
                    x2 = (*SFP2ptr++) * YTmp;
                    x3 = (*SFP3ptr++) * YTmp;
                    x4 = (*SFP4ptr++) * YTmp;
                    x5 = (*SFP5ptr++) * YTmp;
                    x6 = (*SFP6ptr++) * YTmp;
                    x7 = (*SFP7ptr++) * YTmp;
                    x8 = (*SFP8ptr++) * YTmp;
                    x9 = (*SFP9ptr++) * YTmp;
                    x10 = (*SFP10ptr++) * YTmp;
                    x11 = (*SFP11ptr++) * YTmp;
                    x12 = (*SFP12ptr++) * YTmp;
                    x13 = (*SFP13ptr++) * YTmp;
                    x14 = (*SFP14ptr++) * YTmp;
                    x15 = (*SFP15ptr++) * YTmp;
                    x16 = (*SFP16ptr++) * YTmp;
                    x17 = (*SFP17ptr++) * YTmp;
                    x18 = (*SFP18ptr++) * YTmp;
                    x19 = (*SFP19ptr++) * YTmp;
                    while (++YPtr != YPtrEnd)
                    {
                        YTmp = *YPtr;
                        x0 += (*SFP0ptr++) * YTmp;
                        x1 += (*SFP1ptr++) * YTmp;
                        x2 += (*SFP2ptr++) * YTmp;
                        x3 += (*SFP3ptr++) * YTmp;
                        x4 += (*SFP4ptr++) * YTmp;
                        x5 += (*SFP5ptr++) * YTmp;
                        x6 += (*SFP6ptr++) * YTmp;
                        x7 += (*SFP7ptr++) * YTmp;
                        x8 += (*SFP8ptr++) * YTmp;
                        x9 += (*SFP9ptr++) * YTmp;
                        x10 += (*SFP10ptr++) * YTmp;
                        x11 += (*SFP11ptr++) * YTmp;
                        x12 += (*SFP12ptr++) * YTmp;
                        x13 += (*SFP13ptr++) * YTmp;
                        x14 += (*SFP14ptr++) * YTmp;
                        x15 += (*SFP15ptr++) * YTmp;
                        x16 += (*SFP16ptr++) * YTmp;
                        x17 += (*SFP17ptr++) * YTmp;
                        x18 += (*SFP18ptr++) * YTmp;
                        x19 += (*SFP19ptr++) * YTmp;
                    }
                    (*xPtr0++) += x0;
                    (*xPtr1++) += x1;
                    (*xPtr2++) += x2;
                    (*xPtr3++) += x3;
                    (*xPtr4++) += x4;
                    (*xPtr5++) += x5;
                    (*xPtr6++) += x6;
                    (*xPtr7++) += x7;
                    (*xPtr8++) += x8;
                    (*xPtr9++) += x9;
                    (*xPtr10++) += x10;
                    (*xPtr11++) += x11;
                    (*xPtr12++) += x12;
                    (*xPtr13++) += x13;
                    (*xPtr14++) += x14;
                    (*xPtr15++) += x15;
                    (*xPtr16++) += x16;
                    (*xPtr17++) += x17;
                    (*xPtr18++) += x18;
                    (*xPtr19++) += x19;
                    t_v++;
                }
                break;
        }
    }

    pthread_exit( 0 );
}

///////////////////////////////
// Function called by CYTHON //
///////////////////////////////
void COMMIT_At(
    int _nF, int _n, int _nE, int _nV, int _nS, int _ndirs,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICeval, uint32_t *_ICv, uint16_t *_ICo, float *_ICl,
    uint32_t *_ECv, uint16_t *_ECo,
    uint32_t *_ISOv,
    float *_wmrSFP, float *_wmhSFP, float *_isoSFP,
    uint8_t* _ICthreadsT, uint32_t* _ECthreadsT, uint32_t* _ISOthreadsT,
    uint32_t _nIC, uint32_t _nEC, uint32_t _nISO, uint32_t _nThreads
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
    ICeval = _ICeval;
    ICv  = _ICv;
    ICo  = _ICo;
    ICl  = _ICl;
    ECv  = _ECv;
    ECo  = _ECo;
    ISOv = _ISOv;

    nIC = _nIC;
    nEC = _nEC;
    nISO = _nISO;

    switch (nIC)
    {
        case 1:
            wmrSFP0 = _wmrSFP;
            break;
        case 2:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            break;
        case 3:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            break;
        case 4:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            break;
        case 5:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            break;
        case 6:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            wmrSFP5 = wmrSFP4 + _ndirs*_nS;
            break;
        case 7:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            wmrSFP5 = wmrSFP4 + _ndirs*_nS;
            wmrSFP6 = wmrSFP5 + _ndirs*_nS;
            break;
        case 8:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            wmrSFP5 = wmrSFP4 + _ndirs*_nS;
            wmrSFP6 = wmrSFP5 + _ndirs*_nS;
            wmrSFP7 = wmrSFP6 + _ndirs*_nS;
            break;
        case 9:
            wmrSFP0 = _wmrSFP;
            wmrSFP1 = wmrSFP0 + _ndirs*_nS;
            wmrSFP2 = wmrSFP1 + _ndirs*_nS;
            wmrSFP3 = wmrSFP2 + _ndirs*_nS;
            wmrSFP4 = wmrSFP3 + _ndirs*_nS;
            wmrSFP5 = wmrSFP4 + _ndirs*_nS;
            wmrSFP6 = wmrSFP5 + _ndirs*_nS;
            wmrSFP7 = wmrSFP6 + _ndirs*_nS;
            wmrSFP8 = wmrSFP7 + _ndirs*_nS;
            break;
        case 10:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            break;
        case 11:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            break;
        case 12:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            break;
        case 13:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            break;
        case 14:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            break;
        case 15:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            break;
        case 16:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            break;
        case 17:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            wmrSFP16 = wmrSFP15 + _ndirs*_nS;
            break;
        case 18:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0 + _ndirs*_nS;
            wmrSFP2  = wmrSFP1 + _ndirs*_nS;
            wmrSFP3  = wmrSFP2 + _ndirs*_nS;
            wmrSFP4  = wmrSFP3 + _ndirs*_nS;
            wmrSFP5  = wmrSFP4 + _ndirs*_nS;
            wmrSFP6  = wmrSFP5 + _ndirs*_nS;
            wmrSFP7  = wmrSFP6 + _ndirs*_nS;
            wmrSFP8  = wmrSFP7 + _ndirs*_nS;
            wmrSFP9  = wmrSFP8 + _ndirs*_nS;
            wmrSFP10 = wmrSFP9 + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            wmrSFP16 = wmrSFP15 + _ndirs*_nS;
            wmrSFP17 = wmrSFP16 + _ndirs*_nS;
            break;
        case 19:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0  + _ndirs*_nS;
            wmrSFP2  = wmrSFP1  + _ndirs*_nS;
            wmrSFP3  = wmrSFP2  + _ndirs*_nS;
            wmrSFP4  = wmrSFP3  + _ndirs*_nS;
            wmrSFP5  = wmrSFP4  + _ndirs*_nS;
            wmrSFP6  = wmrSFP5  + _ndirs*_nS;
            wmrSFP7  = wmrSFP6  + _ndirs*_nS;
            wmrSFP8  = wmrSFP7  + _ndirs*_nS;
            wmrSFP9  = wmrSFP8  + _ndirs*_nS;
            wmrSFP10 = wmrSFP9  + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            wmrSFP16 = wmrSFP15 + _ndirs*_nS;
            wmrSFP17 = wmrSFP16 + _ndirs*_nS;
            wmrSFP18 = wmrSFP17 + _ndirs*_nS;
            break;
        case 20:
            wmrSFP0  = _wmrSFP;
            wmrSFP1  = wmrSFP0  + _ndirs*_nS;
            wmrSFP2  = wmrSFP1  + _ndirs*_nS;
            wmrSFP3  = wmrSFP2  + _ndirs*_nS;
            wmrSFP4  = wmrSFP3  + _ndirs*_nS;
            wmrSFP5  = wmrSFP4  + _ndirs*_nS;
            wmrSFP6  = wmrSFP5  + _ndirs*_nS;
            wmrSFP7  = wmrSFP6  + _ndirs*_nS;
            wmrSFP8  = wmrSFP7  + _ndirs*_nS;
            wmrSFP9  = wmrSFP8  + _ndirs*_nS;
            wmrSFP10 = wmrSFP9  + _ndirs*_nS;
            wmrSFP11 = wmrSFP10 + _ndirs*_nS;
            wmrSFP12 = wmrSFP11 + _ndirs*_nS;
            wmrSFP13 = wmrSFP12 + _ndirs*_nS;
            wmrSFP14 = wmrSFP13 + _ndirs*_nS;
            wmrSFP15 = wmrSFP14 + _ndirs*_nS;
            wmrSFP16 = wmrSFP15 + _ndirs*_nS;
            wmrSFP17 = wmrSFP16 + _ndirs*_nS;
            wmrSFP18 = wmrSFP17 + _ndirs*_nS;
            wmrSFP19 = wmrSFP18 + _ndirs*_nS;
            break;
    }

    switch (nEC)
    {
        case 1:
            wmhSFP0 = _wmhSFP;
            break;
        case 2:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            break;
        case 3:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            break;
        case 4:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            break;
        case 5:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            break;
        case 6:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            wmhSFP5 = wmhSFP4 + _ndirs*_nS;
            break;
        case 7:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            wmhSFP5 = wmhSFP4 + _ndirs*_nS;
            wmhSFP6 = wmhSFP5 + _ndirs*_nS;
            break;
        case 8:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            wmhSFP5 = wmhSFP4 + _ndirs*_nS;
            wmhSFP6 = wmhSFP5 + _ndirs*_nS;
            wmhSFP7 = wmhSFP6 + _ndirs*_nS;
            break;
        case 9:
            wmhSFP0 = _wmhSFP;
            wmhSFP1 = wmhSFP0 + _ndirs*_nS;
            wmhSFP2 = wmhSFP1 + _ndirs*_nS;
            wmhSFP3 = wmhSFP2 + _ndirs*_nS;
            wmhSFP4 = wmhSFP3 + _ndirs*_nS;
            wmhSFP5 = wmhSFP4 + _ndirs*_nS;
            wmhSFP6 = wmhSFP5 + _ndirs*_nS;
            wmhSFP7 = wmhSFP6 + _ndirs*_nS;
            wmhSFP8 = wmhSFP7 + _ndirs*_nS;
            break;
        case 10:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            break;
        case 11:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            break;
        case 12:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            break;
        case 13:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            break;
        case 14:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            break;
        case 15:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            break;
        case 16:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            break;
        case 17:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            wmhSFP16 = wmhSFP15 + _ndirs*_nS;
            break;
        case 18:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0 + _ndirs*_nS;
            wmhSFP2  = wmhSFP1 + _ndirs*_nS;
            wmhSFP3  = wmhSFP2 + _ndirs*_nS;
            wmhSFP4  = wmhSFP3 + _ndirs*_nS;
            wmhSFP5  = wmhSFP4 + _ndirs*_nS;
            wmhSFP6  = wmhSFP5 + _ndirs*_nS;
            wmhSFP7  = wmhSFP6 + _ndirs*_nS;
            wmhSFP8  = wmhSFP7 + _ndirs*_nS;
            wmhSFP9  = wmhSFP8 + _ndirs*_nS;
            wmhSFP10 = wmhSFP9 + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            wmhSFP16 = wmhSFP15 + _ndirs*_nS;
            wmhSFP17 = wmhSFP16 + _ndirs*_nS;
            break;
        case 19:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0  + _ndirs*_nS;
            wmhSFP2  = wmhSFP1  + _ndirs*_nS;
            wmhSFP3  = wmhSFP2  + _ndirs*_nS;
            wmhSFP4  = wmhSFP3  + _ndirs*_nS;
            wmhSFP5  = wmhSFP4  + _ndirs*_nS;
            wmhSFP6  = wmhSFP5  + _ndirs*_nS;
            wmhSFP7  = wmhSFP6  + _ndirs*_nS;
            wmhSFP8  = wmhSFP7  + _ndirs*_nS;
            wmhSFP9  = wmhSFP8  + _ndirs*_nS;
            wmhSFP10 = wmhSFP9  + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            wmhSFP16 = wmhSFP15 + _ndirs*_nS;
            wmhSFP17 = wmhSFP16 + _ndirs*_nS;
            wmhSFP18 = wmhSFP17 + _ndirs*_nS;
            break;
        case 20:
            wmhSFP0  = _wmhSFP;
            wmhSFP1  = wmhSFP0  + _ndirs*_nS;
            wmhSFP2  = wmhSFP1  + _ndirs*_nS;
            wmhSFP3  = wmhSFP2  + _ndirs*_nS;
            wmhSFP4  = wmhSFP3  + _ndirs*_nS;
            wmhSFP5  = wmhSFP4  + _ndirs*_nS;
            wmhSFP6  = wmhSFP5  + _ndirs*_nS;
            wmhSFP7  = wmhSFP6  + _ndirs*_nS;
            wmhSFP8  = wmhSFP7  + _ndirs*_nS;
            wmhSFP9  = wmhSFP8  + _ndirs*_nS;
            wmhSFP10 = wmhSFP9  + _ndirs*_nS;
            wmhSFP11 = wmhSFP10 + _ndirs*_nS;
            wmhSFP12 = wmhSFP11 + _ndirs*_nS;
            wmhSFP13 = wmhSFP12 + _ndirs*_nS;
            wmhSFP14 = wmhSFP13 + _ndirs*_nS;
            wmhSFP15 = wmhSFP14 + _ndirs*_nS;
            wmhSFP16 = wmhSFP15 + _ndirs*_nS;
            wmhSFP17 = wmhSFP16 + _ndirs*_nS;
            wmhSFP18 = wmhSFP17 + _ndirs*_nS;
            wmhSFP19 = wmhSFP18 + _ndirs*_nS;
            break;
    }

    switch (nISO)
    {
        case 1:
            isoSFP0 = _isoSFP;
            break;
        case 2:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            break;
        case 3:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            break;
        case 4:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            break;
        case 5:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            break;
        case 6:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            isoSFP5 = isoSFP4 + _nS;
            break;
        case 7:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            isoSFP5 = isoSFP4 + _nS;
            isoSFP6 = isoSFP5 + _nS;
            break;
        case 8:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            isoSFP5 = isoSFP4 + _nS;
            isoSFP6 = isoSFP5 + _nS;
            isoSFP7 = isoSFP6 + _nS;
            break;
        case 9:
            isoSFP0 = _isoSFP;
            isoSFP1 = isoSFP0 + _nS;
            isoSFP2 = isoSFP1 + _nS;
            isoSFP3 = isoSFP2 + _nS;
            isoSFP4 = isoSFP3 + _nS;
            isoSFP5 = isoSFP4 + _nS;
            isoSFP6 = isoSFP5 + _nS;
            isoSFP7 = isoSFP6 + _nS;
            isoSFP8 = isoSFP7 + _nS;
            break;
        case 10:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            break;
        case 11:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            break;
        case 12:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            break;
        case 13:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            break;
        case 14:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            break;
        case 15:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            break;
        case 16:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            break;
        case 17:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            isoSFP16 = isoSFP15 + _nS;
            break;
        case 18:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0 + _nS;
            isoSFP2  = isoSFP1 + _nS;
            isoSFP3  = isoSFP2 + _nS;
            isoSFP4  = isoSFP3 + _nS;
            isoSFP5  = isoSFP4 + _nS;
            isoSFP6  = isoSFP5 + _nS;
            isoSFP7  = isoSFP6 + _nS;
            isoSFP8  = isoSFP7 + _nS;
            isoSFP9  = isoSFP8 + _nS;
            isoSFP10 = isoSFP9 + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            isoSFP16 = isoSFP15 + _nS;
            isoSFP17 = isoSFP16 + _nS;
            break;
        case 19:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0  + _nS;
            isoSFP2  = isoSFP1  + _nS;
            isoSFP3  = isoSFP2  + _nS;
            isoSFP4  = isoSFP3  + _nS;
            isoSFP5  = isoSFP4  + _nS;
            isoSFP6  = isoSFP5  + _nS;
            isoSFP7  = isoSFP6  + _nS;
            isoSFP8  = isoSFP7  + _nS;
            isoSFP9  = isoSFP8  + _nS;
            isoSFP10 = isoSFP9  + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            isoSFP16 = isoSFP15 + _nS;
            isoSFP17 = isoSFP16 + _nS;
            isoSFP18 = isoSFP17 + _nS;
            break;
        case 20:
            isoSFP0  = _isoSFP;
            isoSFP1  = isoSFP0  + _nS;
            isoSFP2  = isoSFP1  + _nS;
            isoSFP3  = isoSFP2  + _nS;
            isoSFP4  = isoSFP3  + _nS;
            isoSFP5  = isoSFP4  + _nS;
            isoSFP6  = isoSFP5  + _nS;
            isoSFP7  = isoSFP6  + _nS;
            isoSFP8  = isoSFP7  + _nS;
            isoSFP9  = isoSFP8  + _nS;
            isoSFP10 = isoSFP9  + _nS;
            isoSFP11 = isoSFP10 + _nS;
            isoSFP12 = isoSFP11 + _nS;
            isoSFP13 = isoSFP12 + _nS;
            isoSFP14 = isoSFP13 + _nS;
            isoSFP15 = isoSFP14 + _nS;
            isoSFP16 = isoSFP15 + _nS;
            isoSFP17 = isoSFP16 + _nS;
            isoSFP18 = isoSFP17 + _nS;
            isoSFP19 = isoSFP18 + _nS;
            break;
    }

    ICthreadsT  = _ICthreadsT;
    ECthreadsT  = _ECthreadsT;
    ISOthreadsT = _ISOthreadsT;

    // Run SEPARATE THREADS to perform the multiplication
    pthread_t threads[MAX_THREADS];
    int t;
    for(t=0; t<_nThreads ; t++)
        pthread_create( &threads[t], NULL, COMMIT_At__block, (void *) (long int)t );
    for(t=0; t<_nThreads ; t++)
        pthread_join( threads[t], NULL );
    return;
}


// more global variables
int nSf;
double *icSFB0, *icSFB1, *icSFB2, *icSFB3, *icSFB4, *icSFB5, *icSFB6, *icSFB7, *icSFB8, *icSFB9;
uint32_t *ICp;
uint32_t nICs;

//////////////////////////////////////////////////////////
// Compute a sub-block of the A*x MATRIX-VECTOR product //
//////////////////////////////////////////////////////////
void* COMMIT_A__block_nolut( void *ptr )
{
    int      id = (long)ptr;
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, w;
    uint32_t *eval0, *eval1, *eval2, *eval3, *eval4, *eval5, *eval6, *eval7, *eval8, *eval9;
    double   x0_tmp, x1_tmp, x2_tmp, x3_tmp, x4_tmp, x5_tmp, x6_tmp, x7_tmp, x8_tmp, x9_tmp;
    double   *x_Ptr0, *x_Ptr1, *x_Ptr2, *x_Ptr3, *x_Ptr4, *x_Ptr5, *x_Ptr6, *x_Ptr7, *x_Ptr8, *x_Ptr9;
    double   *Yptr, *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr;
    uint32_t *t_v, *t_vEnd, *t_f, *t_p;
    float    *t_l;
    int      offset;
    double   *xPtr;

    // intra-cellular compartments
    if (nICs > 0)
    {
        // DCT basis functions 
        t_v = ICv + ICthreads[id];
        t_vEnd = ICv + ICthreads[id+1];
        t_l = ICl + ICthreads[id];
        t_f = ICf + ICthreads[id];
        t_p = ICp + ICthreads[id];

        switch (nICs)
        {
            case 1:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    // eval0 points to the same position as x_Ptr0
                    eval0 = ICeval + *t_f;
                    x0 = *x_Ptr0 * (double)(*eval0); 
                    x0_tmp = 0;
                    if (x0 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l); //* (double)(ICeval[*t_f]);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        x0 *= (*SFP0ptr);
                        (*Yptr++) += w * (x0);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    t_p++;

                }
                break;
            case 2:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    if (x0 != 0 || x1 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        (*Yptr++) += w * (x0 + x1_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
            case 3:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x_Ptr2 = x_Ptr1 + nF;
                    eval2 = eval1 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x2 = *x_Ptr2 * (double)(*eval2);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    x2_tmp = 0;
                    if (x0 != 0 || x1 != 0 || x2 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        SFP2ptr = icSFB2 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        x2_tmp += x2 * (*SFP2ptr);
                        (*Yptr++) += w * (x0 + x1_tmp + x2_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
            case 4:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x_Ptr2 = x_Ptr1 + nF;
                    eval2 = eval1 + nF;
                    x_Ptr3 = x_Ptr2 + nF;
                    eval3 = eval2 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x2 = *x_Ptr2 * (double)(*eval2);
                    x3 = *x_Ptr3 * (double)(*eval3);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    x2_tmp = 0;
                    x3_tmp = 0;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        SFP2ptr = icSFB2 + offset;
                        SFP3ptr = icSFB3 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        x2_tmp += x2 * (*SFP2ptr);
                        x3_tmp += x3 * (*SFP3ptr);
                        (*Yptr++) += w * (x0 + x1_tmp + x2_tmp + x3_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
            case 5:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x_Ptr2 = x_Ptr1 + nF;
                    eval2 = eval1 + nF;
                    x_Ptr3 = x_Ptr2 + nF;
                    eval3 = eval2 + nF;
                    x_Ptr4 = x_Ptr3 + nF;
                    eval4 = eval3 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x2 = *x_Ptr2 * (double)(*eval2);
                    x3 = *x_Ptr3 * (double)(*eval3);
                    x4 = *x_Ptr4 * (double)(*eval4);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    x2_tmp = 0;
                    x3_tmp = 0;
                    x4_tmp = 0;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        SFP2ptr = icSFB2 + offset;
                        SFP3ptr = icSFB3 + offset;
                        SFP4ptr = icSFB4 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        x2_tmp += x2 * (*SFP2ptr);
                        x3_tmp += x3 * (*SFP3ptr);
                        x4_tmp += x4 * (*SFP4ptr);
                        (*Yptr++) += w * (x0 + x1_tmp + x2_tmp + x3_tmp + x4_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
            case 6:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x_Ptr2 = x_Ptr1 + nF;
                    eval2 = eval1 + nF;
                    x_Ptr3 = x_Ptr2 + nF;
                    eval3 = eval2 + nF;
                    x_Ptr4 = x_Ptr3 + nF;
                    eval4 = eval3 + nF;
                    x_Ptr5 = x_Ptr4 + nF;
                    eval5 = eval4 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x2 = *x_Ptr2 * (double)(*eval2);
                    x3 = *x_Ptr3 * (double)(*eval3);
                    x4 = *x_Ptr4 * (double)(*eval4);
                    x5 = *x_Ptr5 * (double)(*eval5);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    x2_tmp = 0;
                    x3_tmp = 0;
                    x4_tmp = 0;
                    x5_tmp = 0;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        SFP2ptr = icSFB2 + offset;
                        SFP3ptr = icSFB3 + offset;
                        SFP4ptr = icSFB4 + offset;
                        SFP5ptr = icSFB5 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        x2_tmp += x2 * (*SFP2ptr);
                        x3_tmp += x3 * (*SFP3ptr);
                        x4_tmp += x4 * (*SFP4ptr);
                        x5_tmp += x5 * (*SFP5ptr);
                        (*Yptr++) += w * (x0 + x1_tmp + x2_tmp + x3_tmp + x4_tmp + x5_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
            case 7:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x_Ptr2 = x_Ptr1 + nF;
                    eval2 = eval1 + nF;
                    x_Ptr3 = x_Ptr2 + nF;
                    eval3 = eval2 + nF;
                    x_Ptr4 = x_Ptr3 + nF;
                    eval4 = eval3 + nF;
                    x_Ptr5 = x_Ptr4 + nF;
                    eval5 = eval4 + nF;
                    x_Ptr6 = x_Ptr5 + nF;
                    eval6 = eval5 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x2 = *x_Ptr2 * (double)(*eval2);
                    x3 = *x_Ptr3 * (double)(*eval3);
                    x4 = *x_Ptr4 * (double)(*eval4);
                    x5 = *x_Ptr5 * (double)(*eval5);
                    x6 = *x_Ptr6 * (double)(*eval6);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    x2_tmp = 0;
                    x3_tmp = 0;
                    x4_tmp = 0;
                    x5_tmp = 0;
                    x6_tmp = 0;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        SFP2ptr = icSFB2 + offset;
                        SFP3ptr = icSFB3 + offset;
                        SFP4ptr = icSFB4 + offset;
                        SFP5ptr = icSFB5 + offset;
                        SFP6ptr = icSFB6 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        x2_tmp += x2 * (*SFP2ptr);
                        x3_tmp += x3 * (*SFP3ptr);
                        x4_tmp += x4 * (*SFP4ptr);
                        x5_tmp += x5 * (*SFP5ptr);
                        x6_tmp += x6 * (*SFP6ptr);
                        (*Yptr++) += w * (x0 + x1_tmp + x2_tmp + x3_tmp + x4_tmp + x5_tmp + x6_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
            case 8:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x_Ptr2 = x_Ptr1 + nF;
                    eval2 = eval1 + nF;
                    x_Ptr3 = x_Ptr2 + nF;
                    eval3 = eval2 + nF;
                    x_Ptr4 = x_Ptr3 + nF;
                    eval4 = eval3 + nF;
                    x_Ptr5 = x_Ptr4 + nF;
                    eval5 = eval4 + nF;
                    x_Ptr6 = x_Ptr5 + nF;
                    eval6 = eval5 + nF;
                    x_Ptr7 = x_Ptr6 + nF;
                    eval7 = eval6 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x2 = *x_Ptr2 * (double)(*eval2);
                    x3 = *x_Ptr3 * (double)(*eval3);
                    x4 = *x_Ptr4 * (double)(*eval4);
                    x5 = *x_Ptr5 * (double)(*eval5);
                    x6 = *x_Ptr6 * (double)(*eval6);
                    x7 = *x_Ptr7 * (double)(*eval7);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    x2_tmp = 0;
                    x3_tmp = 0;
                    x4_tmp = 0;
                    x5_tmp = 0;
                    x6_tmp = 0;
                    x7_tmp = 0;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        SFP2ptr = icSFB2 + offset;
                        SFP3ptr = icSFB3 + offset;
                        SFP4ptr = icSFB4 + offset;
                        SFP5ptr = icSFB5 + offset;
                        SFP6ptr = icSFB6 + offset;
                        SFP7ptr = icSFB7 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        x2_tmp += x2 * (*SFP2ptr);
                        x3_tmp += x3 * (*SFP3ptr);
                        x4_tmp += x4 * (*SFP4ptr);
                        x5_tmp += x5 * (*SFP5ptr);
                        x6_tmp += x6 * (*SFP6ptr);
                        x7_tmp += x7 * (*SFP7ptr);
                        (*Yptr++) += w * (x0 + x1_tmp + x2_tmp + x3_tmp + x4_tmp + x5_tmp + x6_tmp + x7_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
            case 9:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x_Ptr2 = x_Ptr1 + nF;
                    eval2 = eval1 + nF;
                    x_Ptr3 = x_Ptr2 + nF;
                    eval3 = eval2 + nF;
                    x_Ptr4 = x_Ptr3 + nF;
                    eval4 = eval3 + nF;
                    x_Ptr5 = x_Ptr4 + nF;
                    eval5 = eval4 + nF;
                    x_Ptr6 = x_Ptr5 + nF;
                    eval6 = eval5 + nF;
                    x_Ptr7 = x_Ptr6 + nF;
                    eval7 = eval6 + nF;
                    x_Ptr8 = x_Ptr7 + nF;
                    eval8 = eval7 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x2 = *x_Ptr2 * (double)(*eval2);
                    x3 = *x_Ptr3 * (double)(*eval3);
                    x4 = *x_Ptr4 * (double)(*eval4);
                    x5 = *x_Ptr5 * (double)(*eval5);
                    x6 = *x_Ptr6 * (double)(*eval6);
                    x7 = *x_Ptr7 * (double)(*eval7);
                    x8 = *x_Ptr8 * (double)(*eval8);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    x2_tmp = 0;
                    x3_tmp = 0;
                    x4_tmp = 0;
                    x5_tmp = 0;
                    x6_tmp = 0;
                    x7_tmp = 0;
                    x8_tmp = 0;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        SFP2ptr = icSFB2 + offset;
                        SFP3ptr = icSFB3 + offset;
                        SFP4ptr = icSFB4 + offset;
                        SFP5ptr = icSFB5 + offset;
                        SFP6ptr = icSFB6 + offset;
                        SFP7ptr = icSFB7 + offset;
                        SFP8ptr = icSFB8 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        x2_tmp += x2 * (*SFP2ptr);
                        x3_tmp += x3 * (*SFP3ptr);
                        x4_tmp += x4 * (*SFP4ptr);
                        x5_tmp += x5 * (*SFP5ptr);
                        x6_tmp += x6 * (*SFP6ptr);
                        x7_tmp += x7 * (*SFP7ptr);
                        x8_tmp += x8 * (*SFP8ptr);
                        (*Yptr++) += w * (x0 + x1_tmp + x2_tmp + x3_tmp + x4_tmp + x5_tmp + x6_tmp + x7_tmp + x8_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
            case 10:
                while (t_v != t_vEnd)
                {
                    x_Ptr0 = x + *t_f;
                    eval0 = ICeval + *t_f;
                    x_Ptr1 = x_Ptr0 + nF;
                    eval1 = eval0 + nF;
                    x_Ptr2 = x_Ptr1 + nF;
                    eval2 = eval1 + nF;
                    x_Ptr3 = x_Ptr2 + nF;
                    eval3 = eval2 + nF;
                    x_Ptr4 = x_Ptr3 + nF;
                    eval4 = eval3 + nF;
                    x_Ptr5 = x_Ptr4 + nF;
                    eval5 = eval4 + nF;
                    x_Ptr6 = x_Ptr5 + nF;
                    eval6 = eval5 + nF;
                    x_Ptr7 = x_Ptr6 + nF;
                    eval7 = eval6 + nF;
                    x_Ptr8 = x_Ptr7 + nF;
                    eval8 = eval7 + nF;
                    x_Ptr9 = x_Ptr8 + nF;
                    eval9 = eval8 + nF;
                    x0 = *x_Ptr0 * (double)(*eval0);
                    x1 = *x_Ptr1 * (double)(*eval1);
                    x2 = *x_Ptr2 * (double)(*eval2);
                    x3 = *x_Ptr3 * (double)(*eval3);
                    x4 = *x_Ptr4 * (double)(*eval4);
                    x5 = *x_Ptr5 * (double)(*eval5);
                    x6 = *x_Ptr6 * (double)(*eval6);
                    x7 = *x_Ptr7 * (double)(*eval7);
                    x8 = *x_Ptr8 * (double)(*eval8);
                    x9 = *x_Ptr9 * (double)(*eval9);
                    x0_tmp = 0;
                    x1_tmp = 0;
                    x2_tmp = 0;
                    x3_tmp = 0;
                    x4_tmp = 0;
                    x5_tmp = 0;
                    x6_tmp = 0;
                    x7_tmp = 0;
                    x8_tmp = 0;
                    x9_tmp = 0;
                    if (x0 != 0 || x1 != 0 || x2 != 0 || x3 != 0 || x4 != 0 || x5 != 0 || x6 != 0 || x7 != 0 || x8 != 0 || x9 != 0)
                    {
                        Yptr = Y + (*t_v);
                        w = (double)(*t_l);
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        SFP1ptr = icSFB1 + offset;
                        SFP2ptr = icSFB2 + offset;
                        SFP3ptr = icSFB3 + offset;
                        SFP4ptr = icSFB4 + offset;
                        SFP5ptr = icSFB5 + offset;
                        SFP6ptr = icSFB6 + offset;
                        SFP7ptr = icSFB7 + offset;
                        SFP8ptr = icSFB8 + offset;
                        SFP9ptr = icSFB9 + offset;
                        x0 *= (*SFP0ptr);
                        x1_tmp += x1 * (*SFP1ptr);
                        x2_tmp += x2 * (*SFP2ptr);
                        x3_tmp += x3 * (*SFP3ptr);
                        x4_tmp += x4 * (*SFP4ptr);
                        x5_tmp += x5 * (*SFP5ptr);
                        x6_tmp += x6 * (*SFP6ptr);
                        x7_tmp += x7 * (*SFP7ptr);
                        x8_tmp += x8 * (*SFP8ptr);
                        x9_tmp += x9 * (*SFP9ptr);
                        (*Yptr++) += w * (x0 + x1_tmp + x2_tmp + x3_tmp + x4_tmp + x5_tmp + x6_tmp + x7_tmp + x8_tmp + x9_tmp);
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                }
                break;
        }
    }

    // isotropic compartments
    if (nISO > 0)
    {
        t_v = ISOv + ISOthreads[id];
        t_vEnd = ISOv + ISOthreads[id+1];
        xPtr = x + nF + ISOthreads[id];

        while (t_v != t_vEnd)
        {
            x0 = *xPtr++;
            if (x0 != 0)
                Y[*t_v] += x0;
            t_v++;
        }
    }

    pthread_exit( 0 );
}

///////////////////////////////
// Function called by CYTHON //
///////////////////////////////
void COMMIT_A_nolut(
    int _nF, int _n, int _nSf,
    double* _vIN, double* _vOUT,
    uint32_t *_ICf, uint32_t *_ICeval, uint32_t *_ICv, float *_ICl, uint32_t *_ICp,
    uint32_t *_ISOv,
    double *_ICmod,
    uint32_t* _ICthreads, uint32_t* _ISOthreads,
    uint32_t _nICs, uint32_t _nISO, uint32_t _nThreads
)
{
    nF = _nF;
    n  = _n;
    nSf = _nSf;

    x = _vIN;
    Y = _vOUT;

    ICf  = _ICf;
    ICeval = _ICeval;
    ICv  = _ICv;
    ICl  = _ICl;
    ICp  = _ICp;
    ISOv = _ISOv;

    nICs = _nICs;
    nISO = _nISO;

    switch (nICs)
    {
        case 1:
            icSFB0 = _ICmod;
            break;
        case 2:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            break;
        case 3:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            break;
        case 4:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            break;
        case 5:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            break;
        case 6:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            icSFB5 = icSFB4 + nSf;
            break;
        case 7:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            icSFB5 = icSFB4 + nSf;
            icSFB6 = icSFB5 + nSf;
            break;
        case 8:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            icSFB5 = icSFB4 + nSf;
            icSFB6 = icSFB5 + nSf;
            icSFB7 = icSFB6 + nSf;
            break;
        case 9:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            icSFB5 = icSFB4 + nSf;
            icSFB6 = icSFB5 + nSf;
            icSFB7 = icSFB6 + nSf;
            icSFB8 = icSFB7 + nSf;
            break;
        case 10:
            icSFB0  = _ICmod;
            icSFB1  = icSFB0 + nSf;
            icSFB2  = icSFB1 + nSf;
            icSFB3  = icSFB2 + nSf;
            icSFB4  = icSFB3 + nSf;
            icSFB5  = icSFB4 + nSf;
            icSFB6  = icSFB5 + nSf;
            icSFB7  = icSFB6 + nSf;
            icSFB8  = icSFB7 + nSf;
            icSFB9  = icSFB8 + nSf;
            break;
    }

    ICthreads  = _ICthreads;
    ISOthreads = _ISOthreads;

    // Run SEPARATE THREADS to perform the multiplication
    pthread_t threads[MAX_THREADS];
    int t;
    for(t=0; t<_nThreads ; t++)
        pthread_create( &threads[t], NULL, COMMIT_A__block_nolut, (void *) (long int)t );
    for(t=0; t<_nThreads ; t++)
        pthread_join( threads[t], NULL );
    return;
}


///////////////////////////////////////////////////////////
// Compute a sub-block of the A'*y MAtRIX-VECTOR product //
///////////////////////////////////////////////////////////
void* COMMIT_At__block_nolut( void *ptr )
{
    int      id = (long)ptr;
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, w, Y_tmp;
    uint32_t *eval0, *eval1, *eval2, *eval3, *eval4, *eval5, *eval6, *eval7, *eval8, *eval9;
    double   *Yptr, *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr;
    uint32_t *t_v, *t_vEnd, *t_f, *t_p;
    float    *t_l;
    uint8_t  *t_t;
    int      offset;
    double   *xPtr;

    // intra-cellular compartments
    if (nICs > 0)
    {
        t_v = ICv;
        t_vEnd = ICv + n;
        t_l = ICl;
        t_f = ICf;
        t_p = ICp;
        t_t = ICthreadsT;

        switch (nICs)
        {
            case 1:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                    t_t++;
                }
                break;
            case 2:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    t_p++;
                    t_t++;
                }
                break;
            case 3:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        SFP2ptr = icSFB2 + offset;
                        eval2 = eval1 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        x2 = (*SFP2ptr) * Y_tmp * (double)(*eval2);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                    t_t++;
                }
                break;
            case 4:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        SFP2ptr = icSFB2 + offset;
                        eval2 = eval1 + nF;
                        SFP3ptr = icSFB3 + offset;
                        eval3 = eval2 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        x2 = (*SFP2ptr) * Y_tmp * (double)(*eval2);
                        x3 = (*SFP3ptr) * Y_tmp * (double)(*eval3);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                    t_t++;
                }
                break;
            case 5:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        SFP2ptr = icSFB2 + offset;
                        eval2 = eval1 + nF;
                        SFP3ptr = icSFB3 + offset;
                        eval3 = eval2 + nF;
                        SFP4ptr = icSFB4 + offset;
                        eval4 = eval3 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        x2 = (*SFP2ptr) * Y_tmp * (double)(*eval2);
                        x3 = (*SFP3ptr) * Y_tmp * (double)(*eval3);
                        x4 = (*SFP4ptr) * Y_tmp * (double)(*eval4);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    t_p++;
                    t_t++;
                }
                break;
            case 6:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        SFP2ptr = icSFB2 + offset;
                        eval2 = eval1 + nF;
                        SFP3ptr = icSFB3 + offset;
                        eval3 = eval2 + nF;
                        SFP4ptr = icSFB4 + offset;
                        eval4 = eval3 + nF;
                        SFP5ptr = icSFB5 + offset;
                        eval5 = eval4 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        x2 = (*SFP2ptr) * Y_tmp * (double)(*eval2);
                        x3 = (*SFP3ptr) * Y_tmp * (double)(*eval3);
                        x4 = (*SFP4ptr) * Y_tmp * (double)(*eval4);
                        x5 = (*SFP5ptr) * Y_tmp * (double)(*eval5);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                    t_t++;
                }
                break;
            case 7:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        SFP2ptr = icSFB2 + offset;
                        eval2 = eval1 + nF;
                        SFP3ptr = icSFB3 + offset;
                        eval3 = eval2 + nF;
                        SFP4ptr = icSFB4 + offset;
                        eval4 = eval3 + nF;
                        SFP5ptr = icSFB5 + offset;
                        eval5 = eval4 + nF;
                        SFP6ptr = icSFB6 + offset;
                        eval6 = eval5 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        x2 = (*SFP2ptr) * Y_tmp * (double)(*eval2);
                        x3 = (*SFP3ptr) * Y_tmp * (double)(*eval3);
                        x4 = (*SFP4ptr) * Y_tmp * (double)(*eval4);
                        x5 = (*SFP5ptr) * Y_tmp * (double)(*eval5);
                        x6 = (*SFP6ptr) * Y_tmp * (double)(*eval6);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                    t_t++;
                }
                break;
            case 8:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        SFP2ptr = icSFB2 + offset;
                        eval2 = eval1 + nF;
                        SFP3ptr = icSFB3 + offset;
                        eval3 = eval2 + nF;
                        SFP4ptr = icSFB4 + offset;
                        eval4 = eval3 + nF;
                        SFP5ptr = icSFB5 + offset;
                        eval5 = eval4 + nF;
                        SFP6ptr = icSFB6 + offset;
                        eval6 = eval5 + nF;
                        SFP7ptr = icSFB7 + offset;
                        eval7 = eval6 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        x2 = (*SFP2ptr) * Y_tmp * (double)(*eval2);
                        x3 = (*SFP3ptr) * Y_tmp * (double)(*eval3);
                        x4 = (*SFP4ptr) * Y_tmp * (double)(*eval4);
                        x5 = (*SFP5ptr) * Y_tmp * (double)(*eval5);
                        x6 = (*SFP6ptr) * Y_tmp * (double)(*eval6);
                        x7 = (*SFP7ptr) * Y_tmp * (double)(*eval7);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                    t_t++;
                }
                break;
            case 9:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        SFP2ptr = icSFB2 + offset;
                        eval2 = eval1 + nF;
                        SFP3ptr = icSFB3 + offset;
                        eval3 = eval2 + nF;
                        SFP4ptr = icSFB4 + offset;
                        eval4 = eval3 + nF;
                        SFP5ptr = icSFB5 + offset;
                        eval5 = eval4 + nF;
                        SFP6ptr = icSFB6 + offset;
                        eval6 = eval5 + nF;
                        SFP7ptr = icSFB7 + offset;
                        eval7 = eval6 + nF;
                        SFP8ptr = icSFB8 + offset;
                        eval8 = eval7 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        x2 = (*SFP2ptr) * Y_tmp * (double)(*eval2);
                        x3 = (*SFP3ptr) * Y_tmp * (double)(*eval3);
                        x4 = (*SFP4ptr) * Y_tmp * (double)(*eval4);
                        x5 = (*SFP5ptr) * Y_tmp * (double)(*eval5);
                        x6 = (*SFP6ptr) * Y_tmp * (double)(*eval6);
                        x7 = (*SFP7ptr) * Y_tmp * (double)(*eval7);
                        x8 = (*SFP8ptr) * Y_tmp * (double)(*eval8);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                    t_t++;
                }
                break;
            case 10:
                while (t_v != t_vEnd)
                {
                    if (*t_t == id)
                    {
                        Yptr = Y + (*t_v);
                        Y_tmp = *Yptr;
                        offset = (*t_p) - 1;
                        SFP0ptr = icSFB0 + offset;
                        eval0 = ICeval + *t_f;
                        SFP1ptr = icSFB1 + offset;
                        eval1 = eval0 + nF;
                        SFP2ptr = icSFB2 + offset;
                        eval2 = eval1 + nF;
                        SFP3ptr = icSFB3 + offset;
                        eval3 = eval2 + nF;
                        SFP4ptr = icSFB4 + offset;
                        eval4 = eval3 + nF;
                        SFP5ptr = icSFB5 + offset;
                        eval5 = eval4 + nF;
                        SFP6ptr = icSFB6 + offset;
                        eval6 = eval5 + nF;
                        SFP7ptr = icSFB7 + offset;
                        eval7 = eval6 + nF;
                        SFP8ptr = icSFB8 + offset;
                        eval8 = eval7 + nF;
                        SFP9ptr = icSFB9 + offset;
                        eval9 = eval8 + nF;
                        x0 = (*SFP0ptr) * Y_tmp * (double)(*eval0);
                        x1 = (*SFP1ptr) * Y_tmp * (double)(*eval1);
                        x2 = (*SFP2ptr) * Y_tmp * (double)(*eval2);
                        x3 = (*SFP3ptr) * Y_tmp * (double)(*eval3);
                        x4 = (*SFP4ptr) * Y_tmp * (double)(*eval4);
                        x5 = (*SFP5ptr) * Y_tmp * (double)(*eval5);
                        x6 = (*SFP6ptr) * Y_tmp * (double)(*eval6);
                        x7 = (*SFP7ptr) * Y_tmp * (double)(*eval7);
                        x8 = (*SFP8ptr) * Y_tmp * (double)(*eval8);
                        x9 = (*SFP9ptr) * Y_tmp * (double)(*eval9);
                        w = (double)(*t_l);
                        x[*t_f] += w * x0;
                        x[*t_f+nF] += w * x1;
                        x[*t_f+2*nF] += w * x2;
                        x[*t_f+3*nF] += w * x3;
                        x[*t_f+4*nF] += w * x4;
                        x[*t_f+5*nF] += w * x5;
                        x[*t_f+6*nF] += w * x6;
                        x[*t_f+7*nF] += w * x7;
                        x[*t_f+8*nF] += w * x8;
                        x[*t_f+9*nF] += w * x9;
                    }
                    t_f++;
                    t_v++;
                    t_l++;
                    
                    t_p++;
                    t_t++;
                }
                break;
        }
    }

    // isotropic compartments
    if (nISO > 0)
    {
        t_v = ISOv + ISOthreadsT[id];
        t_vEnd = ISOv + ISOthreadsT[id+1];
        xPtr = x + nF + ISOthreadsT[id];

        while (t_v != t_vEnd)
            (*xPtr++) += Y[*t_v++];
    }
}

///////////////////////////////
// Function called by CYTHON //
///////////////////////////////
void COMMIT_At_nolut(
    int _nF, int _n, int _nSf,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICeval, uint32_t *_ICv, float *_ICl, uint32_t *_ICp,
    uint32_t *_ISOv,
    double *_ICmod,
    uint8_t* _ICthreadsT, uint32_t* _ISOthreadsT,
    uint32_t _nICs, uint32_t _nISO, uint32_t _nThreads
)
{
    nF = _nF;
    n  = _n;
    nSf = _nSf;

    x = _vOUT;
    Y = _vIN;
    
    ICf  = _ICf;
    ICeval = _ICeval;
    ICv  = _ICv;
    ICl  = _ICl;
    ICp  = _ICp;
    ISOv = _ISOv;

    nICs = _nICs;
    nISO = _nISO;

    switch (nICs)
    {
        case 1:
            icSFB0 = _ICmod;
            break;
        case 2:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            break;
        case 3:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            break;
        case 4:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            break;
        case 5:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            break;
        case 6:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            icSFB5 = icSFB4 + nSf;
            break;
        case 7:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            icSFB5 = icSFB4 + nSf;
            icSFB6 = icSFB5 + nSf;
            break;
        case 8:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            icSFB5 = icSFB4 + nSf;
            icSFB6 = icSFB5 + nSf;
            icSFB7 = icSFB6 + nSf;
            break;
        case 9:
            icSFB0 = _ICmod;
            icSFB1 = icSFB0 + nSf;
            icSFB2 = icSFB1 + nSf;
            icSFB3 = icSFB2 + nSf;
            icSFB4 = icSFB3 + nSf;
            icSFB5 = icSFB4 + nSf;
            icSFB6 = icSFB5 + nSf;
            icSFB7 = icSFB6 + nSf;
            icSFB8 = icSFB7 + nSf;
            break;
        case 10:
            icSFB0  = _ICmod;
            icSFB1  = icSFB0 + nSf;
            icSFB2  = icSFB1 + nSf;
            icSFB3  = icSFB2 + nSf;
            icSFB4  = icSFB3 + nSf;
            icSFB5  = icSFB4 + nSf;
            icSFB6  = icSFB5 + nSf;
            icSFB7  = icSFB6 + nSf;
            icSFB8  = icSFB7 + nSf;
            icSFB9  = icSFB8 + nSf;
            break;
    }

    ICthreadsT  = _ICthreadsT;
    ISOthreadsT = _ISOthreadsT;

    // Run SEPARATE THREADS to perform the multiplication
    pthread_t threads[MAX_THREADS];
    int t;
    for(t=0; t<_nThreads ; t++)
        pthread_create( &threads[t], NULL, COMMIT_At__block_nolut, (void *) (long int)t );
    for(t=0; t<_nThreads ; t++)
        pthread_join( threads[t], NULL );
    return;
}
