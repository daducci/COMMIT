from typing import NoReturn
from os.path import join as path_join


def add_imports() -> str:
    s: str = '''\
#include <pthread.h>
#include <stdint.h>

// max number of threads
#define MAX_THREADS 255\n\n'''
    return s


def add_global_variables() -> str:
    s: str = '''\
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
uint32_t    nIC, nEC, nISO;\n\n\n'''
    return s


def add_commit_a_block() -> str:
    def add_intracellular_compartments() -> str:
        s: str = '''\
    // intra-cellular compartments
    if (nIC > 0)
    {
        t_v = ICv + ICthreads[id];
        t_vEnd = ICv + ICthreads[id+1];
        t_o = ICo + ICthreads[id];
        t_l = ICl + ICthreads[id];
        t_f = ICf + ICthreads[id];
        switch (nIC)
        {'''
        s1: str = '''\
xPtr0 = x + (*t_f);
                    eval0 = ICeval + *t_f;
                    x0 = *xPtr0 * (double)(*eval0);'''
        s2: str = 'if (x0 != 0'
        s3: str = 'SFP0ptr = wmrSFP0 + offset;'
        s4: str = 'x0 * (*SFP0ptr++)'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

                    xPtr{i} = xPtr{i - 1} + nF;
                    eval{i} = eval{i - 1} + nF;
                    x{i} = *xPtr{i} * (double)(*eval{i});'''
                s2 += f' || x{i} != 0'
                s3 += f'''\

                        SFP{i}ptr = wmrSFP{i} + offset;'''
                s4 += f' + x{i} * (*SFP{i}ptr++)'
            s += f'''\

            case {i + 1}:
                while (t_v != t_vEnd)
                {{
                    {s1}
                    {s2})
                    {{
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        {s3}
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += w * ({s4});
                    }}
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                }}
                break;'''

        s += '''\

        }
    }\n\n'''
        return s

    def add_extracellular_compartments() -> str:
        s: str = '''\
    // extra-cellular compartments
    if (nEC > 0)
    {
        t_v = ECv + ECthreads[id];
        t_vEnd = ECv + ECthreads[id+1];
        t_o = ECo + ECthreads[id];
        xPtr0 = x + nIC*nF + ECthreads[id];
        switch (nEC)
        {'''
        s1: str = ''
        s2: str = 'x0 = *xPtr0++;'
        s3: str = 'if (x0 != 0'
        s4: str = 'SFP0ptr = wmhSFP0 + offset;'
        s5: str = 'x0 * (*SFP0ptr++)'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

                xPtr{i} = xPtr{i - 1} + nE;'''
                s2 += f'''\

                    x{i} = *xPtr{i}++;'''
                s3 += f' || x{i} != 0'
                s4 += f'''\

                        SFP{i}ptr = wmhSFP{i} + offset;'''
                s5 += f' + x{i} * (*SFP{i}ptr++)'
            s += f'''\

            case {i + 1}:{s1}
                while (t_v != t_vEnd)
                {{
                    {s2}
                    {s3})
                    {{
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        offset = nS * (*t_o);
                        {s4}
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += ({s5});
                    }}
                    t_v++;
                    t_o++;
                }}
                break;'''

        s += '''\

        }
    }\n\n'''
        return s

    def add_isotropic_compartments() -> str:
        s: str = '''\
    // isotropic compartments
    if (nISO > 0)
    {
        t_v = ISOv + ISOthreads[id];
        t_vEnd = ISOv + ISOthreads[id+1];
        xPtr0 = x + nIC*nF + nEC*nE + ISOthreads[id];
        switch (nISO)
        {'''
        s1: str = ''
        s2: str = 'x0 = *xPtr0++;'
        s3: str = 'if (x0 != 0'
        s4: str = 'SFP0ptr = isoSFP0;'
        s5: str = 'x0 * (*SFP0ptr++)'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

                xPtr{i} = xPtr{i - 1} + nV;'''
                s2 += f'''\

                    x{i} = *xPtr{i}++;'''
                s3 += f' || x{i} != 0'
                s4 += f'''\

                        SFP{i}ptr = isoSFP{i};'''
                s5 += f' + x{i} * (*SFP{i}ptr++)'
            s += f'''\

            case {i + 1}:{s1}
                while (t_v != t_vEnd)
                {{
                    {s2}
                    {s3})
                    {{
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        {s4}
                        while (YPtr != YPtrEnd)
                            (*YPtr++) += ({s5});
                    }}
                    t_v++;
                }}
                break;'''

        s += '''\

        }
    }\n\n'''
        return s

    s: str = '''\
//
// Compute a sub-block of the A*x MATRIX-VECTOR product
// 
void* COMMIT_A__block( void *ptr )
{
    int      id = (long)ptr;
    int      offset;
    uint32_t *eval0, *eval1, *eval2, *eval3, *eval4, *eval5, *eval6, *eval7, *eval8, *eval9, *eval10, *eval11, *eval12, *eval13, *eval14, *eval15, *eval16, *eval17, *eval18, *eval19;
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, w;
    double   *xPtr0, *xPtr1, *xPtr2, *xPtr3, *xPtr4, *xPtr5, *xPtr6, *xPtr7, *xPtr8, *xPtr9, *xPtr10, *xPtr11, *xPtr12, *xPtr13, *xPtr14, *xPtr15, *xPtr16, *xPtr17, *xPtr18, *xPtr19;
    double   *YPtr, *YPtrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr, *SFP10ptr, *SFP11ptr, *SFP12ptr, *SFP13ptr, *SFP14ptr, *SFP15ptr, *SFP16ptr, *SFP17ptr, *SFP18ptr, *SFP19ptr;
    uint32_t *t_v, *t_vEnd, *t_f;
    uint16_t *t_o;
    float    *t_l;\n\n'''
    s += add_intracellular_compartments()
    s += add_extracellular_compartments()
    s += add_isotropic_compartments()
    s += '''\
    pthread_exit( 0 );
}\n\n'''
    return s


def add_commit_a() -> str:
    def add_intracellular_compartments() -> str:
        s: str = '''\
    switch (nIC)
    {'''
        s1: str = 'wmrSFP0 = _wmrSFP;'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

            wmrSFP{i} = wmrSFP{i - 1} + _ndirs*_nS;'''
            s += f'''\

        case {i + 1}:
            {s1}
            break;'''
        s += '''\

    }\n\n'''
        return s

    def add_extracellular_compartments() -> str:
        s: str = '''\
    switch (nEC)
    {'''
        s1: str = 'wmhSFP0 = _wmhSFP;'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

            wmhSFP{i} = wmhSFP{i - 1} + _ndirs*_nS;'''
            s += f'''\

        case {i + 1}:
            {s1}
            break;'''
        s += '''\

    }\n\n'''
        return s

    def add_isotropic_compartments() -> str:
        s: str = '''\
    switch (nISO)
    {'''
        s1: str = 'isoSFP0 = _isoSFP;'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

            isoSFP{i} = isoSFP{i - 1} + _nS;'''
            s += f'''\

        case {i + 1}:
            {s1}
            break;'''
        s += '''\

    }\n\n'''
        return s

    s: str = '''\
//
// Function called by Cython
//
void COMMIT_A(
    int _nF, int _nE, int _nV, int _nS, int _ndirs,
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
    nISO = _nISO;\n\n'''
    s += add_intracellular_compartments()
    s += add_extracellular_compartments()
    s += add_isotropic_compartments()
    s += '''\

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
}\n\n'''
    return s


def add_commit_at_block() -> str:
    def add_intracellular_compartments() -> str:
        s: str = '''\
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
        {'''
        s1: str = 'SFP0ptr = wmrSFP0 + offset;'
        s2: str = '''\
x0 = (*SFP0ptr++) * YTmp;
                        eval0 = ICeval + *t_f;    
'''
        s3: str = 'x0 += (*SFP0ptr++) * YTmp;'
        s4: str = 'x[*t_f] += w * x0 * (double)(*eval0);'
        s5: str = ''

        for i in range(0, 20):
            if i > 0:
                if i > 1:
                    s5 = f'{i}*'
                s1 += f'''\

                        SFP{i}ptr = wmrSFP{i} + offset;'''
                s2 += f'''\

                        x{i} = (*SFP{i}ptr++) * YTmp;
                        eval{i} = eval{i - 1} + nF;'''
                s3 += f'''\

                        x{i} += (*SFP{i}ptr++) * YTmp;'''
                s4 += f'''\

                        x[*t_f+{s5}nF] += w * x{i} * (double)(*eval{i});'''
            s += f'''\

            case {i + 1}:
                while (t_v != t_vEnd)
                {{
                    if (*t_t == id)
                    {{
                        YPtr = Y + nS * (*t_v);
                        YPtrEnd = YPtr + nS;
                        w = (double)(*t_l);
                        offset = nS * (*t_o);
                        YTmp = *YPtr;
                        {s1}
                        {s2}
                        while (++YPtr != YPtrEnd)
                        {{
                            YTmp = *YPtr;
                            {s3}
                        }}
                        {s4}
                    }}
                    t_f++;
                    t_v++;
                    t_o++;
                    t_l++;
                    t_t++;
                }}
                break;'''

        s += '''\

        }
    }\n\n'''
        return s

    def add_extracellular_compartments() -> str:
        s: str = '''\
    // extra-cellular compartments
    if (nEC > 0)
    {
        t_v = ECv + ECthreadsT[id];
        t_vEnd = ECv + ECthreadsT[id+1];
        t_o = ECo + ECthreadsT[id];
        xPtr0 = x + nIC*nF + ECthreadsT[id];
        switch (nEC)
        {'''
        s1: str = ''
        s2: str = 'SFP0ptr = wmhSFP0 + offset;'
        s3: str = 'x0 = (*SFP0ptr++) * YTmp;'
        s4: str = 'x0 += (*SFP0ptr++) * YTmp;'
        s5: str = '(*xPtr0++) += x0;'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

                xPtr{i} = xPtr{i - 1} + nE;'''
                s2 += f'''\

                    SFP{i}ptr = wmhSFP{i} + offset;'''
                s3 += f'''\

                    x{i} = (*SFP{i}ptr++) * YTmp;'''
                s4 += f'''\

                    x{i} += (*SFP{i}ptr++) * YTmp;'''
                s5 += f'''\

                    (*xPtr{i}++) += x{i};'''
            s += f'''\

            case {i + 1}:{s1}
                while (t_v != t_vEnd)
                {{
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    offset = nS * (*t_o);
                    YTmp = *YPtr;
                    {s2}
                    {s3}
                    while (++YPtr != YPtrEnd)
                    {{
                        YTmp = *YPtr;
                        {s4}
                    }}
                    {s5}
                    t_v++;
                    t_o++;
                }}
                break;'''

        s += '''\

        }
    }\n\n'''
        return s

    def add_isotropic_compartments() -> str:
        s: str = '''\
    // isotropic compartments
    if (nISO > 0)
    {
        t_v = ISOv + ISOthreadsT[id];
        t_vEnd = ISOv + ISOthreadsT[id+1];
        xPtr0 = x + nIC*nF + nEC*nE + ISOthreadsT[id];
        switch (nISO)
        {'''
        s1: str = ''
        s2: str = 'SFP0ptr = isoSFP0;'
        s3: str = 'x0 = (*SFP0ptr++) * YTmp;'
        s4: str = 'x0 += (*SFP0ptr++) * YTmp;'
        s5: str = '(*xPtr0++) += x0;'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

                xPtr{i} = xPtr{i - 1} + nV;'''
                s2 += f'''\

                    SFP{i}ptr = isoSFP{i};'''
                s3 += f'''\

                    x{i} = (*SFP{i}ptr++) * YTmp;'''
                s4 += f'''\

                    x{i} += (*SFP{i}ptr++) * YTmp;'''
                s5 += f'''\

                    (*xPtr{i}++) += x{i};'''
            s += f'''\

            case {i + 1}:{s1}
                while (t_v != t_vEnd)
                {{
                    YPtr = Y + nS * (*t_v);
                    YPtrEnd = YPtr + nS;
                    YTmp = *YPtr;
                    {s2}
                    {s3}
                    while (++YPtr != YPtrEnd)
                    {{
                        YTmp = *YPtr;
                        {s4}
                    }}
                    {s5}
                    t_v++;
                }}
                break;'''

        s += '''\

        }
    }\n\n'''
        return s

    s: str = '''\
//
// Compute a sub-block of the At*y MATRIX-VECTOR product
//
void* COMMIT_At__block( void *ptr )
{
    int      id = (long)ptr;
    int      offset;
    uint32_t *eval0, *eval1, *eval2, *eval3, *eval4, *eval5, *eval6, *eval7, *eval8, *eval9, *eval10, *eval11, *eval12, *eval13, *eval14, *eval15, *eval16, *eval17, *eval18, *eval19;
    double   x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, w, YTmp;
    double   *xPtr0, *xPtr1, *xPtr2, *xPtr3, *xPtr4, *xPtr5, *xPtr6, *xPtr7, *xPtr8, *xPtr9, *xPtr10, *xPtr11, *xPtr12, *xPtr13, *xPtr14, *xPtr15, *xPtr16, *xPtr17, *xPtr18, *xPtr19;
    double   *YPtr, *YPtrEnd;
    float    *SFP0ptr, *SFP1ptr, *SFP2ptr, *SFP3ptr, *SFP4ptr, *SFP5ptr, *SFP6ptr, *SFP7ptr, *SFP8ptr, *SFP9ptr, *SFP10ptr, *SFP11ptr, *SFP12ptr, *SFP13ptr, *SFP14ptr, *SFP15ptr, *SFP16ptr, *SFP17ptr, *SFP18ptr, *SFP19ptr;
    uint32_t *t_v, *t_vEnd, *t_f;
    uint16_t *t_o;
    float    *t_l;
    uint8_t  *t_t;\n\n'''
    s += add_intracellular_compartments()
    s += add_extracellular_compartments()
    s += add_isotropic_compartments()
    s += '''\
    pthread_exit( 0 );
}\n\n'''
    return s


def add_commit_at() -> str:
    def add_intracellular_compartments() -> str:
        s: str = '''\
    switch (nIC)
    {'''
        s1: str = 'wmrSFP0 = _wmrSFP;'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

            wmrSFP{i} = wmrSFP{i - 1} + _ndirs*_nS;'''
            s += f'''\

        case {i + 1}:
            {s1}
            break;'''
        s += '''\

    }\n\n'''
        return s

    def add_extracellular_compartments() -> str:
        s: str = '''\
    switch (nEC)
    {'''
        s1: str = 'wmhSFP0 = _wmhSFP;'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

            wmhSFP{i} = wmhSFP{i - 1} + _ndirs*_nS;'''
            s += f'''\

        case {i + 1}:
            {s1}
            break;'''
        s += '''\

    }\n\n'''
        return s

    def add_isotropic_compartments() -> str:
        s: str = '''\
    switch (nISO)
    {'''
        s1: str = 'isoSFP0 = _isoSFP;'

        for i in range(0, 20):
            if i > 0:
                s1 += f'''\

            isoSFP{i} = isoSFP{i - 1} + _nS;'''
            s += f'''\

        case {i + 1}:
            {s1}
            break;'''
        s += '''\

    }\n\n'''
        return s

    s: str = '''\
//
// Function called by Cython
//
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
    nISO = _nISO;\n\n'''
    s += add_intracellular_compartments()
    s += add_extracellular_compartments()
    s += add_isotropic_compartments()
    s += '''\
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
}\n\n'''
    return s


def add_commit_a_block_nolut() -> str:
    def add_intracellular_compartments() -> str:
        s: str = '''\
    // intra-cellular compartments
    t_v    = ICv + ICthreads[id];
    t_vEnd = ICv + ICthreads[id+1];
    t_l    = ICl + ICthreads[id];
    t_f    = ICf + ICthreads[id];

    while( t_v != t_vEnd )
    {
        eval0 = ICeval[*t_f];
        x0 = x[*t_f]* (double)(*eval0);
        if ( x0 != 0 )
            Y[*t_v] += (double)(*t_l) * x0;
        t_f++;
        t_v++;
        t_l++;
    }\n\n'''
        return s

    def add_isotropic_compartments() -> str:
        s: str = '''\
    // isotropic compartments
    if (nISO > 0)
    {
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
    }\n\n'''
        return s

    s: str = '''\
//
// Compute a sub-block of the A*x MATRIX-VECTOR product
//
void* COMMIT_A__block_nolut( void *ptr )
{
    int      id = (long)ptr;
    uint32_t *eval0;
    double   x0;
    double   *xPtr;
    uint32_t *t_v, *t_vEnd, *t_f;
    float    *t_l;\n\n'''
    s += add_intracellular_compartments()
    s += add_isotropic_compartments()
    s += '''\
    pthread_exit( 0 );
}\n\n'''
    return s


def add_commit_a_nolut() -> str:
    s: str = '''\
//
// Function called by Cython
//
void COMMIT_A_nolut(
    int _nF,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICeval, uint32_t *_ICv, float *_ICl,
    uint32_t *_ISOv,
    uint32_t* _ICthreads, uint32_t* _ISOthreads,
    uint32_t _nISO, uint32_t _nThreads
)
{
    nF = _nF;

    x = _vIN;
    Y = _vOUT;

    ICf  = _ICf;
    ICeval = _ICeval;
    ICv  = _ICv;
    ICl  = _ICl;
    ISOv = _ISOv;

    nISO = _nISO;

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
}\n\n'''
    return s


def add_commit_at_block_nolut() -> str:
    def add_intracellular_compartments() -> str:
        s: str = '''\
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
            eval0 = ICeval + *t_f;
            x[*t_f] += (double)(*t_l) * Y[*t_v]* (double)(*eval0);
        t_t++;
        t_f++;
        t_v++;
        t_l++;
    }\n\n'''
        return s

    def add_isotropic_compartments() -> str:
        s: str = '''\
    // isotropic compartments
    if (nISO > 0)
    {
        t_v    = ISOv + ISOthreadsT[id];
        t_vEnd = ISOv + ISOthreadsT[id+1];
        xPtr   = x + nF + ISOthreadsT[id];

        while( t_v != t_vEnd )
            (*xPtr++) += Y[*t_v++];
    }\n\n'''
        return s

    s: str = '''\
//
// Compute a sub-block of the At*y MATRIX-VECTOR product
//
void* COMMIT_At__block_nolut( void *ptr )
{
    int      id = (long)ptr;
    double   *xPtr;
    uint32_t *eval0;
    uint32_t *t_v, *t_vEnd, *t_f;
    float    *t_l;
    uint8_t  *t_t;\n\n'''
    s += add_intracellular_compartments()
    s += add_isotropic_compartments()
    s += '''\

    pthread_exit( 0 );
}\n\n'''
    return s


def add_commit_at_nolut() -> str:
    s: str = '''\
//
// Function called by Cython
//
void COMMIT_At_nolut(
    int _nF, int _n,
    double *_vIN, double *_vOUT,
    uint32_t *_ICf, uint32_t *_ICeval, uint32_t *_ICv, float *_ICl,
    uint32_t *_ISOv,
    uint8_t* _ICthreadsT, uint32_t* _ISOthreadsT,
    uint32_t _nISO, uint32_t _nThreads
)
{
    nF = _nF;
    n  = _n;

    x = _vOUT;
    Y = _vIN;

    ICf  = _ICf;
    ICeval = _ICeval;
    ICv  = _ICv;
    ICl  = _ICl;
    ISOv = _ISOv;

    nISO = _nISO;

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
}\n'''
    return s


def add_commit() -> str:
    s: str = ''
    s += add_commit_a_block()
    s += add_commit_a()
    s += add_commit_at_block()
    s += add_commit_at()
    return s


def add_commit_nolut() -> str:
    s: str = ''
    s += add_commit_a_block_nolut()
    s += add_commit_a_nolut()
    s += add_commit_at_block_nolut()
    s += add_commit_at_nolut()
    return s


def write_operator_c_file() -> NoReturn:
    s: str = ''
    s += add_imports()
    s += add_global_variables()
    s += add_commit()
    s += add_commit_nolut()

    with open(path_join('commit', 'operator', 'operator_c.c'), 'w') as f:
        f.write(s)