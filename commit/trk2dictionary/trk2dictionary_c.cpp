#include <stdio.h>
#include <cstdio>
#include <string>
#include <map>
#include <vector>
#include "Vector.h"
#include "ProgressBar.h"
#include <numpy/arrayobject.h>
#include <math.h>
#include <iostream>
#include <thread>
#include <numeric>
#include <chrono>

#define _FILE_OFFSET_BITS 64
#define MAX_FIB_LEN 10000
#define MAX_THREADS 255

using namespace std;
ProgressBar* PROGRESS = new ProgressBar();

// CLASS to store the segments of one fiber
class segKey
{
    public:
    unsigned short x, y, z;
    unsigned short o;
    segKey(){}

    void set(unsigned short _x, unsigned short _y, unsigned short _z, unsigned short _o)
    {
        x  = _x;
        y  = _y;
        z  = _z;
        o = _o;
    }

    bool const operator <(const segKey& seg) const
    {
        return o < seg.o || (o==seg.o && z<seg.z) || (o==seg.o && z==seg.z && y<seg.y) || (o==seg.o && z==seg.z && y==seg.y && x<seg.x);
    }
};

class segInVoxKey
{
    public:
    unsigned short x, y, z;
    segInVoxKey(){}

    void set(unsigned short _x, unsigned short _y, unsigned short _z)
    {
        x  = _x;
        y  = _y;
        z  = _z;
    }
    bool const operator <(const segInVoxKey& o) const
    {
        return (z<o.z) || (z==o.z && y<o.y) || (z==o.z && y==o.y && x<o.x);
    }
};


// global variables (to avoid passing them at each call)
thread_local map<segKey,float>               FiberSegments;
thread_local float                           FiberLen;      // length of a streamline
thread_local float                           FiberLenTot;   // length of a streamline (considering the blur)
thread_local vector< Vector<double> >        P;

Vector<int>     dim;
Vector<float>   pixdim;
float*          ptrMASK;
float*          ptrISO;
float           fiberShiftXmm, fiberShiftYmm, fiberShiftZmm;
bool            doIntersect;
float           minSegLen, minFiberLen, maxFiberLen;
// Threads variables
vector<thread>  threads;
vector<unsigned int>    totICSegments; 
vector<unsigned int>    totFibers;
unsigned int            totECVoxels = 0;
unsigned int            totECSegments = 0;
bool                    skip_ec = false;

// progressbar verbosity
int verbosity = 0;

// --- Functions Definitions ----
bool rayBoxIntersection( Vector<double>& origin, Vector<double>& direction, Vector<double>& vmin, Vector<double>& vmax, double & t);
void fiberForwardModel( float fiber[3][MAX_FIB_LEN], unsigned int pts, int nReplicas, double* ptrBlurRho, double* ptrBlurAngle, double* ptrBlurWeights, bool doApplyBlur, short* ptrHashTable, float* ptrPEAKS, double* ptrPeaksAffine, int Np, float angle_thr, vector<Vector<double>>& P );
void segmentForwardModel( const Vector<double>& P1, const Vector<double>& P2, int k, double w, short* ptrHashTable, float* ptrPEAKS, double* ptrPeaksAffine, int Np, float angle_thr );
unsigned int read_fiberTRK( FILE* fp, float fiber[3][MAX_FIB_LEN], int ns, int np );
unsigned int read_fiberTCK( FILE* fp, float fiber[3][MAX_FIB_LEN] , float* toVOXMM );


// ---------- Parallel fuction --------------
int ICSegments( char* str_filename, int isTRK, int n_count, int nReplicas, int n_scalars, int n_properties, float* ptrToVOXMM,
double* ptrTDI, float angle_thr, double* ptrBlurRho, double* ptrBlurAngle, double* ptrBlurWeights, bool* ptrBlurApplyTo, short* ptrHashTable, float* ptrPEAKS, double* ptrPeaksAffine, int Np, char* path_out, 
unsigned long long int offset, int idx, unsigned int startpos, unsigned int endpos );

int ECSegments(float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
    double** ptrTDI, short* ptrHashTable, char* path_out, double* ptrPeaksAffine, int idx);

int ISOcompartments(double** ptrTDI, char* path_out, int idx);



// ===========================================
//          Called by Cython
// ===========================================

int trk2dictionary(
    char* str_filename, int data_offset, int Nx, int Ny, int Nz, float Px, float Py, float Pz, int n_count, int n_scalars, int n_properties,
    float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, float min_seg_len, float min_fiber_len, float max_fiber_len,
    float* ptrPEAKS, float angle_thr, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
    float* _ptrMASK, float* _ptrISO, double** ptrTDI, char* path_out, int c, double* ptrPeaksAffine,
    int nReplicas, double* ptrBlurRho, double* ptrBlurAngle, double* ptrBlurWeights, bool* ptrBlurApplyTo,
    float* ptrToVOXMM, short* ptrHashTable, int threads_count, int verbose
)
{

    // Set the global variables
    // ----------------------------------
    dim.Set( Nx, Ny, Nz );
    pixdim.Set( Px, Py, Pz );
    fiberShiftXmm = fiber_shiftX * pixdim.x; // shift in mm for the coordinates
    fiberShiftYmm = fiber_shiftY * pixdim.y;
    fiberShiftZmm = fiber_shiftZ * pixdim.z;
    ptrMASK       = _ptrMASK;
    ptrISO        = _ptrISO;
    doIntersect   = c > 0;
    minSegLen     = min_seg_len;
    minFiberLen   = min_fiber_len;
    maxFiberLen   = max_fiber_len;
    totICSegments.resize( threads_count, 0 );
    totFibers.resize( threads_count, 0 );
    totECVoxels   = 0;
    totECSegments = 0;
    verbosity     = verbose;


    // Compute the batch size for each thread
    // ---------------------------------------
    int p_size = n_count > threads_count ? threads_count : n_count;
    int* batch_size = new int[p_size]();

    for (int i = 0; i < n_count; i++)
    {
        batch_size[i%p_size] += 1;
    }

    
    // Compute the starting position
    // -----------------------------------------
    unsigned int elements = threads_count + 1;
    unsigned int *Pos = new unsigned int[elements]();

    for (int i = 1; i < threads_count; i++)
    {
        Pos[i] = Pos[i-1] + batch_size[i];
    }

    Pos[threads_count] = n_count;


    // Check the file extension
    // -------------------------------------
    int isTRK;

    char *ext = strrchr(str_filename, '.');
    if (strcmp(ext,".trk")==0) // for .trk file
        isTRK = 1;
    else if (strcmp(ext,".tck")==0) // for .tck file
        isTRK = 0;
        else
        return 0;    


    // Open tractogram file and compute the offset for each thread
    // -----------------------------------------------------------------
    unsigned long long int current;
    unsigned long long int OffsetArr[threads_count];
    int f = 0;
    float Buff[3];
    int N;

    FILE* fpTractogram = fopen(str_filename,"rb");
    if (fpTractogram == NULL) return 0;
    fseek( fpTractogram, data_offset, SEEK_SET ); // skip the header   

    OffsetArr[0] = ftell( fpTractogram );

    if(isTRK) {
        while( f != n_count) {
            fread( (char*)&N, 1, 4, fpTractogram ); // read the number of points in each streamline
            if( N >= MAX_FIB_LEN || N <= 0 ) return 0;   // check
            for( int k=0; k<N; k++){
                fread((char*)Buff, 1, 12, fpTractogram);
            }
            fseek(fpTractogram,4*n_properties,SEEK_CUR);
            f++;
            current = ftell( fpTractogram );
            for( int i = 1; i < threads_count; i++ ){
                    if( f == Pos[i] ) 
                        OffsetArr[i] = current;
                }
        }
    } else {

        while( f != n_count ) {

            fread((char*)Buff, 1, 12, fpTractogram );

            if( isnan(Buff[0]) ){
                f++;
                current = ftell( fpTractogram );

                for( int i = 1; i < threads_count; i++ ){
                    if( f == Pos[i] ) 
                        OffsetArr[i] = current;
                }

            }

        }

    }

    fclose(fpTractogram);



    // ==========================================
    //          Parallel IC compartments
    // ==========================================

    printf( "\n   \033[0;32m* Exporting IC compartments:\033[0m\n" );
    // unsigned int width = 25;
    // PROGRESS = new ProgressBar( (unsigned int) n_count, (unsigned int) width);
    if (verbosity > 0)
    {
        PROGRESS->reset((unsigned int) n_count);
        PROGRESS->setPrefix("     ");
    }
    // ---- Original ------
    for( int i = 0; i<threads_count; i++ ){
        threads.push_back( thread( ICSegments, str_filename, isTRK, n_count, nReplicas, n_scalars, n_properties, ptrToVOXMM,
        ptrTDI[i], angle_thr, ptrBlurRho, ptrBlurAngle, ptrBlurWeights, ptrBlurApplyTo, ptrHashTable, ptrPEAKS, ptrPeaksAffine, Np, path_out, OffsetArr[i], 
        i, Pos[i], Pos[i+1]  ) );
    }


    for( int i = 0; i<threads_count; i++ ) {
        threads[i].join();
    }

    if (verbosity > 0)
        PROGRESS->close();

    printf( "     [ %d streamlines kept, %d segments in total ]\n", std::accumulate(totFibers.begin(), totFibers.end(), 0), std::accumulate( totICSegments.begin(), totICSegments.end(), 0) );
    totFibers.clear();
    threads.clear();

    // ==========================================
    //          Parallel EC compartments
    // ==========================================

    printf( "\n   \033[0;32m* Exporting EC compartments:\033[0m\n" );

    int EC = ECSegments( ptrPEAKS, Np, vf_THR, ECix, ECiy, ECiz, ptrTDI, ptrHashTable, path_out, ptrPeaksAffine, threads_count );

    printf("     [ %d voxels, %d segments ]\n", totECVoxels, totECSegments );

    /*=========================*/
    /*     Restricted ISO compartments     */
    /*=========================*/
    printf( "\n   \033[0;32m* Exporting ISO compartments:\033[0m\n" );

    int totISO = ISOcompartments(ptrTDI, path_out, threads_count);

    printf("     [ %d voxels ]\n", totISO );

    return 1;

}


int ECSegments(float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
    double** ptrTDI, short* ptrHashTable, char* path_out, double* ptrPeaksAffine, int threads){

    // Variables definition
    string    filename;
    string    OUTPUT_path(path_out);
    std::size_t pos = OUTPUT_path.find("/temp");
    OUTPUT_path = OUTPUT_path.substr (0,pos);

    unsigned short o;
    unsigned int v;
    unsigned int temp_totECSegments = 0, temp_totECVoxels = 0;


    filename = OUTPUT_path+"/dictionary_EC_v.dict";        FILE* pDict_EC_v       = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_EC_o.dict";        FILE* pDict_EC_o       = fopen(filename.c_str(),"wb");

    if ( ptrPEAKS != NULL &&  skip_ec )
    {
        Vector<double> dir;
        double         longitude, colatitude;
        segKey         ec_seg;
        int            ix, iy, iz, id, atLeastOne;
        float          peakMax;
        float          norms[ Np ];
        float          *ptr;
        int            ox, oy;
        int            skip = 0;

        for(iz=0; iz<dim.z; iz++)
        {
            for(iy=0; iy<dim.y; iy++)
                for(ix=0; ix<dim.x; ix++)
                {
                    // check if in mask previously computed from IC segments
                    for(int i =0; i<threads; i++){
                        if ( ptrTDI[i][ iz + dim.z * ( iy + dim.y * ix ) ] == 0 ){
                            skip += 1;
                        }
                    }
                    if(skip==threads){
                        skip = 0;
                        continue;
                    }
                    skip = 0;
                    peakMax = -1;
                    for(id=0; id<Np ;id++)
                    {
                        ptr = ptrPEAKS + 3*(id + Np * ( iz + dim.z * ( iy + dim.y * ix ) ));
                        dir.x = ptr[0];
                        dir.y = ptr[1];
                        dir.z = ptr[2];
                        norms[id] = dir.norm();
                        if ( norms[id] > peakMax )
                            peakMax = norms[id];
                    }

                    if ( peakMax > 0 )
                    {
                        ec_seg.x  = ix;
                        ec_seg.y  = iy;
                        ec_seg.z  = iz;
                        atLeastOne = 0;
                        for(id=0; id<Np ;id++)
                        {
                        if ( norms[id]==0 || norms[id] < vf_THR*peakMax ) continue; // peak too small, don't consider it

                        // get the orientation of the current peak
                        ptr = ptrPEAKS + 3*(id + Np * ( iz + dim.z * ( iy + dim.y * ix ) ));

                        // multiply by the affine matrix
                        dir.x = ptr[0] * ptrPeaksAffine[0] + ptr[1] * ptrPeaksAffine[1] + ptr[2] * ptrPeaksAffine[2];
                        dir.y = ptr[0] * ptrPeaksAffine[3] + ptr[1] * ptrPeaksAffine[4] + ptr[2] * ptrPeaksAffine[5];
                        dir.z = ptr[0] * ptrPeaksAffine[6] + ptr[1] * ptrPeaksAffine[7] + ptr[2] * ptrPeaksAffine[8];

                        // flip axes if requested
                        dir.x *= ECix;
                        dir.y *= ECiy;
                        dir.z *= ECiz;
                        if ( dir.y < 0 )
                        {
                            // ensure to be in the right hemisphere (the one where kernels were pre-computed)
                            dir.x = -dir.x;
                            dir.y = -dir.y;
                            dir.z = -dir.z;
                        }
                        colatitude = atan2( sqrt(dir.x*dir.x + dir.y*dir.y), dir.z );
                        longitude  = atan2( dir.y, dir.x );
                        ox = (int)round(colatitude/M_PI*180.0);
                        oy = (int)round(longitude/M_PI*180.0);

                        v = ec_seg.x + dim.x * ( ec_seg.y + dim.y * ec_seg.z );
                        o = ptrHashTable[ox*181 + oy];
                        fwrite( &v, 4, 1, pDict_EC_v );
                        fwrite( &o, 2, 1, pDict_EC_o );
                        temp_totECSegments++;
                        atLeastOne = 1;
                        }
                    if ( atLeastOne>0 )
                        temp_totECVoxels++;
                    }
                }
        }
    }
    totECSegments = temp_totECSegments;
    totECVoxels = temp_totECVoxels;

    fclose( pDict_EC_v );
    fclose( pDict_EC_o );

    return 1;

}



int ISOcompartments(double** ptrTDI, char* path_out, int threads){
    // Variables definition
    string    filename;
    string    OUTPUT_path(path_out);
    std::size_t pos = OUTPUT_path.find("/temp");
    OUTPUT_path = OUTPUT_path.substr (0,pos);
    unsigned int totISOVoxels = 0, v=0;

    filename = OUTPUT_path+"/dictionary_ISO_v.dict";        FILE* pDict_ISO_v   = fopen( filename.c_str(),   "wb" );

    int            ix, iy, iz, id, atLeastOne;
    int            skip = 0;

    for(iz=0; iz<dim.z ;iz++){
        for(iy=0; iy<dim.y ;iy++)
        for(ix=0; ix<dim.x ;ix++){
            // check if ptrISO and ptrMASK are not NULL
            if ( ptrISO != NULL ){
                if ( ptrISO[ iz + dim.z * ( iy + dim.y * ix ) ] == 0 ) continue;
            }
            if ( ptrMASK != NULL ){
                if ( ptrMASK[ iz + dim.z * ( iy + dim.y * ix ) ] == 0 ) continue;
            }
            // check if in mask previously computed from IC segments
            for(int i =0; i<threads; i++){
                if ( ptrTDI[i][ iz + dim.z * ( iy + dim.y * ix ) ] == 0 ){
                    skip += 1;
                }
            }
            if(skip==threads){
                skip = 0;
                continue;
            }
            skip = 0;
            v = ix + dim.x * ( iy + dim.y * iz );
            fwrite( &v, 4, 1, pDict_ISO_v );    
            totISOVoxels++; 
        }
    }
    fclose( pDict_ISO_v );

    return totISOVoxels;
}




/********************************************************************************************************************/
/*                                                Parallel Function                                                 */
/********************************************************************************************************************/

int ICSegments( char* str_filename, int isTRK, int n_count, int nReplicas, int n_scalars, int n_properties, float* ptrToVOXMM, double* ptrTDI, float angle_thr, double* ptrBlurRho, 
double* ptrBlurAngle, double* ptrBlurWeights, bool* ptrBlurApplyTo, short* ptrHashTable, float* ptrPEAKS, double* ptrPeaksAffine, int Np, char* path_out, 
unsigned long long int offset, int idx, unsigned int startpos, unsigned int endpos ) 
{

    // Variables definition
    float           fiber[3][MAX_FIB_LEN] = {0} ;
    float           fiberNorm;   // normalization
    unsigned int    N, v, tempTotFibers, temp_totICSegments;
    unsigned short  o;
    unsigned char   kept;
    string          filename;
    string          OUTPUT_path(path_out);

    unsigned int sumFibers = startpos;

    map<segKey,float>::iterator it;
    map<segInVoxKey,float> FiberNorm;
    map<segInVoxKey,float>::iterator itNorm;

    segInVoxKey inVoxKey; 

    P.resize(nReplicas);


    // Creo e apro i files per i risultati
    filename = OUTPUT_path+"/dictionary_TRK_norm_" + std::to_string(idx) + ".dict";   FILE* pDict_TRK_norm = fopen(filename.c_str(),"wb");
    if ( !pDict_TRK_norm )
    {
        printf( "\n[trk2dictionary] Unable to create output files" );
    }
    filename = OUTPUT_path+"/dictionary_IC_f_" + std::to_string(idx) + ".dict";        FILE* pDict_IC_f       = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_v_" + std::to_string(idx) + ".dict";        FILE* pDict_IC_v       = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_o_" + std::to_string(idx) + ".dict";        FILE* pDict_IC_o       = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_len_" + std::to_string(idx) + ".dict";      FILE* pDict_IC_len     = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_TRK_len_" + std::to_string(idx) + ".dict";     FILE* pDict_TRK_len    = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_TRK_lenTot_" + std::to_string(idx) + ".dict";  FILE* pDict_TRK_lenTot = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_TRK_kept_" + std::to_string(idx) + ".dict";    FILE* pDict_TRK_kept   = fopen(filename.c_str(),"wb");


    // ---- Original -----
    // Open tractogram file
    FILE* fpTractogram1 = fopen( str_filename,"rb" );
    if ( fpTractogram1 == NULL ) return 0; // if there's no tractogram file, then return 0
    fseek(fpTractogram1, offset, SEEK_SET);

    tempTotFibers = 0;
    temp_totICSegments = 0;
    int incr_new = 0;
    int incr_old = 0;
    // Iterate over streamlines

    for(int f=startpos; f<endpos; f++) 
    {        

        if ( isTRK )
            N = read_fiberTRK( fpTractogram1, fiber, n_scalars, n_properties );
        else
            N = read_fiberTCK( fpTractogram1, fiber , ptrToVOXMM );


        fiberForwardModel( fiber, N, nReplicas, ptrBlurRho, ptrBlurAngle, ptrBlurWeights, ptrBlurApplyTo[f], ptrHashTable, ptrPEAKS, ptrPeaksAffine, Np, angle_thr, P );

        kept = 0;

        if ( FiberSegments.size() > 0 )
        {
            if ( FiberLen > minFiberLen && FiberLen < maxFiberLen )
            {
                // add segments to files
                for (it=FiberSegments.begin(); it!=FiberSegments.end(); it++)
                {
                    // NB: plese note inverted ordering for 'v'
                    v = it->first.x + dim.x * ( it->first.y + dim.y * it->first.z );
                    o = it->first.o;       

                    fwrite( &sumFibers,      4, 1, pDict_IC_f );
                    fwrite( &v,              4, 1, pDict_IC_v );
                    fwrite( &o,              2, 1, pDict_IC_o );
                    fwrite( &(it->second),   4, 1, pDict_IC_len );       
                    
                    ptrTDI[ it->first.z + dim.z * ( it->first.y + dim.y * it->first.x ) ] += it->second;

                    inVoxKey.set( it->first.x, it->first.y, it->first.z );
                    FiberNorm[inVoxKey] += it->second;
                }

                for (fiberNorm=0, itNorm=FiberNorm.begin(); itNorm!=FiberNorm.end(); itNorm++)
                    // std::cout << "FiberNorm: " << itNorm->second << std::endl;
                    fiberNorm += pow(itNorm->second,2);
                fiberNorm = sqrt(fiberNorm);
                FiberNorm.clear();

                fwrite( &fiberNorm,   1, 4, pDict_TRK_norm );   // actual length considered in optimization
                fwrite( &FiberLen,    1, 4, pDict_TRK_len );    // length of the streamline
                fwrite( &FiberLenTot, 1, 4, pDict_TRK_lenTot ); // length of the streamline (considering the blur)

                temp_totICSegments += FiberSegments.size();
                sumFibers ++;
                tempTotFibers ++;
                kept = 1;
            }

        }

        fwrite( &kept, 1, 1, pDict_TRK_kept );
        totFibers[idx] = tempTotFibers;
        totICSegments[idx] = temp_totICSegments;

        if (verbosity > 0)
        {
            if (idx == 0)
            {
                incr_new = std::accumulate(totFibers.begin(), totFibers.end(), 0);
                for(int i=incr_old; i<incr_new; i++)
                    PROGRESS->inc();
                incr_old = incr_new;
            }
        }
    }
    fclose( fpTractogram1 );
    fclose( pDict_TRK_norm );
    fclose( pDict_IC_f );
    fclose( pDict_IC_v );
    fclose( pDict_IC_o );
    fclose( pDict_IC_len );
    fclose( pDict_TRK_len );
    fclose( pDict_TRK_lenTot );
    fclose( pDict_TRK_kept );

    return 1;
}


/********************************************************************************************************************/
/*                                                 fiberForwardModel                                                */
/********************************************************************************************************************/
void fiberForwardModel( float fiber[3][MAX_FIB_LEN], unsigned int pts, int nReplicas, double* ptrBlurRho, double* ptrBlurAngle, double* ptrBlurWeights, bool doApplyBlur, short* ptrHashTable,
                        float* ptrPEAKS, double* ptrPeaksAffine, int Np, float angle_thr, vector<Vector<double>>& P )
{
    thread_local static Vector<double> S1, S2, S1m, S2m, P_old, P_int, q, n, nr, qxn, qxqxn;
    thread_local static Vector<double> vox, vmin, vmax, dir1, dir2;
    thread_local static double         len, t, alpha, w, R, dot;
    thread_local static int            i, j, k;

    FiberLen = 0.0;
    FiberLenTot = 0.0;
    FiberSegments.clear();
    if ( pts <= 2 )
        return;

    // create duplicate points on circles
    S1.x = fiber[0][0]+fiberShiftXmm;
    S1.y = fiber[1][0]+fiberShiftYmm;
    S1.z = fiber[2][0]+fiberShiftZmm;
    dir2.x = (fiber[0][1]+fiberShiftXmm) - S1.x;
    dir2.y = (fiber[1][1]+fiberShiftYmm) - S1.y;
    dir2.z = (fiber[2][1]+fiberShiftZmm) - S1.z;
    dir2.Normalize();
    n.x = dir2.y-dir2.z;
    n.y = dir2.z-dir2.x;
    n.z = dir2.x-dir2.y;
    n.Normalize();

    // duplicate first point and move to corresponding grid locations
    for(k=0; k<nReplicas ;k++)
    {
        if ( !doApplyBlur && k>0 )
            continue;
        R = ptrBlurRho[k];
        alpha = ptrBlurAngle[k];

        // quaternion (q.x, q.y, q.z, w) for rotation
        w = sin(alpha/2.0);
        q.x = dir2.x * w;
        q.y = dir2.y * w;
        q.z = dir2.z * w;
        w = cos(alpha/2.0);

        // rotate the segment's normal
        qxn.x = 2.0 * ( q.y * n.z - q.z * n.y );
        qxn.y = 2.0 * ( q.z * n.x - q.x * n.z );
        qxn.z = 2.0 * ( q.x * n.y - q.y * n.x );
        qxqxn.x = q.y * qxn.z - q.z * qxn.y;
        qxqxn.y = q.z * qxn.x - q.x * qxn.z;
        qxqxn.z = q.x * qxn.y - q.y * qxn.x;
        nr.x = n.x + w * qxn.x + qxqxn.x;
        nr.y = n.y + w * qxn.y + qxqxn.y;
        nr.z = n.z + w * qxn.z + qxqxn.z;
        nr.Normalize();

        // move first point to corresponding grid location
        S2.x = S1.x + R*nr.x;
        S2.y = S1.y + R*nr.y;
        S2.z = S1.z + R*nr.z;
        P[k] = S2;
    }

    // move all remaining points
    for(i=1; i<pts ;i++)
    {
        /* get the intersection plane */
        // S2 = point on plane
        S2.x = fiber[0][i]+fiberShiftXmm;
        S2.y = fiber[1][i]+fiberShiftYmm;
        S2.z = fiber[2][i]+fiberShiftZmm;

        // n = normal to plane
        dir1.x = S2.x - (fiber[0][i-1]+fiberShiftXmm);
        dir1.y = S2.y - (fiber[1][i-1]+fiberShiftYmm);
        dir1.z = S2.z - (fiber[2][i-1]+fiberShiftZmm);
        dir1.Normalize();
        if ( i == pts-1 )
        {
            dir2.x = dir1.x;
            dir2.y = dir1.y;
            dir2.z = dir1.z;
        } else {
            dir2.x = (fiber[0][i+1]+fiberShiftXmm) - S2.x;
            dir2.y = (fiber[1][i+1]+fiberShiftYmm) - S2.y;
            dir2.z = (fiber[2][i+1]+fiberShiftZmm) - S2.z;
            dir2.Normalize();
        }
        n.x = 0.5*(dir1.x+dir2.x);
        n.y = 0.5*(dir1.y+dir2.y);
        n.z = 0.5*(dir1.z+dir2.z);

        // normalize to avoid computations later on
        dot = dir1.x*n.x + dir1.y*n.y + dir1.z*n.z;
        n.x /= dot;
        n.y /= dot;
        n.z /= dot;

        /* translate points */
        for(k=0; k<nReplicas ;k++)
        {
            if ( !doApplyBlur && k>0 )
                continue;

            if ( ptrBlurWeights[k] < 1e-3 )
                continue;

            P_old.x = P[k].x;
            P_old.y = P[k].y;
            P_old.z = P[k].z;
            len = (S2.x-P_old.x)*n.x + (S2.y-P_old.y)*n.y + (S2.z-P_old.z)*n.z;
            if ( len>0 )
            {
                P[k].x += dir1.x * len;
                P[k].y += dir1.y * len;
                P[k].z += dir1.z * len;

                /* save segment */
                if ( doIntersect==false )
                    segmentForwardModel( P_old, P[k], k, ptrBlurWeights[k], ptrHashTable, ptrPEAKS, ptrPeaksAffine, Np, angle_thr );
                else
                {
                    S1m.x = P_old.x;
                    S1m.y = P_old.y;
                    S1m.z = P_old.z;
                    S2m.x = P[k].x;
                    S2m.y = P[k].y;
                    S2m.z = P[k].z;
                    while( 1 )
                    {
                        len = sqrt( pow(S2m.x-S1m.x,2) + pow(S2m.y-S1m.y,2) + pow(S2m.z-S1m.z,2) ); // in mm
                        if ( len <= minSegLen )
                            break;

                        if ( floor(S1m.x/pixdim.x)==floor(S2m.x/pixdim.x) &&
                            floor(S1m.y/pixdim.y)==floor(S2m.y/pixdim.y) &&
                            floor(S1m.z/pixdim.z)==floor(S2m.z/pixdim.z)
                            )
                        {
                            // same voxel, no need to compute intersections
                            segmentForwardModel( S1m, S2m, k, ptrBlurWeights[k], ptrHashTable, ptrPEAKS, ptrPeaksAffine, Np, angle_thr);
                            break;
                        }

                        // compute AABB of the first point (in mm)
                        vmin.x = floor( (S1m.x + 1e-6*dir1.x)/pixdim.x ) * pixdim.x;
                        vmin.y = floor( (S1m.y + 1e-6*dir1.y)/pixdim.y ) * pixdim.y;
                        vmin.z = floor( (S1m.z + 1e-6*dir1.z)/pixdim.z ) * pixdim.z;
                        vmax.x = vmin.x + pixdim.x;
                        vmax.y = vmin.y + pixdim.y;
                        vmax.z = vmin.z + pixdim.z;

                        if ( rayBoxIntersection( S1m, dir1, vmin, vmax, t ) && t>0 && t<len )
                        {
                            // add the portion S1P, and then reiterate
                            P_int.x = S1m.x + t*dir1.x;
                            P_int.y = S1m.y + t*dir1.y;
                            P_int.z = S1m.z + t*dir1.z;
                            segmentForwardModel( S1m, P_int, k, ptrBlurWeights[k], ptrHashTable, ptrPEAKS, ptrPeaksAffine, Np, angle_thr );
                            S1m.x = P_int.x;
                            S1m.y = P_int.y;
                            S1m.z = P_int.z;
                        }
                        else
                        {
                            // add the segment S1S2 and stop iterating
                            segmentForwardModel( S1m, S2m, k, ptrBlurWeights[k], ptrHashTable, ptrPEAKS, ptrPeaksAffine, Np, angle_thr );
                            break;
                        }
                    }
                }
            }
        }
    }
}


/********************************************************************************************************************/
/*                                                segmentForwardModel                                               */
/********************************************************************************************************************/
void segmentForwardModel( const Vector<double>& P1, const Vector<double>& P2, int k, double w, short* ptrHashTable, float* ptrPEAKS, double* ptrPeaksAffine, int Np, float angle_thr )
{
    thread_local static Vector<int>    vox;
    thread_local static Vector<double> dir, dirpeak;
    thread_local static double         longitude, colatitude, len;
    thread_local static segKey         key;
    thread_local static int            ox, oy;
    thread_local static int            skip = 0;
    thread_local static float          *ptr_peak;
    thread_local static int            id;
    thread_local static float          peakMax;
    thread_local static float          norm, angle;
    thread_local static unsigned int   tempt;
    thread_local static float          dot_angle=0.0;
    thread_local static int            id_peak=0;


    // direction of the segment
    dir.y = P2.y-P1.y;
    if ( dir.y >= 0 )
    {
        dir.x = P2.x-P1.x;
        dir.z = P2.z-P1.z;
    }
    else
    {
        dir.x = P1.x-P2.x;
        dir.y = P1.y-P2.y;
        dir.z = P1.z-P2.z;
    }

    // length of the segment
    len = dir.norm();
    if ( w*len <= minSegLen )
        return;
    dir.Normalize();

    // voxel of the segment is the centroid
    vox.x = floor( 0.5 * (P1.x + P2.x) / pixdim.x );
    vox.y = floor( 0.5 * (P1.y + P2.y) / pixdim.y );
    vox.z = floor( 0.5 * (P1.z + P2.z) / pixdim.z );

    if ( vox.x>=dim.x || vox.x<0 || vox.y>=dim.y || vox.y<0 || vox.z>=dim.z || vox.z<0 )
        return;
    if ( ptrMASK && ptrMASK[ vox.z + dim.z * ( vox.y + dim.y * vox.x ) ]==0 )
        return;


    // get the orientation of the current peak
    // ptr = ptrPEAKS + 3*(id + Np * ( iz + dim.z * ( iy + dim.y * ix ) ));
    tempt = 0;
    peakMax = -1;
    for(id=0; id<Np ;id++)
    {   
        ptr_peak = ptrPEAKS + 3*(id + Np * ( vox.z + dim.z * ( vox.y + dim.y * vox.x ) ));
        dirpeak.x = ptr_peak[0];
        dirpeak.y = ptr_peak[1];
        dirpeak.z = ptr_peak[2];
        norm = dirpeak.norm();

        if (norm > peakMax)
            peakMax = norm;
            id_peak = id;
    }
        
    for (id=0; id<Np ;id++)
    {
        ptr_peak = ptrPEAKS + 3*(id + Np * ( vox.z + dim.z * ( vox.y + dim.y * vox.x ) ));
        if (ptr_peak[0] == 0 && ptr_peak[1] == 0 && ptr_peak[2] == 0)
        {
            break;
        }
        dirpeak.x = ptr_peak[0];
        dirpeak.y = ptr_peak[1];
        dirpeak.z = ptr_peak[2];
        norm = dirpeak.norm();
        if ( norm < 0.1*peakMax ){
            tempt++;
            continue;
        }

        dirpeak.Normalize();
    
        if ( dir.y < 0 )
        {
            dirpeak.x = -dirpeak.x;
            dirpeak.y = -dirpeak.y;
            dirpeak.z = -dirpeak.z;
        }

        dot_angle = dir.x*dirpeak.x + dir.y*dirpeak.y + dir.z*dirpeak.z;

        if ( dot_angle < 0 ){
            dot_angle = -dot_angle;
            angle = acos( dot_angle );
        }
        else {
            angle = acos( dot_angle );
        }

        if ( angle < angle_thr )
            // break;
            break;
        tempt++;
    }
    if (tempt == Np)
        return;
    

    // add the segment to the data structure
    longitude  = atan2(dir.y, dir.x);
    colatitude = atan2( sqrt(dir.x*dir.x + dir.y*dir.y), dir.z );
    ox = (int)round(colatitude/M_PI*180.0); // theta // i1
    oy = (int)round(longitude/M_PI*180.0);  // phi   // i2
    key.set( vox.x, vox.y, vox.z, (unsigned short) ptrHashTable[ox*181 + oy] );
    FiberSegments[key] += w*len;
    FiberLenTot += w*len;
    if ( k==0 ) // fiber length computed only from original segments
        FiberLen += len;
}


/********************************************************************************************************************/
/*                                                rayBoxIntersection                                                */
/********************************************************************************************************************/
bool rayBoxIntersection( Vector<double>& origin, Vector<double>& direction, Vector<double>& vmin, Vector<double>& vmax, double & t)
{
    thread_local static double tmin, tmax, tymin, tymax, tzmin, tzmax;
    thread_local static Vector<double> invrd;

    // inverse direction to catch float problems
    invrd.x = 1.0 / direction.x;
    invrd.y = 1.0 / direction.y;
    invrd.z = 1.0 / direction.z;

    if (invrd.x >= 0)
    {
      tmin = (vmin.x - origin.x) * invrd.x;
      tmax = (vmax.x - origin.x) * invrd.x;
    }
    else
    {
      tmin = (vmax.x - origin.x) * invrd.x;
      tmax = (vmin.x - origin.x) * invrd.x;
    }

    if (invrd.y >= 0)
    {
      tymin = (vmin.y - origin.y) * invrd.y;
      tymax = (vmax.y - origin.y) * invrd.y;
    }
    else
    {
      tymin = (vmax.y - origin.y) * invrd.y;
      tymax = (vmin.y - origin.y) * invrd.y;
    }

    if ( (tmin > tymax) || (tymin > tmax) ) return false;
    if ( tymin > tmin) tmin = tymin;
    if ( tymax < tmax) tmax = tymax;

    if (invrd.z >= 0)
    {
      tzmin = (vmin.z - origin.z) * invrd.z;
      tzmax = (vmax.z - origin.z) * invrd.z;
    }else
    {
      tzmin = (vmax.z - origin.z) * invrd.z;
      tzmax = (vmin.z - origin.z) * invrd.z;
    }

    if ( (tmin > tzmax) || (tzmin > tmax) ) return false;
    if ( tzmin > tmin) tmin = tzmin;
    if ( tzmax < tmax) tmax = tzmax;

    // check if values are valid
    t = tmin;
    if (t <= 0) t = tmax;

    return true;
}


// Read a fiber from file .trk
unsigned int read_fiberTRK( FILE* fp, float fiber[3][MAX_FIB_LEN], int ns, int np )
{
    int N;
    fread((char*)&N, 1, 4, fp);

    if ( N >= MAX_FIB_LEN || N <= 0 )
        return 0;

    float P[3];
    for(int i=0; i<N; i++)
    {
        fread((char*)P, 1, 12, fp);
        fiber[0][i] = P[0];
        fiber[1][i] = P[1];
        fiber[2][i] = P[2];
        fseek(fp,4*ns,SEEK_CUR);
    }
    fseek(fp,4*np,SEEK_CUR);

    return N;
}

// Read a fiber from file .tck
unsigned int read_fiberTCK( FILE* fp, float fiber[3][MAX_FIB_LEN], float* ptrToVOXMM )
{
    int i = 0;
    thread_local float P[3];
    fread((char*)P, 1, 12, fp);
    while( !(isnan(P[0])) && !(isnan(P[1])) &&  !(isnan(P[2])) )
    {
        fiber[0][i] = P[0] * ptrToVOXMM[0] + P[1] * ptrToVOXMM[1] + P[2] * ptrToVOXMM[2]  + ptrToVOXMM[3];
        fiber[1][i] = P[0] * ptrToVOXMM[4] + P[1] * ptrToVOXMM[5] + P[2] * ptrToVOXMM[6]  + ptrToVOXMM[7];
        fiber[2][i] = P[0] * ptrToVOXMM[8] + P[1] * ptrToVOXMM[9] + P[2] * ptrToVOXMM[10] + ptrToVOXMM[11];
        i++;
        fread((char*)P, 1, 12, fp);
    }

    return i;
}
