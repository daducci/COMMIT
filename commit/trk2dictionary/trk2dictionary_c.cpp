#include <stdio.h>
#include <cstdio>
#include <string>
#include <map>
#include "Vector.h"
#include "ProgressBar.h"

#define MAX_FIB_LEN 10000


// CLASS to store the segments of one fiber
class segKey
{
    public:
    unsigned char x, y, z, ox, oy;
    segKey(){}

    void set(unsigned char _x, unsigned char _y, unsigned char _z, unsigned char _ox, unsigned char _oy)
    {
        x  = _x;
        y  = _y;
        z  = _z;
        ox = _ox;
        oy = _oy;
    }

    bool const operator <(const segKey& o) const
    {
        return oy<o.oy || (oy==o.oy && ox<o.ox) || (oy==o.oy && ox==o.ox && z<o.z) || (oy==o.oy && ox==o.ox && z==o.z && y<o.y) || (oy==o.oy && ox==o.ox && z==o.z && y==o.y && x<o.x);
    }
    bool const operator ==(const segKey& o) const
    {
        return oy<o.oy && ox==o.ox && z==o.z && y==o.y && x==o.x;
    }
};

// global variables (to avoid passing them at each call)
std::map<segKey,float> FiberSegments;

Vector<int>     dim;
Vector<float>   pixdim;
float*          ptrMASK;
unsigned int    nPointsToSkip;
float           fiberShiftXmm, fiberShiftYmm, fiberShiftZmm;
bool            doIntersect;

bool rayBoxIntersection( Vector<double>& origin, Vector<double>& direction, Vector<double>& vmin, Vector<double>& vmax, double & t);
void fiberForwardModel( float fiber[3][MAX_FIB_LEN], unsigned int pts );
void segmentForwardModel( const Vector<double>& P1, const Vector<double>& P2 );
unsigned int read_fiber( FILE* fp, float fiber[3][MAX_FIB_LEN], int ns, int np );


// =========================
// Function called by CYTHON
// =========================
int trk2dictionary(
    char* strTRKfilename, int Nx, int Ny, int Nz, float Px, float Py, float Pz, int n_count, int n_scalars, int n_properties, float fiber_shift, int points_to_skip,
    float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
    float* _ptrMASK, float* ptrTDI, char* path_out, int c
)
{
    /*=========================*/
    /*     IC compartments     */
    /*=========================*/
    float          fiber[3][MAX_FIB_LEN];
    float          fiberLen;
    unsigned int   N, totICSegments = 0, totFibers = 0;
    Vector<double> P;
    std::string    filename;
    std::string    OUTPUT_path(path_out);
    std::map<segKey,float>::iterator it;

    printf( "\t* Exporting IC compartments:\n" );

    FILE* fpTRK = fopen(strTRKfilename,"r+b");
    if (fpTRK == NULL) return 0;
    fseek(fpTRK,1000,SEEK_SET);

    // set global variables
    dim.Set( Nx, Ny, Nz );
    pixdim.Set( Px, Py, Pz );
    nPointsToSkip = points_to_skip;
    fiberShiftXmm = fiber_shift * pixdim.x; // shift in mm for the coordinates
    fiberShiftYmm = fiber_shift * pixdim.y;
    fiberShiftZmm = fiber_shift * pixdim.z;
    ptrMASK       = _ptrMASK;
    doIntersect   = c > 0;

    // open files
    filename = OUTPUT_path+"/dictionary_IC_trkLen.dict";   FILE* pDict_IC_trkLen = fopen(filename.c_str(),"wb");
    if ( !pDict_IC_trkLen )
    {
        printf( "\n[trk2dictionary] Unable to create output files" );
        return 0;
    }
    filename = OUTPUT_path+"/dictionary_IC_f.dict";        FILE* pDict_IC_f      = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_vx.dict";       FILE* pDict_IC_vx     = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_vy.dict";       FILE* pDict_IC_vy     = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_vz.dict";       FILE* pDict_IC_vz     = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_ox.dict";       FILE* pDict_IC_ox     = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_oy.dict";       FILE* pDict_IC_oy     = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_len.dict";      FILE* pDict_IC_len    = fopen(filename.c_str(),"wb");

    // iterate over fibers
    ProgressBar PROGRESS( n_count );
    PROGRESS.setPrefix("\t  ");
    for(int f=0; f<n_count ;f++)
    {
        PROGRESS.inc();
        N = read_fiber( fpTRK, fiber, n_scalars, n_properties );
        fiberForwardModel( fiber, N );

        if ( FiberSegments.size() > 0 )
        {
            // store data to files
            fiberLen = 0;
            for (it=FiberSegments.begin(); it!=FiberSegments.end(); it++)
            {
                fwrite( &totFibers,      4, 1, pDict_IC_f );
                fwrite( &(it->first.x),  1, 1, pDict_IC_vx );
                fwrite( &(it->first.y),  1, 1, pDict_IC_vy );
                fwrite( &(it->first.z),  1, 1, pDict_IC_vz );
                fwrite( &(it->first.ox), 1, 1, pDict_IC_ox );
                fwrite( &(it->first.oy), 1, 1, pDict_IC_oy );
                fwrite( &(it->second),   4, 1, pDict_IC_len );
                ptrTDI[ it->first.z + dim.z * ( it->first.y + dim.y * it->first.x ) ] += it->second;
                fiberLen += it->second;
            }
            fwrite( &fiberLen,  1, 4, pDict_IC_trkLen );
            totICSegments += FiberSegments.size();
            totFibers++;
        }
    }
    PROGRESS.close();

    fclose( fpTRK );
    fclose( pDict_IC_trkLen );
    fclose( pDict_IC_f );
    fclose( pDict_IC_vx );
    fclose( pDict_IC_vy );
    fclose( pDict_IC_vz );
    fclose( pDict_IC_ox );
    fclose( pDict_IC_oy );
    fclose( pDict_IC_len );

    printf("\t  [ %d fibers, %d segments ]\n", totFibers, totICSegments );


    /*=========================*/
    /*     EC compartments     */
    /*=========================*/
    unsigned int totECSegments = 0, totECVoxels = 0;

    printf( "\t* Exporting EC compartments:\n" );

    filename = OUTPUT_path+"/dictionary_EC_vx.dict";       FILE* pDict_EC_vx  = fopen( filename.c_str(),   "wb" );
    filename = OUTPUT_path+"/dictionary_EC_vy.dict";       FILE* pDict_EC_vy  = fopen( filename.c_str(),   "wb" );
    filename = OUTPUT_path+"/dictionary_EC_vz.dict";       FILE* pDict_EC_vz  = fopen( filename.c_str(),   "wb" );
    filename = OUTPUT_path+"/dictionary_EC_ox.dict";       FILE* pDict_EC_ox  = fopen( filename.c_str(),   "wb" );
    filename = OUTPUT_path+"/dictionary_EC_oy.dict";       FILE* pDict_EC_oy  = fopen( filename.c_str(),   "wb" );

    if ( ptrPEAKS != NULL )
    {
        Vector<double> dir;
        double         longitude, colatitude;
        segKey         ec_seg;
        int            ix, iy, iz, id, atLeastOne;
        float          peakMax;
        float          norms[ Np ];
        float          *ptr;

        PROGRESS.reset( dim.z );
        for(iz=0; iz<dim.z ;iz++)
        {
            PROGRESS.inc();
            for(iy=0; iy<dim.y ;iy++)
            for(ix=0; ix<dim.x ;ix++)
            {
                // check if in mask previously computed from IC segments
                if ( ptrTDI[ iz + dim.z * ( iy + dim.y * ix ) ] == 0 ) continue;

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

                        // store this orientation (invert axes if needed)
                        ptr = ptrPEAKS + 3*(id + Np * ( iz + dim.z * ( iy + dim.y * ix ) ));
                        dir.x = ECix * ptr[0];
                        dir.y = ECiy * ptr[1];
                        dir.z = ECiz * ptr[2];
                        if ( dir.y < 0 )
                        {
                            // ensure to be in the right hemisphere (the one where kernels were pre-computed)
                            dir.x = -dir.x;
                            dir.y = -dir.y;
                            dir.z = -dir.z;
                        }
                        colatitude = atan2( sqrt(dir.x*dir.x + dir.y*dir.y), dir.z );
                        longitude  = atan2( dir.y, dir.x );
                        ec_seg.ox = (int)round(colatitude/M_PI*180.0);
                        ec_seg.oy = (int)round(longitude/M_PI*180.0);
                        fwrite( &ec_seg.x,   1, 1, pDict_EC_vx );
                        fwrite( &ec_seg.y,   1, 1, pDict_EC_vy );
                        fwrite( &ec_seg.z,   1, 1, pDict_EC_vz );
                        fwrite( &ec_seg.ox,  1, 1, pDict_EC_ox );
                        fwrite( &ec_seg.oy,  1, 1, pDict_EC_oy );
                        totECSegments++;
                        atLeastOne = 1;
                    }
                    if ( atLeastOne>0 )
                        totECVoxels++;
                }
            }
        }
        PROGRESS.close();
    }

    fclose( pDict_EC_vx );
    fclose( pDict_EC_vy );
    fclose( pDict_EC_vz );
    fclose( pDict_EC_ox );
    fclose( pDict_EC_oy );

    printf("\t  [ %d voxels, %d segments ]\n", totECVoxels, totECSegments );

    return 1;
}


/********************************************************************************************************************/
/*                                                 fiberForwardModel                                                */
/********************************************************************************************************************/
void fiberForwardModel( float fiber[3][MAX_FIB_LEN], unsigned int pts )
{
    static Vector<double> S1, S2, S1m, S2m, P;
    static Vector<double> vox, vmin, vmax, dir;
    static double         len, t;
    static int            i, j, k;

    FiberSegments.clear();
    for(i=nPointsToSkip; i<pts-1-nPointsToSkip ;i++)
    {
        // original segment to be processed
        S1.Set( fiber[0][i]   + fiberShiftXmm, fiber[1][i]   + fiberShiftYmm, fiber[2][i]   + fiberShiftZmm );
        S2.Set( fiber[0][i+1] + fiberShiftXmm, fiber[1][i+1] + fiberShiftYmm, fiber[2][i+1] + fiberShiftZmm );

        // get a normal to the vector to move
        dir.x = S2.x-S1.x;
        dir.y = S2.y-S1.y;
        dir.z = S2.z-S1.z;
        dir.Normalize();
        if ( doIntersect==false )
            segmentForwardModel( S1, S2 );
        else
            while( 1 )
            {
                len = sqrt( pow(S2.x-S1.x,2) + pow(S2.y-S1.y,2) + pow(S2.z-S1.z,2) ); // in mm
                if ( len < 1e-3 )
                    break;

                // compute AABB of the first point (in mm)
                vmin.x = floor( (S1.x + 1e-6*dir.x)/pixdim.x ) * pixdim.x;
                vmin.y = floor( (S1.y + 1e-6*dir.y)/pixdim.y ) * pixdim.y;
                vmin.z = floor( (S1.z + 1e-6*dir.z)/pixdim.z ) * pixdim.z;
                vmax.x = vmin.x + pixdim.x;
                vmax.y = vmin.y + pixdim.y;
                vmax.z = vmin.z + pixdim.z;

                if ( rayBoxIntersection( S1, dir, vmin, vmax, t ) && t>0 && t<len )
                {
                    // add the portion S1P, and then reiterate
                    P.Set( S1.x + t*dir.x, S1.y + t*dir.y, S1.z + t*dir.z );
                    segmentForwardModel( S1, P );
                    S1.Set( P.x, P.y, P.z );
                }
                else
                {
                    // add the segment S1S2 and stop iterating
                    segmentForwardModel( S1, S2 );
                    break;
                }
            }
    }
}


/********************************************************************************************************************/
/*                                                segmentForwardModel                                               */
/********************************************************************************************************************/
void segmentForwardModel( const Vector<double>& P1, const Vector<double>& P2 )
{
    static Vector<int>    vox;
    static Vector<double> dir, dirTrue;
    static double         longitude, colatitude, len;
    static segKey         key;

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

    // length of segment
    len = dir.norm();
    if ( len<1e-4 ) return; // in mm
    dir.Normalize();

    // voxel of the segment is the centroid
    vox.x = floor( 0.5 * (P1.x + P2.x) / pixdim.x );
    vox.y = floor( 0.5 * (P1.y + P2.y) / pixdim.y );
    vox.z = floor( 0.5 * (P1.z + P2.z) / pixdim.z );
    if ( vox.x>=dim.x || vox.x<0 || vox.y>=dim.y || vox.y<0 || vox.z>=dim.z || vox.z<0 )
        return;
    if ( ptrMASK && ptrMASK[ vox.z + dim.z * ( vox.y + dim.y * vox.x ) ]==0 )
        return;

    // add the segment to the data structure
    longitude  = atan2(dir.y, dir.x);
    colatitude = atan2( sqrt(dir.x*dir.x + dir.y*dir.y), dir.z );
    key.set( vox.x, vox.y, vox.z, (int)round(colatitude/M_PI*180.0), (int)round(longitude/M_PI*180.0) );
    FiberSegments[key] += len;
}


/********************************************************************************************************************/
/*                                                rayBoxIntersection                                                */
/********************************************************************************************************************/
bool rayBoxIntersection( Vector<double>& origin, Vector<double>& direction, Vector<double>& vmin, Vector<double>& vmax, double & t)
{
    static double tmin, tmax, tymin, tymax, tzmin, tzmax;
    static Vector<double> invrd;

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


// Read a fiber from file
unsigned int read_fiber( FILE* fp, float fiber[3][MAX_FIB_LEN], int ns, int np )
{
    int N;
    fread((char*)&N, 1, 4, fp);

    if ( N >= MAX_FIB_LEN || N <= 0 )
        return 0;

    float tmp[3];
    for(int i=0; i<N; i++)
    {
        fread((char*)tmp, 1, 12, fp);
        fiber[0][i] = tmp[0];
        fiber[1][i] = tmp[1];
        fiber[2][i] = tmp[2];
        fseek(fp,4*ns,SEEK_CUR);
    }
    fseek(fp,4*np,SEEK_CUR);

    return N;
}
