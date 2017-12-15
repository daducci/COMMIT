#include <stdio.h>
#include <cstdio>
#include <string>
#include <map>
#include <vector>
#include "Vector.h"
#include "ProgressBar.h"

#define MAX_FIB_LEN 10000


// CLASS to store the segments of one fiber
class segKey
{
    public:
    unsigned short x, y, z;
    unsigned char ox, oy;
    segKey(){}

    void set(unsigned short _x, unsigned short _y, unsigned short _z, unsigned char _ox, unsigned char _oy)
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
std::map<segKey,float> FiberSegments;

Vector<int>     dim;
Vector<float>   pixdim;
float*          ptrMASK;
unsigned int    nPointsToSkip;
float           fiberShiftXmm, fiberShiftYmm, fiberShiftZmm;
bool            doIntersect;
float           minSegLen;

std::vector<double> radii;         // radii for the extrusion
std::vector<double> weights;       // damping weight
std::vector<int>    sectors;       // number of duplicates across the extrusion circle
double              radiusSigma;   // modulates the impact of each segment as function of radius


bool rayBoxIntersection( Vector<double>& origin, Vector<double>& direction, Vector<double>& vmin, Vector<double>& vmax, double & t);
void fiberForwardModel( float fiber[3][MAX_FIB_LEN], unsigned int pts, std::vector<int> sectors, std::vector<double> radii, std::vector<double> weight );
void segmentForwardModel( const Vector<double>& P1, const Vector<double>& P2, double w );
unsigned int read_fiber( FILE* fp, float fiber[3][MAX_FIB_LEN], int ns, int np );


// =========================
// Function called by CYTHON
// =========================
int trk2dictionary(
    char* strTRKfilename, int Nx, int Ny, int Nz, float Px, float Py, float Pz, int n_count, int n_scalars, int n_properties,
    float fiber_shiftX, float fiber_shiftY, float fiber_shiftZ, int points_to_skip, float min_seg_len,
    float* ptrPEAKS, int Np, float vf_THR, int ECix, int ECiy, int ECiz,
    float* _ptrMASK, float* ptrTDI, char* path_out, int c, double* ptrAFFINE,
    int nBlurRadii, double blurSigma, double* ptrBlurRadii, int* ptrBlurSamples, double* ptrBlurWeights
)
{
    /*=========================*/
    /*     IC compartments     */
    /*=========================*/
    float          fiber[3][MAX_FIB_LEN];
    float          fiberNorm, fiberLen;
    unsigned int   N, totICSegments = 0, totFibers = 0, v;
    unsigned short o;
    unsigned char  kept;
    Vector<double> P;
    std::string    filename;
    std::string    OUTPUT_path(path_out);
    std::map<segKey,float>::iterator it;

    std::map<segInVoxKey,float> FiberNorm;
    std::map<segInVoxKey,float>::iterator itNorm;
    segInVoxKey         inVoxKey;

    printf( "\t* Exporting IC compartments:\n" );

    FILE* fpTRK = fopen(strTRKfilename,"rb");
    if (fpTRK == NULL) return 0;
    fseek(fpTRK,1000,SEEK_SET);

    // set global variables
    dim.Set( Nx, Ny, Nz );
    pixdim.Set( Px, Py, Pz );
    nPointsToSkip = points_to_skip;
    fiberShiftXmm = fiber_shiftX * pixdim.x; // shift in mm for the coordinates
    fiberShiftYmm = fiber_shiftY * pixdim.y;
    fiberShiftZmm = fiber_shiftZ * pixdim.z;
    ptrMASK       = _ptrMASK;
    doIntersect   = c > 0;
    minSegLen     = min_seg_len;

    radii.clear();
    sectors.clear();
    weights.clear();
    for(int i=0; i<nBlurRadii ;i++)
    {
        radii.push_back( ptrBlurRadii[i] );
        sectors.push_back( ptrBlurSamples[i] );
        weights.push_back( ptrBlurWeights[i] );
    }
    radiusSigma = blurSigma;

    // open files
    filename = OUTPUT_path+"/dictionary_TRK_norm.dict";   FILE* pDict_TRK_norm = fopen(filename.c_str(),"wb");
    if ( !pDict_TRK_norm )
    {
        printf( "\n[trk2dictionary] Unable to create output files" );
        return 0;
    }
    filename = OUTPUT_path+"/dictionary_IC_f.dict";        FILE* pDict_IC_f      = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_v.dict";        FILE* pDict_IC_v      = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_o.dict";        FILE* pDict_IC_o      = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_IC_len.dict";      FILE* pDict_IC_len    = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_TRK_len.dict";     FILE* pDict_TRK_len   = fopen(filename.c_str(),"wb");
    filename = OUTPUT_path+"/dictionary_TRK_kept.dict";    FILE* pDict_TRK_kept  = fopen(filename.c_str(),"wb");

    // iterate over fibers
    ProgressBar PROGRESS( n_count );
    PROGRESS.setPrefix("\t  ");
    for(int f=0; f<n_count ;f++)
    {
        PROGRESS.inc();
        N = read_fiber( fpTRK, fiber, n_scalars, n_properties );
        fiberForwardModel( fiber, N, sectors, radii, weights );

        kept = 0;
        if ( FiberSegments.size() > 0 )
        {
            // add segments to files
            fiberNorm = 0;
            fiberLen = 0;
            for (it=FiberSegments.begin(); it!=FiberSegments.end(); it++)
            {
                // NB: plese note inverted ordering for 'v'
                v = it->first.x + dim.x * ( it->first.y + dim.y * it->first.z );
                o = it->first.oy + 181 * it->first.ox;
                fwrite( &totFibers,      4, 1, pDict_IC_f );
                fwrite( &v,              4, 1, pDict_IC_v );
                fwrite( &o,              2, 1, pDict_IC_o );
                fwrite( &(it->second),   4, 1, pDict_IC_len );
                ptrTDI[ it->first.z + dim.z * ( it->first.y + dim.y * it->first.x ) ] += it->second;
                inVoxKey.set( it->first.x, it->first.y, it->first.z );
                FiberNorm[inVoxKey] += it->second;
                fiberLen += it->second;
            }
            for (itNorm=FiberNorm.begin(); itNorm!=FiberNorm.end(); itNorm++)
            {
                fiberNorm += pow(itNorm->second,2);
            }
            fiberNorm = sqrt(fiberNorm);
            FiberNorm.clear();
            fwrite( &fiberNorm,  1, 4, pDict_TRK_norm ); // actual length considered in optimization
            fwrite( &fiberLen,   1, 4, pDict_TRK_len );
            totICSegments += FiberSegments.size();
            totFibers++;
            kept = 1;
        }
        fwrite( &kept, 1, 1, pDict_TRK_kept );
    }
    PROGRESS.close();

    fclose( fpTRK );
    fclose( pDict_TRK_norm );
    fclose( pDict_IC_f );
    fclose( pDict_IC_v );
    fclose( pDict_IC_o );
    fclose( pDict_IC_len );
    fclose( pDict_TRK_len );
    fclose( pDict_TRK_kept );

    printf("\t  [ %d fibers kept, %d segments in total ]\n", totFibers, totICSegments );


    /*=========================*/
    /*     EC compartments     */
    /*=========================*/
    unsigned int totECSegments = 0, totECVoxels = 0;

    printf( "\t* Exporting EC compartments:\n" );

    filename = OUTPUT_path+"/dictionary_EC_v.dict";        FILE* pDict_EC_v   = fopen( filename.c_str(),   "wb" );
    filename = OUTPUT_path+"/dictionary_EC_o.dict";        FILE* pDict_EC_o   = fopen( filename.c_str(),   "wb" );

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

                        // get the orientation of the current peak
                        ptr = ptrPEAKS + 3*(id + Np * ( iz + dim.z * ( iy + dim.y * ix ) ));

                        // multiply by the affine matrix
                        dir.x = ptr[0] * ptrAFFINE[0] + ptr[1] * ptrAFFINE[1] + ptr[2] * ptrAFFINE[2];
                        dir.y = ptr[0] * ptrAFFINE[3] + ptr[1] * ptrAFFINE[4] + ptr[2] * ptrAFFINE[5];
                        dir.z = ptr[0] * ptrAFFINE[6] + ptr[1] * ptrAFFINE[7] + ptr[2] * ptrAFFINE[8];

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
                        ec_seg.ox = (int)round(colatitude/M_PI*180.0);
                        ec_seg.oy = (int)round(longitude/M_PI*180.0);

                        v = ec_seg.x + dim.x * ( ec_seg.y + dim.y * ec_seg.z );
                        o = ec_seg.oy + 181 * ec_seg.ox;
                        fwrite( &v, 4, 1, pDict_EC_v );
                        fwrite( &o, 2, 1, pDict_EC_o );
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

    fclose( pDict_EC_v );
    fclose( pDict_EC_o );

    printf("\t  [ %d voxels, %d segments ]\n", totECVoxels, totECSegments );

    return 1;
}


/********************************************************************************************************************/
/*                                                 fiberForwardModel                                                */
/********************************************************************************************************************/
void fiberForwardModel( float fiber[3][MAX_FIB_LEN], unsigned int pts, std::vector<int> sectors, std::vector<double> radii, std::vector<double> weights )
{
    static Vector<double> S1, S2, S1m, S2m, P, q, n, qxn, qxqxn;
    static Vector<double> vox, vmin, vmax, dir;
    static double         len, t, alpha, w, R;
    static int            i, j, k;

    FiberSegments.clear();
    for(i=nPointsToSkip; i<pts-1-nPointsToSkip ;i++)
    {
        // original segment to be processed
        S1.Set( fiber[0][i]   + fiberShiftXmm, fiber[1][i]   + fiberShiftYmm, fiber[2][i]   + fiberShiftZmm );
        S2.Set( fiber[0][i+1] + fiberShiftXmm, fiber[1][i+1] + fiberShiftYmm, fiber[2][i+1] + fiberShiftZmm );
        dir.x = S2.x-S1.x;
        dir.y = S2.y-S1.y;
        dir.z = S2.z-S1.z;
        dir.Normalize();

        // get a normal to the vector to move
        n.x = dir.y-dir.z;
        n.y = dir.z-dir.x;
        n.z = dir.x-dir.y;
        n.Normalize();

        /* assign contribution(s) */
        for(k=0; k<(int)radii.size() ;k++)
        {
            if ( weights[k] < 1e-3 )
                continue;

            R = radii[k];

            // quaternion (q.x, q.y, q.z, w) for rotation
            alpha = 2.0*M_PI/sectors[k];
            w = sin(alpha/2.0);
            q.x = dir.x * w;
            q.y = dir.y * w;
            q.z = dir.z * w;
            w = cos(alpha/2.0);


            for(j=0; j<sectors[k] ;j++)
            {
                // rotate the segment's normal
                qxn.x = 2.0 * ( q.y * n.z - q.z * n.y );
                qxn.y = 2.0 * ( q.z * n.x - q.x * n.z );
                qxn.z = 2.0 * ( q.x * n.y - q.y * n.x );
                qxqxn.x = q.y * qxn.z - q.z * qxn.y;
                qxqxn.y = q.z * qxn.x - q.x * qxn.z;
                qxqxn.z = q.x * qxn.y - q.y * qxn.x;
                n.x += w * qxn.x + qxqxn.x;
                n.y += w * qxn.y + qxqxn.y;
                n.z += w * qxn.z + qxqxn.z;
                // n /= np.linalg.norm(n)

                // move the segment
                S1m.x = S1.x + R*n.x;
                S1m.y = S1.y + R*n.y;
                S1m.z = S1.z + R*n.z;
                S2m.x = S2.x + R*n.x;
                S2m.y = S2.y + R*n.y;
                S2m.z = S2.z + R*n.z;

                if ( doIntersect==false )
                    segmentForwardModel( S1m, S2m, weights[k] );
                else
                    while( 1 )
                    {
                        len = sqrt( pow(S2m.x-S1m.x,2) + pow(S2m.y-S1m.y,2) + pow(S2m.z-S1m.z,2) ); // in mm
                        if ( len <= minSegLen )
                            break;

                        // compute AABB of the first point (in mm)
                        vmin.x = floor( (S1m.x + 1e-6*dir.x)/pixdim.x ) * pixdim.x;
                        vmin.y = floor( (S1m.y + 1e-6*dir.y)/pixdim.y ) * pixdim.y;
                        vmin.z = floor( (S1m.z + 1e-6*dir.z)/pixdim.z ) * pixdim.z;
                        vmax.x = vmin.x + pixdim.x;
                        vmax.y = vmin.y + pixdim.y;
                        vmax.z = vmin.z + pixdim.z;

                        if ( rayBoxIntersection( S1m, dir, vmin, vmax, t ) && t>0 && t<len )
                        {
                            // add the portion S1P, and then reiterate
                            P.Set( S1m.x + t*dir.x, S1m.y + t*dir.y, S1m.z + t*dir.z );
                            segmentForwardModel( S1m, P, weights[k] );
                            S1m.Set( P.x, P.y, P.z );
                        }
                        else
                        {
                            // add the segment S1S2 and stop iterating
                            segmentForwardModel( S1m, S2m, weights[k] );
                            break;
                        }
                    }
            }
        }
    }
}


/********************************************************************************************************************/
/*                                                segmentForwardModel                                               */
/********************************************************************************************************************/
void segmentForwardModel( const Vector<double>& P1, const Vector<double>& P2, double w )
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
    if ( len <= minSegLen )
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

    // add the segment to the data structure
    longitude  = atan2(dir.y, dir.x);
    colatitude = atan2( sqrt(dir.x*dir.x + dir.y*dir.y), dir.z );
    key.set( vox.x, vox.y, vox.z, (int)round(colatitude/M_PI*180.0), (int)round(longitude/M_PI*180.0) );
    FiberSegments[key] += w * len;
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
