/*
 * Create a dictionary structure to be used with COMMIT from an input .trk file.
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

#include "tclap/CmdLine.h"
#include "COLOR_ui.h"
#include "ProgressBar.h"
#include "NIFTI.h"
#include "VECTOR.h"
#include "TrackVis.h"
using namespace std;
using namespace blitz;

TrackVis       TRK_input;
NIFTI          niiWM;
bool           hasWM = false;
NIFTI          niiPEAKS;
NIFTI          niiMASK;
VECTOR<int>    dim;
VECTOR<float>  pixdim;


// to store the segments of one fiber
#include <map>

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

map<segKey,float> FiberSegments;
map<segKey,float>::iterator it;
unsigned int    nTrueDirs;
unsigned int    nPointsToSkip;
float           fiberShift, fiberShiftXmm, fiberShiftYmm, fiberShiftZmm;
bool			computeIntersections;
float			ECix, ECiy, ECiz;		// invert axes for EC segments
float			vf_THR;					// threshold to consider the peaks

bool rayBoxIntersection( VECTOR<double>& origin, VECTOR<double>& direction, VECTOR<double>& vmin, VECTOR<double>& vmax, double & t);
void fiberForwardModel( blitz::Array<float,2> & fiber, unsigned short pts );
void segmentForwardModel( const VECTOR<double>& P1, const VECTOR<double>& P2 );



/********************************************************************************************************************/
/*                                                       MAIN                                                       */
/********************************************************************************************************************/
int main(int argc, char** argv)
{
    TCLAP::CmdLine cmd("Takes as input any tractogram in TrackVis format (.trk), and optionally a peaks file in Nifti format (.nii), and creates the linear operator (i.e. dictionary) used by COMMIT. The fiber tracts will generate the intra-cellular (IC) contributions, whereas the extra-cellular (EC) ones will be taken from the peaks file, if present.", ' ', "1.0");

    TCLAP::SwitchArg 			argINTER(   "c", "intersections" , "Compute fiber-lattice intersections", cmd, false );
    TCLAP::ValueArg<float> 		argSHIFT(   "s", "shift"         , "Offset to be applied to coordinates (voxel units)", false, 0, "offset", cmd );
    TCLAP::ValueArg<int> 		argSKIP(    "n", "skip"          , "Points to skip at begin/end of fiber", false, 0, "number", cmd );

    TCLAP::ValueArg<float> 		argVF(      "t", "thr"           , "Threshold for the EC peaks [0..1]", false, 0.1, "float", cmd );
    TCLAP::SwitchArg 			argINVERTZ( "z", "iz"            , "Invert z-axis in EC segments", cmd, false );
    TCLAP::SwitchArg 			argINVERTY( "y", "iy"            , "Invert y-axis in EC segments", cmd, false );
    TCLAP::SwitchArg 			argINVERTX( "x", "ix"            , "Invert x-axis in EC segments", cmd, false );


    TCLAP::ValueArg<string> 	argOUT(     "o", "out"           , "Output path", true, "", "path", cmd );
    TCLAP::ValueArg<string> 	argWM(      "w", "wm"            , "Binary mask to filter segments [.nii]", false, "", "path/filename", cmd );
    TCLAP::ValueArg<string> 	argPEAKS(   "p", "peaks"         , "Diffusion peaks for EC contributions [.nii]", false, "", "path/filename", cmd );
    TCLAP::ValueArg<string> 	argTRK(     "i", "trk"           , "Fiber tracts for IC contributions [.trk]", true, "", "path/filename", cmd );

    try	{ cmd.parse( argc, argv ); }
    catch (TCLAP::ArgException &e) { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }

    string TRK_filename( argTRK.getValue() );
    string PEAKS_filename( argPEAKS.getValue() ), WM_filename( argWM.getValue() ), OUTPUT_path( argOUT.getValue() ), filename;
    fiberShift           = argSHIFT.getValue();
    nPointsToSkip        = argSKIP.getValue();
    computeIntersections = argINTER.getValue();
    ECix = argINVERTX.getValue()==true ? -1 : 1;
    ECiy = argINVERTY.getValue()==true ? -1 : 1;
    ECiz = argINVERTZ.getValue()==true ? -1 : 1;
    vf_THR = argVF.getValue();

    if ( vf_THR<0 || vf_THR>1 )
    {
        COLOR_error("'t' parameter must be in the range [0..1]");
        return EXIT_FAILURE;
    }
    if ( nPointsToSkip<0 )
    {
        COLOR_error("'n' parameter must be non-negative");
        return EXIT_FAILURE;
    }


    /*=========================*/
    /*          SETUP          */
    /*=========================*/
    time_t tStart = time(0);

    printf("\n-> SETTINGS:\n");
    printf( "\t- Segment position = %s\n", computeIntersections ? "COMPUTE INTERSECTIONS" : "CENTROID" );
    printf( "\t- Fiber shift      = %.3f (voxel-size units)\n", fiberShift );
    printf( "\t- Points to skip   = %d\n", nPointsToSkip );


    // open INPUT .TRK file
    // ====================
    printf("\n-> %s:\n", TRK_filename.c_str());
    TRK_input = TrackVis();
    if ( !TRK_input.open( TRK_filename.c_str() ) )
    {
        COLOR_error( "Unable to open input .trk file" );
        return EXIT_FAILURE;
    }

    dim.x 		= TRK_input.hdr.dim[0];
    dim.y 		= TRK_input.hdr.dim[1];
    dim.z 		= TRK_input.hdr.dim[2];
    pixdim.x 	= TRK_input.hdr.voxel_size[0];
    pixdim.y 	= TRK_input.hdr.voxel_size[1];
    pixdim.z 	= TRK_input.hdr.voxel_size[2];
    printf( "\t- %d x %d x %d\n",       dim.x,    dim.y,    dim.z );
    printf( "\t- %.4f x %.4f x %.4f\n", pixdim.x, pixdim.y, pixdim.z );
    printf( "\t- %d fibers\n",     TRK_input.hdr.n_count );
    printf( "\t- %d scalars\n",    TRK_input.hdr.n_scalars );
    printf( "\t- %d properties\n", TRK_input.hdr.n_properties );

    fiberShiftXmm = fiberShift * pixdim.x;
    fiberShiftYmm = fiberShift * pixdim.y;
    fiberShiftZmm = fiberShift * pixdim.z;

    // open WHITE MATTER mask (used to discard segments not in it)
    // ===========================================================
    if ( !WM_filename.empty() )
    {
        printf("\n-> %s:\n", WM_filename.c_str());
        niiWM.open( WM_filename.c_str(), true );
        if ( !niiWM.isValid() ) {
            COLOR_error( "Unable to open file!" );
            return EXIT_FAILURE;
        }
        printf("\t- %d x %d x %d\n",          niiWM.hdr->dim[1],    niiWM.hdr->dim[2],    niiWM.hdr->dim[3]);
        printf("\t- %.4f x %.4f x %.4f\n",    niiWM.hdr->pixdim[1], niiWM.hdr->pixdim[2], niiWM.hdr->pixdim[3]);

        if ( dim.x!=niiWM.hdr->dim[1] || dim.y!=niiWM.hdr->dim[2] || dim.z!=niiWM.hdr->dim[3] ||
             abs(pixdim.x-niiWM.hdr->pixdim[1])>1e-3 || abs(pixdim.y-niiWM.hdr->pixdim[2])>1e-3 || abs(pixdim.z-niiWM.hdr->pixdim[3])>1e-3 )
        {
            COLOR_error( "WM mask and TRK file do not have the same geometry" );
            return EXIT_FAILURE;
        }
        hasWM = true;
    }
    else
    {
        printf("\n-> No mask specified to filter IC compartments\n");
    }


    // open dataset with main DIFFUSION PEAKS
    // ======================================
    if ( !PEAKS_filename.empty() )
    {
        printf("\n-> %s:\n", PEAKS_filename.c_str());
        niiPEAKS.open( PEAKS_filename.c_str(), true );
        if ( !niiPEAKS.isValid() ) {
            COLOR_error( "Unable to open file!" );
            return EXIT_FAILURE;
        }
        printf("\t- %d x %d x %d x %d\n", niiPEAKS.hdr->dim[1], niiPEAKS.hdr->dim[2], niiPEAKS.hdr->dim[3], niiPEAKS.hdr->dim[4]);
        printf("\t- %.4f x %.4f x %.4f\n", niiPEAKS.hdr->pixdim[1], niiPEAKS.hdr->pixdim[2], niiPEAKS.hdr->pixdim[3]);
        printf("\t- Ignoring peaks < %.2f * MaxPeak\n",vf_THR);
        if ( niiPEAKS.hdr->dim[4] % 3 )
        {
            COLOR_error( "DIR dataset must have 3*k volumes!" );
            return EXIT_FAILURE;
        }
        nTrueDirs = niiPEAKS.hdr->dim[4]/3;
    }
    else
    {
        printf("\n-> No dataset specified for EC compartments\n");
    }

    // mask to keep track of which voxels are traversed
    {
        int   dim_vec[3]    = { dim.x,    dim.y,    dim.z };
        float pixdim_vec[3] = { pixdim.x, pixdim.y, pixdim.z };
        niiMASK.make( 3, dim_vec, pixdim_vec, DT_FLOAT32 );
        (*niiMASK.img) = 0;
    }

    niiMASK.copyHeader( niiPEAKS.hdr );
    niiMASK.hdr->dim[0] 	= 3;
    niiMASK.hdr->dim[4]     = 1;
    niiMASK.hdr->datatype 	= DT_FLOAT32;
    niiMASK.hdr->nbyper 	= 4;
    niiMASK.hdr->cal_min	= 0;
    niiMASK.hdr->cal_max	= 1;
    niiMASK.hdr->xyz_units	= 10;
    nifti_update_dims_from_array(niiMASK.hdr);



    /*=========================*/
    /*     IC compartments     */
    /*=========================*/
    printf( "\n-> Exporting IC compartments from TRK file...\n" );

    filename = OUTPUT_path+"/dictionary_IC_trkLen.dict";   FILE* pDict_IC_trkLen = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_IC_f.dict";        FILE* pDict_IC_f      = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_IC_vx.dict";       FILE* pDict_IC_vx     = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_IC_vy.dict";       FILE* pDict_IC_vy     = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_IC_vz.dict";       FILE* pDict_IC_vz     = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_IC_ox.dict";       FILE* pDict_IC_ox     = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_IC_oy.dict";       FILE* pDict_IC_oy     = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_IC_len.dict";      FILE* pDict_IC_len    = fopen( filename.c_str(),   "w" );
    if ( !pDict_IC_trkLen )
    {
        COLOR_error( "Unable to create output files" );
        return EXIT_FAILURE;
    }

    blitz::Array<float,2>   fiber(3,10000);
    unsigned int            pts, totICSegments = 0, totFibers = 0;
    float                   fiberLen;
    VECTOR<double>          P;

    ProgressBar PROGRESS( TRK_input.hdr.n_count, 3 );
    for(int f=0; f<TRK_input.hdr.n_count ;f++)
    {
        PROGRESS.inc();

        // read a fiber
        pts = TRK_input.read( &fiber );

        // compute length (in mm) to discard fibers too short
        fiberLen = 0;
        for(int i=0; i<(int)pts-1 ;i++)
        {
            P.x   = fiber(0,i+1) - fiber(0,i);
            P.y   = fiber(1,i+1) - fiber(1,i);
            P.z   = fiber(2,i+1) - fiber(2,i);
            fiberLen += P.norm();
        }
        // if ( fiberLen < minFiberLen )
            // continue;

        // process each segment separately and add them to the structure "FiberSegments"
        fiberForwardModel( fiber, pts );

        // store data on file
        if ( FiberSegments.size() > 0 )
        {
            fiberLen = 0;
            for (it=FiberSegments.begin(); it!=FiberSegments.end(); it++)
                fiberLen += it->second;
            for (it=FiberSegments.begin(); it!=FiberSegments.end(); it++)
            {
                fwrite( &totFibers,      4, 1, pDict_IC_f );
                fwrite( &(it->first.x),  1, 1, pDict_IC_vx );
                fwrite( &(it->first.y),  1, 1, pDict_IC_vy );
                fwrite( &(it->first.z),  1, 1, pDict_IC_vz );
                fwrite( &(it->first.ox), 1, 1, pDict_IC_ox );
                fwrite( &(it->first.oy), 1, 1, pDict_IC_oy );
                fwrite( &(it->second),   4, 1, pDict_IC_len );
                (*niiMASK.img)( (int)it->first.x, (int)it->first.y, (int)it->first.z ) += it->second;
            }
            totICSegments += FiberSegments.size();
            fwrite( &fiberLen,  1, 4, pDict_IC_trkLen );
            totFibers++;
        }
    }
    PROGRESS.close();

    TRK_input.close();
    fclose( pDict_IC_trkLen );
    fclose( pDict_IC_f );
    fclose( pDict_IC_vx );
    fclose( pDict_IC_vy );
    fclose( pDict_IC_vz );
    fclose( pDict_IC_ox );
    fclose( pDict_IC_oy );
    fclose( pDict_IC_len );

    printf("   [ %d fibers, %d segments, %d voxels ]\n", totFibers, totICSegments, count((*niiMASK.img)>0) );



    /*=========================*/
    /*     EC compartments     */
    /*=========================*/
    printf( "\n-> Exporting EC compartments from PEAKS file...\n" );

    filename = OUTPUT_path+"/dictionary_EC_vx.dict";       FILE* pDict_EC_vx  = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_EC_vy.dict";       FILE* pDict_EC_vy  = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_EC_vz.dict";       FILE* pDict_EC_vz  = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_EC_ox.dict";       FILE* pDict_EC_ox  = fopen( filename.c_str(),   "w" );
    filename = OUTPUT_path+"/dictionary_EC_oy.dict";       FILE* pDict_EC_oy  = fopen( filename.c_str(),   "w" );

    unsigned int totECSegments = 0, totECVoxels = 0;

    if ( !PEAKS_filename.empty() )
    {
        VECTOR<double>      dir;
        double              longitude, colatitude;
        segKey              ec_seg;
        int atLeastOne;
        float peakMax;
        float norms[ nTrueDirs ];

        if ( ECix == -1 )
            printf( "\t- x-axis will be flipped\n" );
        if ( ECiy == -1 )
            printf( "\t- y-axis will be flipped\n" );
        if ( ECiz == -1 )
            printf( "\t- z-axis will be flipped\n" );

        PROGRESS.reset( dim.z );
        for(int iz=0; iz<dim.z ;iz++)
        {
            PROGRESS.inc();
            for(int iy=0; iy<dim.y ;iy++)
            for(int ix=0; ix<dim.x ;ix++)
            {
                if ( (*niiMASK.img)(ix,iy,iz) == 0 ) continue; // not in MASK previously computed from IC segments

                peakMax = -1;
                for(int id=0; id<nTrueDirs ;id++)
                {
                    dir.x = (*niiPEAKS.img)(ix,iy,iz,id*3+0);
                    dir.y = (*niiPEAKS.img)(ix,iy,iz,id*3+1);
                    dir.z = (*niiPEAKS.img)(ix,iy,iz,id*3+2);
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
                    for(int id=0; id<nTrueDirs ;id++)
                    {
                        if ( norms[id]==0 || norms[id] < vf_THR*peakMax ) continue; // peak too small, don't consider it

                        // store this orientation
                        dir.x = ECix * (*niiPEAKS.img)(ix,iy,iz,id*3+0);
                        dir.y = ECiy * (*niiPEAKS.img)(ix,iy,iz,id*3+1);
                        dir.z = ECiz * (*niiPEAKS.img)(ix,iy,iz,id*3+2);
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
                        if ( ec_seg.ox<0 || ec_seg.ox>180 || ec_seg.oy < 0 || ec_seg.oy>180 )
                        {
                            COLOR_error( "segment orientation (ox,oy) out of bound" );
                            return EXIT_FAILURE;
                        }

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

    printf("   [ %d segments, %d voxels ]\n", totECSegments, totECVoxels );



    /*=========================*/
    /*     Save voxel MASK     */
    /*=========================*/
    printf( "\n-> Saving voxel mask...\n" );

    niiMASK.save( OUTPUT_path + "/dictionary_mask.nii" );

    printf( "   [ OK ]\n" );


    printf("\n[ Total time: %.0f seconds]\n\n", difftime(time(0),tStart)  );
    return EXIT_SUCCESS;
}



/********************************************************************************************************************/
/*                                                 fiberForwardModel                                                */
/********************************************************************************************************************/
void fiberForwardModel( blitz::Array<float,2> & fiber, unsigned short pts )
{
    static VECTOR<double>        S1, S2, S1m, S2m, P;
    static VECTOR<double>        vox, vmin, vmax, dir;
    static double                len, t;
    static int                   i, j, k;

    FiberSegments.clear();
    for(i=nPointsToSkip; i<pts-1-nPointsToSkip ;i++)
    {
        // original segment to be processed
        S1.Set( fiber(0,i)   + fiberShiftXmm, fiber(1,i)   + fiberShiftYmm, fiber(2,i)   + fiberShiftZmm );
        S2.Set( fiber(0,i+1) + fiberShiftXmm, fiber(1,i+1) + fiberShiftYmm, fiber(2,i+1) + fiberShiftZmm );

        // get a normal to the vector to move
        dir.x = S2.x-S1.x;
        dir.y = S2.y-S1.y;
        dir.z = S2.z-S1.z;
        dir.Normalize();
        if ( computeIntersections==false )
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
void segmentForwardModel( const VECTOR<double>& P1, const VECTOR<double>& P2 )
{
    static VECTOR<double>    vox;
    static VECTOR<double>    dir, dirTrue;
    static double            longitude, colatitude, len;
    static                   segKey key;

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
    if ( hasWM && (*niiWM.img)((int)vox.x,(int)vox.y,(int)vox.z)==0 )
        return;

    // add the segment to the data structure
    longitude  = atan2(dir.y, dir.x);
    colatitude = atan2( sqrt(dir.x*dir.x + dir.y*dir.y), dir.z );

    key.set( vox.x, vox.y, vox.z, (int)round(colatitude/M_PI*180.0), (int)round(longitude/M_PI*180.0) );
    FiberSegments[key] += len;
    if ( key.ox<0 || key.ox>180 || key.oy < 0 || key.oy>180 )
    {
        COLOR_error( "ox or oy out of LUT bound [0,180]" );
        return;
    }
}



/********************************************************************************************************************/
/*                                                rayBoxIntersection                                                */
/********************************************************************************************************************/
bool rayBoxIntersection( VECTOR<double>& origin, VECTOR<double>& direction, VECTOR<double>& vmin, VECTOR<double>& vmax, double & t)
{
    static double tmin, tmax, tymin, tymax, tzmin, tzmax;
    static VECTOR<double> invrd;

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
