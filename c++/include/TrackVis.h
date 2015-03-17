#ifndef __TRACKVIS_H__
#define __TRACKVIS_H__

#include <COLOR_ui.h>
#include <blitz/array.h>


#define		TRACKVIS_SAVE_ALL		0
#define		TRACKVIS_SAVE_HALF		1
#define		TRACKVIS_SAVE_UNIQUE	2

#define		TRACKVIS_VOXEL_OFFSET	0



// Structure to hold metadata of a TrackVis file
// ---------------------------------------------
struct TrackVis_header
{
    char                id_string[6];
    short int           dim[3];
    float               voxel_size[3];
    float               origin[3];
    short int           n_scalars;
    char                scalar_name[10][20];
    short int           n_properties;
    char                property_name[10][20];
    char                reserved[508];
    char                voxel_order[4];
    char                pad2[4];
    float               image_orientation_patient[6];
    char                pad1[2];
    unsigned char       invert_x;
    unsigned char       invert_y;
    unsigned char       invert_z;
    unsigned char       swap_xy;
    unsigned char       swap_yz;
    unsigned char       swap_zx;
    int                 n_count;
    int                 version;
    int                 hdr_size;
};



// Class to handle TrackVis files.
// -------------------------------
class TrackVis
{
    private:
        string				filename;
        FILE* 				fp;
        int 				maxSteps;  // [TODO] should be related to the variable defined for fiber-tracking

    public:
        TrackVis_header		hdr;

        short 	create( string filename, short int* dim, float* pixdim );
        short 	open( string filename );
        short 	read( blitz::Array<float,2>* fiber, blitz::Array<float,2>* scalars=NULL, blitz::Array<float,2>* properties=NULL );
        short	append( blitz::Array<float,2>* fiber, int numPoints, short saveMethod=TRACKVIS_SAVE_UNIQUE );
        void	writeHdr();
        void	updateTotal( int totFibers );
        void	close();
        FILE*   getFilePtr();

        void 	mm2vox( blitz::Array<float,1> P_mm, blitz::Array<int,1>* P );

        TrackVis();
        ~TrackVis();
};


TrackVis::TrackVis()  { filename = ""; fp = NULL; maxSteps = 20000; }
TrackVis::~TrackVis() { if (fp) fclose( fp ); }


// Create a TrackVis file and store standard metadata. The file is ready to append fibers.
// ---------------------------------------------------------------------------------------
short TrackVis::create( string filename, short int* dim, float* pixdim )
{
    // prepare the header
    for(int i=0; i<3 ;i++)
    {
        if ( dim[i]<=0 || pixdim[i]<=0 ) return 0;
        hdr.dim[i] 			= dim[i];
        hdr.voxel_size[i] 	= pixdim[i];
        hdr.origin[i] 		= 0;
    }
    hdr.n_scalars = 0;
    hdr.n_properties = 0;
    sprintf(hdr.voxel_order,"LPS");
    sprintf(hdr.pad2,"RAS");
    hdr.image_orientation_patient[0] = 1.0;
    hdr.image_orientation_patient[1] = 0.0;
    hdr.image_orientation_patient[2] = 0.0;
    hdr.image_orientation_patient[3] = 0.0;
    hdr.image_orientation_patient[4] = 1.0;
    hdr.image_orientation_patient[5] = 0.0;
    hdr.pad1[0] = 0;
    hdr.pad1[1] = 0;
    hdr.invert_x = 0;
    hdr.invert_y = 0;
    hdr.invert_z = 0;
    hdr.swap_xy = 0;
    hdr.swap_yz = 0;
    hdr.swap_zx = 0;
    hdr.n_count = 0;
    hdr.version = 1;
    hdr.hdr_size = 1000;

    // write the header to the file
    fp = fopen(filename.c_str(),"w+b");
    if (fp == NULL) { printf("\n\n[ERROR] Unable to create file '%s'\n\n",filename.c_str()); return 0; }
    sprintf(hdr.id_string,"TRACK");
    fwrite((char*)&hdr, 1, 1000, fp);

    this->filename = filename;

    return 1;
}



// Open an existing TrackVis file and read metadata information.
// The file pointer is positiond at the beginning of fibers data
// -------------------------------------------------------------
short TrackVis::open( string filename )
{
    fp = fopen(filename.c_str(),"r+b");
    if (fp == NULL) { printf("\n\n[ERROR] Unable to open file '%s'\n\n",filename.c_str()); return 0; }
    this->filename = filename;

    return fread((char*)(&hdr), 1, 1000, fp);
}



// Append a fiber to the file
// --------------------------
short TrackVis::append( blitz::Array<float,2>* fiber, int numPoints, short saveMethod )
{
    unsigned int 	numSaved, pos = 0;
    float 	tmp[3*maxSteps];

    if ( numPoints > maxSteps )
    {
        cerr <<COLOR(30,41,0)<<"[ERROR]"<<COLOR(31,48,0)<<" Trying to write a fiber too long!"<<COLOR_reset<<"\n\n";
        return 0;
    }


    if ( saveMethod == TRACKVIS_SAVE_HALF )
    {
        // Save only 1 POINT OUT OF 2 (in reversed order), but include always the endpoints
        numSaved = ceil((float)(numPoints-1)/2)+1;
        for(int i=numPoints-1; i>0 ;i-=2)
        {
            tmp[pos++] = ( (*fiber)(0,i)+TRACKVIS_VOXEL_OFFSET );
            tmp[pos++] = ( (*fiber)(1,i)+TRACKVIS_VOXEL_OFFSET );
            tmp[pos++] = ( (*fiber)(2,i)+TRACKVIS_VOXEL_OFFSET );
        }
        tmp[pos++] = ( (*fiber)(0,0)+TRACKVIS_VOXEL_OFFSET );
        tmp[pos++] = ( (*fiber)(1,0)+TRACKVIS_VOXEL_OFFSET );
        tmp[pos++] = ( (*fiber)(2,0)+TRACKVIS_VOXEL_OFFSET );
    }
    else if ( saveMethod == TRACKVIS_SAVE_UNIQUE )
    {
        // Save UNIQUE points (discard consecutive points inside the same voxel)
        numSaved = 0;
        int oldX = 0, oldY = 0, oldZ = 0;
        int    X = 0,    Y = 0,    Z = 0;
         for(int i=0; i<numPoints ;i++)
        {
            X = round( (*fiber)(0,i) );
            Y = round( (*fiber)(1,i) );
            Z = round( (*fiber)(2,i) );
            if ( pos==0 || X!=oldX || Y!=oldY || Z!=oldZ )
            {
                tmp[pos++] = ( (*fiber)(0,i)+TRACKVIS_VOXEL_OFFSET );
                tmp[pos++] = ( (*fiber)(1,i)+TRACKVIS_VOXEL_OFFSET );
                tmp[pos++] = ( (*fiber)(2,i)+TRACKVIS_VOXEL_OFFSET );
                oldX = X; oldY = Y; oldZ = Z;
                numSaved++;
            }
        }
    }
    else
    {
        // Save ALL points
        numSaved = numPoints;
         for(unsigned int i=0; i<numSaved ;i++)
        {
            tmp[pos++] = ( (*fiber)(0,(int)i)+TRACKVIS_VOXEL_OFFSET );
            tmp[pos++] = ( (*fiber)(1,(int)i)+TRACKVIS_VOXEL_OFFSET );
            tmp[pos++] = ( (*fiber)(2,(int)i)+TRACKVIS_VOXEL_OFFSET );
        }
    }

    // write the coordinates to the file
    if ( fwrite((char*)&numSaved, 1, 4, fp) != 4 )
    {
        COLOR_error( " Problems saving the fiber!" );
        return 1;
    }
    if ( fwrite((char*)&tmp, 1, 4*pos, fp) != 4*pos )
    {
        COLOR_error( " Problems saving the fiber!" );
        return 1;
    }

    return 0;
}



// Read one fiber from the file
// ----------------------------
short TrackVis::read( blitz::Array<float,2>* fiber, blitz::Array<float,2>* scalars, blitz::Array<float,2>* properties )
{
    int numPoints;
    fread((char*)&numPoints, 1, 4, fp);

    if ( numPoints >= maxSteps || numPoints <= 0 )
    {
        cerr <<COLOR(30,41,0)<<"[ERROR]"<<COLOR(31,48,0)<<" Trying to read a fiber with "<< numPoints <<" points!\n\n"<<COLOR_reset<<"\n\n";
        return -1;
    }

    fiber->resize( 3, numPoints );
    if ( scalars!=NULL )
        scalars->resize( numPoints, hdr.n_scalars );
    if ( properties!=NULL )
        properties->resize( 1, hdr.n_properties );


    float tmp[3];
     for(int i=0; i<numPoints; i++)
     {
        fread((char*)tmp, 1, 12, fp);
         (*fiber)(0,i) = tmp[0];
        (*fiber)(1,i) = tmp[1];
        (*fiber)(2,i) = tmp[2];
         for(int j=0; j<hdr.n_scalars; j++)
         {
             fread((char*)tmp, 1, 4, fp);
             if ( scalars!=NULL )
                 (*scalars)(i,j) = tmp[0];
         }
     }
     for(int j=0; j<hdr.n_properties; j++)
     {
         fread((char*)tmp, 1, 4, fp);
         if ( properties!=NULL )
            (*properties)(j) = tmp[0];
     }

     return numPoints;
}



// Update the field in the header to the new FIBER TOTAL.
// ------------------------------------------------------
void TrackVis::updateTotal( int totFibers )
{
    fseek(fp, 1000-12, SEEK_SET);
    fwrite((char*)&totFibers, 1, 4, fp);
}


void TrackVis::writeHdr()
{
    fseek(fp, 0, SEEK_SET);
    fwrite((char*)&hdr, 1, 1000, fp);
}


// Close the TrackVis file, but keep the metadata in the header.
// -------------------------------------------------------------
void TrackVis::close()
{
    fclose(fp);
    fp = NULL;
}



// Convert a COORDINATE from MILLIMITER to VOXEL.
// ----------------------------------------------
void TrackVis::mm2vox( blitz::Array<float,1> P_mm, blitz::Array<int,1>* P )
{
    (*P)(0) = round( P_mm(0) / hdr.voxel_size[0] - TRACKVIS_VOXEL_OFFSET );
    (*P)(1) = round( P_mm(1) / hdr.voxel_size[1] - TRACKVIS_VOXEL_OFFSET );
    (*P)(2) = round( P_mm(2) / hdr.voxel_size[2] - TRACKVIS_VOXEL_OFFSET );
    if ( (*P)(0)<0 ) (*P)(0) = 0;
    if ( (*P)(1)<0 ) (*P)(1) = 0;
    if ( (*P)(2)<0 ) (*P)(2) = 0;
    if( (*P)(0) > hdr.dim[0]-1 ) (*P)(0) = hdr.dim[0]-1;
    if( (*P)(1) > hdr.dim[1]-1 ) (*P)(1) = hdr.dim[1]-1;
    if( (*P)(2) > hdr.dim[2]-1 ) (*P)(2) = hdr.dim[2]-1;
}


// Return the file pointer
// -----------------------
FILE* TrackVis::getFilePtr()
{
    return fp;
}

#endif
