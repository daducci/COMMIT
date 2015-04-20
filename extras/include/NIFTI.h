#ifndef __NIFTI_H__
#define __NIFTI_H__

#include <blitz/array.h>
using namespace blitz;
#include <nifti1_io.h>


#ifdef INT_2_BYTES
    typedef	char				INT8;
    typedef	unsigned char		UINT8;
    typedef	int					INT16;
    typedef	unsigned int		UINT16;
    typedef	long				INT32;
    typedef	unsigned long		UINT32;
    typedef	double				FLOAT32;
#else
    typedef	char				INT8;
    typedef	unsigned char		UINT8;
    typedef	short				INT16;
    typedef	unsigned short		UINT16;
    typedef	int					INT32;
    typedef	unsigned int		UINT32;
    typedef	float				FLOAT32;
#endif


/* Errorcode used for the field "errorCode" */
#define		NIFTI_ERROR_NOERROR				0
#define		NIFTI_ERROR_WRONGFILETYPE		1
#define		NIFTI_ERROR_DATANOTLOADED		2
#define		NIFTI_ERROR_UNKNOWN				9



/*****************************************************
  *****               NIFTI class                *****
  ****************************************************/
class NIFTI
{
    public:
        nifti_image* 	  hdr;
        Array<FLOAT32,7>* img;

        short			isValid() { return errorCode==NIFTI_ERROR_NOERROR; };
        short			getErrorCode() { return errorCode; };

        short			make( const int ndims, const int* dim, const float* pixdim, const short datatype = DT_FLOAT32 );
        short			open( string filename, bool loadData = true );
        short			load();
        short			unload();
        short			save( string NEW_filename = "" );
        void			copyHeader( const nifti_image* src );

        NIFTI( string filename, bool loadData = true );
        NIFTI();
        ~NIFTI();

    private:
        short			errorCode;
        // short			getDatatypeCode();
};


/* Constructor/destructor */
NIFTI::NIFTI( void )
{
    img			= NULL;
    hdr			= NULL;
    errorCode 	= NIFTI_ERROR_DATANOTLOADED;
};


NIFTI::NIFTI( string filename, bool loadData )
{
    img		= NULL;
    hdr		= NULL;
    this->open( filename, loadData );
}


NIFTI::~NIFTI()
{
     if ( hdr ) nifti_image_unload( hdr );
     if ( img ) img->free();
}


/* OPEN a nifti file (only the header is loaded) */
short NIFTI::open( string filename, bool loadData )
{
    unload();

    try
    {
        // not a NIFTI file
        if ( is_nifti_file(filename.c_str()) < 1 ) { errorCode = NIFTI_ERROR_WRONGFILETYPE; return 0; }

        hdr = nifti_image_read( filename.c_str(), 0 );
        if ( hdr==NULL ) { errorCode = NIFTI_ERROR_DATANOTLOADED; return 0; }
    }
    catch(exception& ex)
    {
        errorCode = NIFTI_ERROR_UNKNOWN;
        return 0;
    }

    // correct the "dim" field if it contains zeros
    for(int i=1; i<8 ;i++)
        if ( hdr->dim[i]==0 ) hdr->dim[i] = 1;

    errorCode = NIFTI_ERROR_NOERROR;
    if ( loadData ) return this->load();
    return 1;
}


/*  MAKE a new dataset  */
short NIFTI::make( const int ndims, const int* dim, const float* pixdim, const short datatype )
{
    if ( ndims<1 || ndims>7 ) return 0;
    if ( datatype!=DT_INT8 && datatype!=DT_UINT8 && datatype!=DT_INT16 && datatype!=DT_UINT16 && datatype!=DT_INT32 && datatype!=DT_UINT32 && datatype!=DT_FLOAT32 ) return 0;

    int   d[8] = {0,1,1,1,1,1,1,1};
    float p[8] = {1,1,1,1,1,1,1,1};
    for(int i=0; i<ndims ;i++)
        { d[i+1] = dim[i]; p[i+1] = pixdim[i]; }
    d[0] = ndims;

    nifti_image_unload( hdr );
    hdr = nifti_make_new_nim(d, datatype, 1);
    for(int i=0; i<ndims ;i++)
        hdr->pixdim[i+1] = p[i+1];
     nifti_update_dims_from_array( hdr );

    try
    {
        if ( img ) img->free();
        img = new Array<FLOAT32,7>(
            shape(hdr->dim[1], hdr->dim[2], hdr->dim[3], hdr->dim[4], hdr->dim[5], hdr->dim[6], hdr->dim[7]),
            ColumnMajorArray<7>()
        );

        // cast data to FLOAT 32 (internal format)
        FLOAT32* outPtr    = (FLOAT32*)( img->data() );
        FLOAT32* outPtrEnd = outPtr + hdr->nvox;
        switch( hdr->datatype )
        {
            case DT_INT8:
            {
                INT8* inPtr = (INT8*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_UINT8:
            {
                UINT8* inPtr = (UINT8*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_INT16:
            {
                INT16* inPtr = (INT16*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_UINT16:
            {
                UINT16* inPtr = (UINT16*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_INT32:
            {
                INT32* inPtr = (INT32*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_UINT32:
            {
                UINT32* inPtr = (UINT32*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_FLOAT32:
            {
                FLOAT32* inPtr = (FLOAT32*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            default:
                return 0;
        }
    }
    catch(exception& ex)
    {
        return 0;
    }

     return (img==NULL?0:1);
}



/*  LOAD/UNLOAD data  */
short NIFTI::load()
{
    if ( errorCode>0 ) return 0;
    if ( nifti_image_load(hdr) < 0 ) return 0;

    try
    {
        if ( img ) img->free();
        img = new Array<FLOAT32,7>(
            shape(hdr->dim[1], hdr->dim[2], hdr->dim[3], hdr->dim[4], hdr->dim[5], hdr->dim[6], hdr->dim[7]),
            ColumnMajorArray<7>()
        );

        // cast data to FLOAT 32 (internal format)
        FLOAT32* outPtr    = (FLOAT32*)( img->data() );
        FLOAT32* outPtrEnd = outPtr + hdr->nvox;
        switch( hdr->datatype )
        {
            case DT_INT8:
            {
                INT8* inPtr = (INT8*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_UINT8:
            {
                UINT8* inPtr = (UINT8*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_INT16:
            {
                INT16* inPtr = (INT16*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_UINT16:
            {
                UINT16* inPtr = (UINT16*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_INT32:
            {
                INT32* inPtr = (INT32*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_UINT32:
            {
                UINT32* inPtr = (UINT32*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            case DT_FLOAT32:
            {
                FLOAT32* inPtr = (FLOAT32*)(hdr->data);
                while( outPtr != outPtrEnd )
                    *(outPtr++) = (FLOAT32)(*inPtr++);
                break;
            }
            default:
                return 0;
        }
    }
    catch(exception& ex)
    {
        return 0;
    }

    return (img==NULL?0:1);
}


/* Unload the data but keep the metadata */
short NIFTI::unload( )
{
    nifti_image_unload( hdr );
    if ( img ) img->free();
        return 1;
}


/*  SAVE data  */
short NIFTI::save( string NEW_filename )
{
    if ( !nifti_validfilename( NEW_filename.c_str() ) ) return 0;
    nifti_set_filenames( hdr, NEW_filename.c_str(), 0, hdr->byteorder );

    try
    {
        // cast data from FLOAT 32 (internal format)
        FLOAT32* inPtr    = (FLOAT32*)( img->data() );
        FLOAT32* inPtrEnd = inPtr + hdr->nvox;
        switch( hdr->datatype )
        {
            case DT_INT8:
            {
                INT8* outPtr = (INT8*)(hdr->data);
                while( inPtr != inPtrEnd )
                    *(outPtr++) = (INT8)(*inPtr++);
                break;
            }
            case DT_UINT8:
            {
                UINT8* outPtr = (UINT8*)(hdr->data);
                while( inPtr != inPtrEnd )
                    *(outPtr++) = (UINT8)(*inPtr++);
                break;
            }
            case DT_INT16:
            {
                INT16* outPtr = (INT16*)(hdr->data);
                while( inPtr != inPtrEnd )
                    *(outPtr++) = (INT16)(*inPtr++);
                break;
            }
            case DT_UINT16:
            {
                UINT16* outPtr = (UINT16*)(hdr->data);
                while( inPtr != inPtrEnd )
                    *(outPtr++) = (UINT16)(*inPtr++);
                break;
            }
            case DT_INT32:
            {
                INT32* outPtr = (INT32*)(hdr->data);
                while( inPtr != inPtrEnd )
                    *(outPtr++) = (INT32)(*inPtr++);
                break;
            }
            case DT_UINT32:
            {
                UINT32* outPtr = (UINT32*)(hdr->data);
                while( inPtr != inPtrEnd )
                    *(outPtr++) = (UINT32)(*inPtr++);
                break;
            }
            case DT_FLOAT32:
            {
                FLOAT32* outPtr = (FLOAT32*)(hdr->data);
                while( inPtr != inPtrEnd )
                    *(outPtr++) = (FLOAT32)( *inPtr++);
                break;
            }
            default:
                return 0;
        }
    }
    catch(exception& ex)
    {
        return 0;
    }

    nifti_image_write( hdr );
    return 1;
}


void NIFTI::copyHeader( const nifti_image* src )
{
    if ( !src ) return;

    void* tmp = hdr->data;
    hdr = nifti_copy_nim_info( src );
    hdr->data = tmp;
}

#endif
