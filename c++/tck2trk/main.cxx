/*
 * Title
 *
 * @param
 *
 */
#include <iostream>
#include <fstream>
#include <string>

#include "tclap/CmdLine.h"
#include "COLOR_ui.h"
#include "NIFTI.h"
#include "TrackVis.h"
#include "VECTOR.h"

#include <map>

using namespace std;
using namespace blitz;


// ==============================
// Trim a string (left and right)
// ==============================
void trim( string & S )
{
	int k1, k2;
	for(k1=0; k1<S.size() ;k1++)
		if ( S.at(k1)!=' ' ) break;
	for(k2=S.size()-1; k2>=0 ;k2--)
		if ( S.at(k2)!=' ' ) break;
	if (k1<=k2)
		S = S.substr(k1,k2-k1+1);
	else
		S = "";
}

int main(int argc, char** argv)
{
	map<string,string> HDR;


/*--------------------------*/
/*  Check INPUT parameters  */
/*--------------------------*/

	TCLAP::CmdLine cmd("alessandro.daducci@epfl.ch", ' ', "0.2");

	TCLAP::ValueArg<string> 	argOUT(  "o", "out"  , "Output fibers [.trk]", true, "", "path/filename", cmd);
	TCLAP::ValueArg<string> 	argREF(  "r", "ref"  , "White-matter binary mask [.nii]", true, "", "path/filename", cmd);
	TCLAP::ValueArg<string> 	argIN(   "i", "in"   , "Input fibers [.tck]", true, "", "path/filename", cmd);

	try	{ cmd.parse( argc, argv ); }
	catch (TCLAP::ArgException &e) { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }

	string IN_filename(  argIN.getValue()  );
	string OUT_filename( argOUT.getValue() );
	string REF_filename( argREF.getValue() );


/*----------------------*/
/*  Opening 'TCK' file  */
/*----------------------*/
	cout <<COLOR(37,37,1)<< "\n-> Opening MRTRIX input file...\n" <<COLOR_reset;

	FILE* fp = fopen(IN_filename.c_str(),"r+b");
	if (fp == NULL)
	{
		COLOR_error("Unable to open the input file");
		return EXIT_FAILURE;
	}

	char buffer[1000];
	fscanf( fp, "%100[^\n]\n", buffer );
	if ( strncmp("mrtrix tracks",buffer,13)!=0 )
	{
		COLOR_error("The file is not a valid MRTRIX .tck file");
		fclose(fp);
		return EXIT_FAILURE;
	}

	string key, value;
	size_t pos;
	while ( 1 )
	{
		fscanf( fp, "%100[^\n]\n", buffer );
		if ( strncmp("END",buffer,3)==0 ) break;

		key = buffer;
		pos = key.find( ":" );
		value = key.substr( pos+1 );
		key = key.substr(0, pos );

		trim(key);
		trim(value);

		HDR[ key ] = value;
		printf( "\t - %s : '%s'\n", key.c_str(), value.c_str() );
	}
	cout <<COLOR(37,37,1)<< "   [OK]\n\n" <<COLOR_reset;


	// checking input data
	if ( HDR.count( "file" )==0 )
	{
		COLOR_error("The field 'file' is missing in the header");
		return EXIT_FAILURE;
	}
	;
	int FileOffset = atoi( HDR["file"].substr( HDR["file"].find_last_of( ' ' )+1 ).c_str());

	if ( HDR.count( "count" )==0 )
	{
		COLOR_error("The field 'count' is missing in the header");
		return EXIT_FAILURE;
	}
	int nFIBERS = atoi( HDR["count"].c_str() );


/*------------------------------------------------*/
/*  Reading the NIFTI header to use as reference  */
/*------------------------------------------------*/
	cout <<COLOR(37,37,1)<< "-> Reading 'REFERENCE' dataset...\n" <<COLOR_reset;

	NIFTI<UINT8> niiREF( REF_filename, true );
	if ( !niiREF.isValid() ) {
		if ( niiREF.getErrorCode() == NIFTI_ERROR_WRONGDATATYPE )
			COLOR_error( "Datatype should be UINT8!" );
		else
			COLOR_error( "Unable to open file!" );
		return EXIT_FAILURE;
	}

	printf("\t - dim   : %d x %d x %d x %d\n", niiREF.hdr->dim[1],niiREF.hdr->dim[2],niiREF.hdr->dim[3],niiREF.hdr->dim[4]);
	printf("\t - pixdim: %.4f x %.4f x %.4f\n", niiREF.hdr->pixdim[1],niiREF.hdr->pixdim[2],niiREF.hdr->pixdim[3]);

	// Reading QFORM
	printf("\t - qform : ");
	float REF_M[4][4];
	if ( niiREF.hdr->qform_code )
	{
		mat44 A = nifti_mat44_inverse( niiREF.hdr->qto_xyz ); //mat44 B = niiREF.hdr->qto_ijk;
		for(int j=0; j<4 ;j++)
		for(int i=0; i<4 ;i++)
			REF_M[i][j] = A.m[i][j];

		printf("\n");
		printf("\t\t%8.3f %8.3f %8.3f %8.3f\n", REF_M[0][0], REF_M[1][0], REF_M[2][0], REF_M[3][0] );
		printf("\t\t%8.3f %8.3f %8.3f %8.3f\n", REF_M[0][1], REF_M[1][1], REF_M[2][1], REF_M[3][1] );
		printf("\t\t%8.3f %8.3f %8.3f %8.3f\n", REF_M[0][2], REF_M[1][2], REF_M[2][2], REF_M[3][2] );
		printf("\t\t%8.3f %8.3f %8.3f %8.3f\n", REF_M[0][3], REF_M[1][3], REF_M[2][3], REF_M[3][3] );
	}
	else
	{
		COLOR_error( "not found!\n" );
		return EXIT_FAILURE;
	}
	cout <<COLOR(37,37,1)<< "   [OK]\n\n" <<COLOR_reset;



/*--------------------------------------*/
/*  Reading and converting streamlines  */
/*--------------------------------------*/
	cout <<COLOR(37,37,1)<< "-> Copying "<< nFIBERS <<" fibers...\n" <<COLOR_reset;

	short int TRK_dim[3] = { niiREF.hdr->dim[1],niiREF.hdr->dim[2],niiREF.hdr->dim[3] };
	float TRK_pixdim[3]  = { niiREF.hdr->pixdim[1],niiREF.hdr->pixdim[2],niiREF.hdr->pixdim[3] };

	TrackVis TRK = TrackVis();

	if ( !TRK.create( OUT_filename, TRK_dim, TRK_pixdim ) )
	{
		COLOR_error("Unable to create the .trk file");
		return EXIT_FAILURE;
	}
	sprintf(TRK.hdr.voxel_order,"RAS");
	sprintf(TRK.hdr.pad2,"RAS");
	TRK.writeHdr();

	fseek( fp, FileOffset, SEEK_SET );

	float xyz[3];
	int len, totStored = 0;
	blitz::Array<float,2> fiber(3,10000);	// tmp fiber to store coordinates

	for(int i=0; i<nFIBERS ;i++)
	{
		len = 0;
		while( 1 )
		{
			fread( xyz, 4, 3, fp );
			if ( xyz[0] != xyz[0] || xyz[1] != xyz[1] || xyz[2] != xyz[2] )		// end of fiber
				break;
			fiber(0,len) = TRK_pixdim[0] * ( REF_M[0][0]*xyz[0] + REF_M[0][1]*xyz[1] + REF_M[0][2]*xyz[2] + REF_M[0][3] );
			fiber(1,len) = TRK_pixdim[1] * ( REF_M[1][0]*xyz[0] + REF_M[1][1]*xyz[1] + REF_M[1][2]*xyz[2] + REF_M[1][3] );
			fiber(2,len) = TRK_pixdim[2] * ( REF_M[2][0]*xyz[0] + REF_M[2][1]*xyz[1] + REF_M[2][2]*xyz[2] + REF_M[2][3] );
			len++;
		}

		TRK.append( &fiber, len, TRACKVIS_SAVE_ALL );
		totStored++;

		printf("\r   [ %5.1f%% ]", 100.0 * i / nFIBERS);
	}
	cout <<COLOR(37,37,1)<< "\r   [OK]      \n\n" <<COLOR_reset;

	fclose(fp);
	TRK.updateTotal( totStored );


	cout <<COLOR(37,37,1)<< "[ "<< totStored <<" fibers written ]\n\n" <<COLOR_reset;
	return 0;
}
