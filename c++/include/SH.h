#ifndef __SPHERICAL_HARMONICS_H__
#define __SPHERICAL_HARMONICS_H__

#include <blitz/array.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>



namespace SH
{


inline int index(int l, int m) { return (l*(l+1)/2 + m); }


/**
	Aux function to multiply blitz++ Arrays.
*/
// #define 	mult(A,B) 	sum( A(blitz::tensor::i,blitz::tensor::j) * B(blitz::tensor::j,blitz::tensor::k), blitz::tensor::j );

template <typename T>
inline blitz::Array<T,2> mult( blitz::Array<T,2> & A, blitz::Array<T,2> & B )
{
	return (blitz::Array<T,2>)sum( A(blitz::tensor::i,blitz::tensor::k) * B(blitz::tensor::k,blitz::tensor::j), blitz::tensor::k );
}

template <typename T>
inline blitz::Array<T,1> mult( blitz::Array<T,2> & FIT, blitz::Array<T,1> & B )
{
	return (blitz::Array<T,1>)( sum( FIT(blitz::tensor::i,blitz::tensor::j) * B(blitz::tensor::j), blitz::tensor::j ) );
}


const unsigned char		BASIS_DESCOTEAUX	=	1;
const unsigned char		BASIS_MRTRIX		=	2;



/**
	Compute SH basis.

	@param[in]		lmax		Spherical Harmonics order
	@param[in]		dirs		Array with the directions on the sphere [ndirs x 3]
	@param[out]		Ylm			Spherical Harmonics basis [ndirs x (lmax+1)*(lmax+2)/2]
	@param[in]		type		[BASIS_DESCOTEAUX, BASIS_MRTRIX]
*/
template <typename T>
short createBasis( int lmax, const blitz::Array<T,2> & dirs, blitz::Array<T,2> & Ylm, unsigned char type = BASIS_DESCOTEAUX )
{
	if ( lmax<0 || lmax%2==1 )
		return 0;
	if ( dirs.extent(0)<1 || dirs.extent(1)!=3 )
		return 0;

	int ndirs = dirs.extent(0), ncoeffs = (lmax+1)*(lmax+2)/2;

	// compute the spherical coordinates of the directions
	blitz::Array<T,1> el(ndirs), az(ndirs);
	T x, y, z, r;
	int d;
	for(d=0; d<ndirs ;d++)
	{
		x = dirs(d,0);
		y = dirs(d,1);
		z = dirs(d,2);
		r = sqrt(x*x + y*y + z*z);
		az(d) = atan2( y/r, x/r );
		el(d) = acos( z/r );
	}


	// create the matrix of the SH basis
	Ylm.resize( ndirs, ncoeffs );
	int l, m, j = 0;
// 	double Plm[lmax+1], dPlm[lmax+1];
// 	float cEL, cAZ, sAZ;

	switch ( type )
	{
		// Descoteaux et al (2007)
		case BASIS_DESCOTEAUX:
		for( d=0; d < ndirs ;d++ )
		{
			for( l = 0, j = 0; l <= lmax ; l+=2 )
			for( m = -l; m <= l ;m++, j++ )
				if ( m > 0 )
					Ylm(d,j) = M_SQRT2 * gsl_sf_legendre_sphPlm( l, m, cos(el(d)) ) * sin(m*az(d));
				else if( m < 0 )
					Ylm(d,j) = M_SQRT2 * gsl_sf_legendre_sphPlm( l, -m, cos(el(d)) ) * cos(-m*az(d));
				else
					Ylm(d,j) = gsl_sf_legendre_sphPlm( l, 0, cos(el(d)) );
		}
		break;

		// mrtrix software
		case BASIS_MRTRIX:
		for( d = 0; d < ndirs ;d++ )
		for( l = 0; l <= lmax ; l+=2 )
		for( m = 0; m <= l ;m++ )
		{
			T s = gsl_sf_legendre_sphPlm( l, m, cos(el(d)) );
			if ( m )
			{
                Ylm(d,index(l, m)) = s * cos(m*az(d));
                Ylm(d,index(l,-m)) = s * sin(m*az(d));
			}
			else
				Ylm(d,index(l,0)) = s;
		}
		break;


		// testing (not to be used)
// 		case -1:
// 		for( d=0; d < ndirs ;d++ )
// 		{
// 			cEL = cos(el(d));
//
// 			// m = 0
// 			gsl_sf_legendre_sphPlm_deriv_array( lmax, 0, cEL, Plm, dPlm );
// 			for (l = 0; l <= lmax; l+=2)
//  				Ylm(d,index(l,0)) = Plm[l];
//
// 			// m = 1
//
// 			// m > 1
// 			for( m = 1; m <= lmax ;m++ )
// 			{
//  				gsl_sf_legendre_sphPlm_deriv_array( lmax, m, cEL, Plm, dPlm );
// 				cAZ = cos(m*az(d));
// 				sAZ = sin(m*az(d));
// 				for (l = ((m&1) ? m+1 : m); l <= lmax; l+=2)
// 				{
// 					Ylm(d,index(l, m)) = Plm[l-m] * cAZ;
// 					Ylm(d,index(l,-m)) = Plm[l-m] * sAZ;
// 				}
// 			}
// 		}
// 		break;

		default:
		return 0;
	}

	return 1;
}


/**
	Pre-compute the matrix needed to fit the SH coeff of a function

	@param[in]		Ylm			Spherical Harmonics basis
	@param[in]		lambda		Regularization term
	@param[out]		FIT			Matrix to estimate the SH coeff of the signal
*/
template <typename T>
short precomputeFitMatrix( blitz::Array<T,2> & Ylm, double lambda, blitz::Array<T,2> & FIT )
{
	int i, j, l, m;
	int nD = Ylm.extent(0);
	int nC = Ylm.extent(1);
	int lmax = 0.5 * ( -3 + sqrt(1+8*nC) );

	FIT.resize( nC, nD );

	blitz::Array<T,2> B( nC, nC );
	blitz::Array<T,2> YlmT( nC, nD );

	YlmT = Ylm.transpose(1,0);
	B = mult( YlmT, Ylm );

	T diag1, diag2, diag3;
	for( l = 0, j=0; l<=lmax ; l+=2 )
	{
		diag1 = lambda * l*l * (l+1)*(l+1);
		for( m = -l; m<=l ;m++, j++ )
			B(j,j)   += diag1;
	}

	// compute the FIT matrix
	gsl_matrix* M = gsl_matrix_alloc( nC, nC );
	gsl_matrix* Minv = gsl_matrix_alloc( nC, nC );
	for(i = 0; i < nC; i++)
	for(j = 0; j < nC; j++)
		gsl_matrix_set( M, i, j, B(i,j) );

	gsl_permutation* p = gsl_permutation_alloc( nC );
	gsl_linalg_LU_decomp( M, p, &i);
	gsl_linalg_LU_invert( M, p, Minv);

	for(i = 0; i < nC; i++)
	for(j = 0; j < nC; j++)
		B(i,j) = gsl_matrix_get( Minv, i, j );

	FIT = mult( B, YlmT );
	return 1;
}


}

#endif
