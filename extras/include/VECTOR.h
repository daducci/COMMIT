#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <math.h>	// [CHECK]

template <class T>
class VECTOR
{
	public:
		T 		x, y, z;

		void 	Set( T _x, T _y, T _z );
		float 	norm();
		void 	Normalize();
		void	Multiply( T k );
		void 	VectorProduct( const VECTOR<T> & v1, const VECTOR<T> & v2 );
		T		ScalarProduct( const VECTOR<T> & v );
		float	DistanceTo( const VECTOR<T> & v );

        void    operator=( const VECTOR<T> & v );

		VECTOR();
		VECTOR( T _x, T _y, T _z );
};


template <class T>
VECTOR<T>::VECTOR()
{ x = y = z = 0; }


template <class T>
VECTOR<T>::VECTOR( T _x, T _y, T _z )
{ x = _x; y = _y; z = _z; }


template <class T>
void VECTOR<T>::Set( T _x, T _y, T _z )
{ x = _x; y = _y; z = _z; }


template <class T>
float VECTOR<T>::norm()
{
	return sqrt(x*x+y*y+z*z);
}


template <class T>
void VECTOR<T>::Normalize()
{
	float len = sqrt(x*x+y*y+z*z);
	if (len==0) return;
	x = x/len;
	y = y/len;
	z = z/len;
}


template <class T>
void VECTOR<T>::Multiply( T k )
{
	x *= k;
	y *= k;
	z *= k;
}


template <class T>
T VECTOR<T>::ScalarProduct( const VECTOR<T> & v )
{
	return x*v.x + y*v.y + z*v.z;
}


template <class T>
void VECTOR<T>::VectorProduct( const VECTOR<T> & v1, const VECTOR<T> & v2 )
{
	x = v1.y*v2.z - v2.y*v1.z;
	y = v1.z*v2.x - v2.z*v1.x;
	z = v1.x*v2.y - v2.x*v1.y;
}

template <class T>
float VECTOR<T>::DistanceTo( const VECTOR<T> & v )
{
	return sqrt( pow(x-v.x,2) + pow(y-v.y,2) + pow(z-v.z,2) );
}

template <class T>
void VECTOR<T>::operator=( const VECTOR<T> & v )
{
    x = v.x;
    y = v.y;
    z = v.z;
}


typedef VECTOR<float> 		POINT;


#endif
