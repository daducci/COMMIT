#ifndef __Vector_H__
#define __Vector_H__

#include <math.h>


template <class T>
class Vector
{
    public:
        T 		x, y, z;

        void 	Set( T _x, T _y, T _z );
        float 	norm();
        void 	Normalize();
        void	Multiply( T k );
        void 	VectorProduct( const Vector<T> & v1, const Vector<T> & v2 );
        T		ScalarProduct( const Vector<T> & v );
        float	DistanceTo( const Vector<T> & v );

        void    operator=( const Vector<T> & v );

        Vector();
        Vector( T _x, T _y, T _z );
};


template <class T>
Vector<T>::Vector()
{ x = y = z = 0; }


template <class T>
Vector<T>::Vector( T _x, T _y, T _z )
{ x = _x; y = _y; z = _z; }


template <class T>
void Vector<T>::Set( T _x, T _y, T _z )
{ x = _x; y = _y; z = _z; }


template <class T>
float Vector<T>::norm()
{
    return sqrt(x*x+y*y+z*z);
}


template <class T>
void Vector<T>::Normalize()
{
    float len = sqrt(x*x+y*y+z*z);
    if (len==0) return;
    x = x/len;
    y = y/len;
    z = z/len;
}


template <class T>
void Vector<T>::Multiply( T k )
{
    x *= k;
    y *= k;
    z *= k;
}


template <class T>
T Vector<T>::ScalarProduct( const Vector<T> & v )
{
    return x*v.x + y*v.y + z*v.z;
}


template <class T>
void Vector<T>::VectorProduct( const Vector<T> & v1, const Vector<T> & v2 )
{
    x = v1.y*v2.z - v2.y*v1.z;
    y = v1.z*v2.x - v2.z*v1.x;
    z = v1.x*v2.y - v2.x*v1.y;
}

template <class T>
float Vector<T>::DistanceTo( const Vector<T> & v )
{
    return sqrt( pow(x-v.x,2) + pow(y-v.y,2) + pow(z-v.z,2) );
}

template <class T>
void Vector<T>::operator=( const Vector<T> & v )
{
    x = v.x;
    y = v.y;
    z = v.z;
}

#endif
