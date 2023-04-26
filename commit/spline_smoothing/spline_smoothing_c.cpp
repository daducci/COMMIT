#include "Catmull.h"
#include "psimpl.h"
#include "Vector.h"
#include <vector>
#include <iostream>


// =========================
// Function called by CYTHON
// =========================

int do_spline_smoothing( float* ptr_npaFiberI, int nP, float* ptr_npaFiberO, float ratio, float segment_len)
{
    std::vector<float>          polyline_simplified;
    std::vector<Vector<float>>  CPs;
    Catmull                     FIBER;
    

    // simplify input polyline down to nP*ratio points
    psimpl::simplify_douglas_peucker_n <3>( ptr_npaFiberI, ptr_npaFiberI+3*nP, nP*ratio, std::back_inserter(polyline_simplified) );
    // psimpl::simplify_douglas_peucker <3> ( ptr_npaFiberI, ptr_npaFiberI+3*nP, tolerance, std::back_inserter(polyline_simplified) );
    CPs.resize( polyline_simplified.size()/3 );
    for( int j=0,index=0; j < polyline_simplified.size(); j=j+3 )
        CPs[index++].Set( polyline_simplified[j], polyline_simplified[j+1], polyline_simplified[j+2] );


    // perform interpolation
    FIBER.set( CPs );
    FIBER.eval( FIBER.L/segment_len );
    FIBER.arcLengthReparametrization( segment_len );
    // copy coordinates of the smoothed streamline back to python
    for( int j=0; j < FIBER.P.size(); j++ )
    {
        *(ptr_npaFiberO++) = FIBER.P[j].x;
        *(ptr_npaFiberO++) = FIBER.P[j].y;
        *(ptr_npaFiberO++) = FIBER.P[j].z;
    }
    return FIBER.P.size();
}

int compute_control_p(float* ptr_origFib, int nP, float* ptr_knotsFib, float ratio, float segment_len, float tolerance) {
    std::vector<float>          polyline_simplified;
    std::vector<Vector<float>>  CPs;
    std::vector<Vector<float> >	KNOTs;
    Catmull                     FIBER;
    

    // simplify input polyline down to nP*ratio points
    psimpl::simplify_douglas_peucker_n <3>( ptr_origFib, ptr_origFib+3*nP, nP*ratio, std::back_inserter(polyline_simplified) );
    // psimpl::simplify_douglas_peucker <3> ( ptr_npaFiberI, ptr_npaFiberI+3*nP, tolerance, std::back_inserter(polyline_simplified) );

    CPs.resize( polyline_simplified.size()/3 );
    for( int j=0,index=0; j < polyline_simplified.size(); j=j+3 ){
        CPs[index++].Set( polyline_simplified[j], polyline_simplified[j+1], polyline_simplified[j+2] );
    }
    
    KNOTs.resize(CPs.size()+2);
    KNOTs[0] = CPs[0];
    for( int i = 0; i< CPs.size(); i++)
        KNOTs[i+1] = CPs[i];
    KNOTs.back() = CPs[CPs.size()-1];

    // copy coordinates of the knots back to python
    for( int j=0; j < KNOTs.size(); j++ ){
        *(ptr_knotsFib++) = KNOTs[j].x;
        *(ptr_knotsFib++) = KNOTs[j].y;
        *(ptr_knotsFib++) = KNOTs[j].z;
    }
    
    return KNOTs.size();
}
