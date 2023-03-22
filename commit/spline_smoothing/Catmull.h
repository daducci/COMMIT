#ifndef __CATMULL_H__
#define __CATMULL_H__

#include "Vector.h"
#include <vector>


/* ------------------------------------------------------------------------------------------------------ */
/*                                    Class to manage Catmull curves                                       */
/* ------------------------------------------------------------------------------------------------------ */
class Catmull
{
    public:

    std::vector<Vector<float> >	KNOTs;
    std::vector<Vector<float> >	P, dP, ddP;
    float				        L, weight, curv;

    Catmull(): L(0), weight(0), curv(0) { }
    Catmull(const Catmull& c) {L = c.L; weight = c.weight; curv = c.curv; KNOTs = c.KNOTs; P = c.P; dP = c.dP; ddP = c.ddP; }

    void		set(std::vector<Vector<float> >);
    float       degrees( int );
    float		approxLen( );
    void		eval( int );
    void		arcLengthReparametrization( float );
    float       curvature ( int );
    void 		MoveKnot(int, POINT, float );
    void		Translate(int, int, POINT, float );
    void 		MoveKnotOrthogonal(int,int, float);
};


/* =============================== */
/* Set the position of the KNOTs   */
/* =============================== */
void Catmull::set(std::vector<Vector<float> > points)
{
    KNOTs.resize(points.size()+2);
    KNOTs[0] = points[0];
    for( int i = 0; i< points.size(); i++)
        KNOTs[i+1] = points[i];
    KNOTs.back() = points[points.size()-1];
    L = approxLen();
}

/* ================================ */
/* Evaluate the curve as a polyline */
/* ================================ */
void Catmull::eval( int nSegments )
{
    float subdiv_step, subdiv_step2, subdiv_step3;
    float ax, ay, az;
    float bx, by, bz;
    float cx, cy, cz;
    float fx, fy, fz;
    float dfx,  dfy,  dfz;
    float ddfx, ddfy, ddfz;
    float dddfx, dddfy, dddfz;

    int idx = 0;

    if ( nSegments<1 ) return;

    P.resize(  nSegments*(KNOTs.size()-3)+1 );
    dP.resize( nSegments*(KNOTs.size()-3)+1 );
    ddP.resize( nSegments*(KNOTs.size()-3)+1 );
    for( int k=0; k<(KNOTs.size()-3); k++ )
    {
        ax = KNOTs[k+3].x - 3.0 * KNOTs[k+2].x + 3.0 * KNOTs[k+1].x - KNOTs[k+0].x;
        ay = KNOTs[k+3].y - 3.0 * KNOTs[k+2].y + 3.0 * KNOTs[k+1].y - KNOTs[k+0].y;
        az = KNOTs[k+3].z - 3.0 * KNOTs[k+2].z + 3.0 * KNOTs[k+1].z - KNOTs[k+0].z;

        bx = 2.0 * KNOTs[k+0].x - 5.0 * KNOTs[k+1].x + 4.0 * KNOTs[k+2].x - KNOTs[k+3].x;
        by = 2.0 * KNOTs[k+0].y - 5.0 * KNOTs[k+1].y + 4.0 * KNOTs[k+2].y - KNOTs[k+3].y;
        bz = 2.0 * KNOTs[k+0].z - 5.0 * KNOTs[k+1].z + 4.0 * KNOTs[k+2].z - KNOTs[k+3].z;

        cx = KNOTs[k+2].x - KNOTs[k+0].x;
        cy = KNOTs[k+2].y - KNOTs[k+0].y;
        cz = KNOTs[k+2].z - KNOTs[k+0].z;

        subdiv_step  = 1.0 / float(nSegments);
        subdiv_step2 = subdiv_step * subdiv_step;
        subdiv_step3 = subdiv_step * subdiv_step2;

        dfx =  0.5 * subdiv_step3 * ax + 0.5 * subdiv_step2 * bx + 0.5 * subdiv_step * cx;
        dfy =  0.5 * subdiv_step3 * ay + 0.5 * subdiv_step2 * by + 0.5 * subdiv_step * cy;
        dfz =  0.5 * subdiv_step3 * az + 0.5 * subdiv_step2 * bz + 0.5 * subdiv_step * cz;

        ddfx = 3.0 * subdiv_step3 * ax + subdiv_step2 * bx;
        ddfy = 3.0 * subdiv_step3 * ay + subdiv_step2 * by;
        ddfz = 3.0 * subdiv_step3 * az + subdiv_step2 * bz;

        dddfx = 3.0 * subdiv_step3 * ax;
        dddfy = 3.0 * subdiv_step3 * ay;
        dddfz = 3.0 * subdiv_step3 * az;

        fx = KNOTs[k+1].x;
        fy = KNOTs[k+1].y;
        fz = KNOTs[k+1].z;

        for(int i=0; i<nSegments ;i++)
        {
            P[idx].x = fx;
            P[idx].y = fy;
            P[idx].z = fz;
            dP[idx].x = dfx;
            dP[idx].y = dfy;
            dP[idx].z = dfz;
            ddP[idx].x = ddfx;
            ddP[idx].y = ddfy;
            ddP[idx].z = ddfz;
            idx++;
            fx += dfx; dfx += ddfx; ddfx += dddfx;
            fy += dfy; dfy += ddfy; ddfy += dddfy;
            fz += dfz; dfz += ddfz; ddfz += dddfz;
        }
    }

    P[idx].x = KNOTs[KNOTs.size()-2].x;
    P[idx].y = KNOTs[KNOTs.size()-2].y;
    P[idx].z = KNOTs[KNOTs.size()-2].z;
    dP[idx].x = dfx;
    dP[idx].y = dfy;
    dP[idx].z = dfz;
    ddP[idx].x = ddfx;
    ddP[idx].y = dfy;
    ddP[idx].z = dfz;
}

/* =============================================================================== */
/* Calculate the curvature (restriction)                                           */
/* Source: Efficient Energies and Algorithms for parametric snakes, Mathews Jacob  */
/* =============================================================================== */
float Catmull::curvature( int nSegments )
{
    float c;
    curv = 0;
    for(int i=0;i<nSegments;i++)
        curv=curv+(sqrt(pow(ddP[i].z*dP[i].y - ddP[i].y*dP[i].z,2) + pow(ddP[i].x*dP[i].z - ddP[i].z*dP[i].x,2) + pow(ddP[i].x*dP[i].y - ddP[i].y*dP[i].x,2)))/(sqrt( pow(pow(dP[i].x,2) + pow(dP[i].y,2) + pow(dP[i].z,2),3)));
    curv = curv / ( nSegments );
    return curv;
}

/* ============================================================== */
/* Estimate an approximation of the length of the curve           */
/* ============================================================== */
float Catmull::approxLen( )
{
    float length = 0;
    for(int i = 0; i< KNOTs.size()-3; i++)
        length+= KNOTs[i+1].DistanceTo( KNOTs[i+2] );
    return length;
}

/* ============================================================== */
/* Reparametrize the curve by the arc-length                      */
/* ============================================================== */
void Catmull::arcLengthReparametrization( float SEGMENT_len )
{
    int N = P.size();
    int ii = 0;
    float* arcLength = new float[N];
    std::vector<Vector<float> > P_old;
    std::vector<Vector<float> > dP_old;
    std::vector<Vector<float> > ddP_old;

    dP_old = dP;
    P_old  = P;
    ddP_old = ddP;

    //Compute the total length of arc
    arcLength[0] = 0;
    for( int i=1; i<N; i++)
        arcLength[i] = arcLength[i-1] + P_old[i].DistanceTo( P_old[i-1] );
    this->L = arcLength[N-1];
    int N2 = ceil( this->L / SEGMENT_len );
    float delta = this->L / (N2-1);
    P.resize( N2 );
    dP.resize( N2 );
    ddP.resize( N2 );
    //Redefine the points
    int index = 0, found;
    float s, k1, k2, k3;

    for(int i = 0; i< P.size(); i++)
    {
        s = delta * ii;
        ii++;
        found = 0;
        for (; index < (N-1) && !found; index++)
            if ( arcLength[index] <= s && arcLength[index+1] >= s )
                found = 1;
        index--;
        k1 = arcLength[index+1] - s;
        k2 = s - arcLength[index];
        k3 = arcLength[index+1] - arcLength[index];

        P[i].x  = ( k1 *  P_old[index].x + k2 *  P_old[index+1].x ) / k3;
        P[i].y  = ( k1 *  P_old[index].y + k2 *  P_old[index+1].y ) / k3;
        P[i].z  = ( k1 *  P_old[index].z + k2 *  P_old[index+1].z ) / k3;
        dP[i].x = ( k1 * dP_old[index].x + k2 * dP_old[index+1].x ) / k3;
        dP[i].y = ( k1 * dP_old[index].y + k2 * dP_old[index+1].y ) / k3;
        dP[i].z = ( k1 * dP_old[index].z + k2 * dP_old[index+1].z ) / k3;

        ddP[i].x = ( k1 * ddP_old[index].x + k2 * ddP_old[index+1].x ) / k3;
        ddP[i].y = ( k1 * ddP_old[index].y + k2 * ddP_old[index+1].y ) / k3;
        ddP[i].z = ( k1 * ddP_old[index].z + k2 * ddP_old[index+1].z ) / k3;
    }

    delete arcLength;
}

void Catmull::MoveKnot(int knot, POINT move_direction, float magnitude) {
	POINT direction;
	direction.Set(move_direction.x, move_direction.y, move_direction.z);
	direction.Multiply(magnitude);
	KNOTs[knot].Add(direction);
}

void Catmull::Translate(int StartKnot, int EndKnot, POINT move_direction, float magnitude) {
	POINT direction;
	direction.Set(move_direction.x, move_direction.y, move_direction.z);
	direction.Multiply(magnitude);
	for(int knot = StartKnot; knot <= EndKnot; knot++) {
		this->MoveKnot(knot, direction, magnitude);
	}
}

void Catmull::MoveKnotOrthogonal(int StartKnot,int EndKnot, float magnitude) {								//cannot be edge knots!
	std::vector<POINT> OldKnots;
	POINT baseline;
	POINT knotvector;
	POINT straightenvector;
	OldKnots.resize(this->KNOTs.size());
	for (int i = 0; i< KNOTs.size(); i++) {
		OldKnots[i].Set(KNOTs[i].x,KNOTs[i].y,KNOTs[i].z);
	}
	for( int knot = StartKnot; knot <= EndKnot; knot++) {
		baseline.Set(OldKnots[knot+1].x - OldKnots[knot-1].x,  OldKnots[knot+1].y - OldKnots[knot-1].y, OldKnots[knot+1].z - OldKnots[knot-1].z);
		knotvector.Set(OldKnots[knot].x - OldKnots[knot-1].x, OldKnots[knot].y - OldKnots[knot-1].y, OldKnots[knot].z - OldKnots[knot-1].z);
		straightenvector.Perp_Projection(knotvector, baseline);
		straightenvector.Normalize();
		this->MoveKnot(knot,straightenvector,magnitude);
		
	}
}

#endif
