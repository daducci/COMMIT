

void COMMIT_L(
    int nF, int nIC, int nV, int nS,
    double *vIN, double *vOUT)
{
    for(int f = 0; f < nF; f++){

        vOUT[nV*nS] += -2*vIN[f] + x[nF + f];

        for(int r = 1; r < nIC-1; r++){
            vOUT[nV*nS + r] += vIN[(r-1)*nF + f] -2*vIN[r*nF + f] + vIN[(r+1)*nF + f];
        }

        vOUT[nV*nS + nIC - 1] += vIN[(nIC-2)*nF + f] - 2*vIN[(nIC-1)*nF + f];
    }
}

void COMMIT_Lt(
    int nF, int nIC, int nV, int nS,
    double *vIN, double *vOUT)
{
    
}