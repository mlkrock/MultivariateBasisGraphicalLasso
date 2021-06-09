#include "mgl.h"
/* global screening for fused multiple graphical lasso
Author: Sen Yang
Mar 23 2013, Arizona State University*/
void mexFunction (int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
	const mwSize *sDims;
	double *G, *adj;
    double *ind = NULL;
	double lam, rho;
	size_t Num,M,N,K, penalty,nodes;
	sDims = mxGetDimensions(prhs[0]);

	G = mxGetPr(prhs[0]);

	lam = mxGetScalar(prhs[1]);
    rho = mxGetScalar(prhs[2]);
    penalty = (size_t) mxGetScalar(prhs[3]);
    
    if(penalty==2){
        ind = mxGetPr(prhs[4]);
        nodes = (size_t) mxGetScalar(prhs[5]);
    }

	M =(size_t) sDims[0];
	N =(size_t) sDims[1];
	K = (size_t) sDims[2];
    plhs[0] = mxCreateDoubleMatrix(M, M, mxREAL);
    adj = mxGetPr(plhs[0]);
	Num = globalScreening(G, adj, (int)M, (int)K, lam, rho, (int)penalty,ind,(int)nodes);
	plhs[1] = mxCreateDoubleScalar((double)Num);
}
