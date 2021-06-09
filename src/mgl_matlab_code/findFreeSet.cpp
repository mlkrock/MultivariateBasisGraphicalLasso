#include "mgl.h"
/* find free set 
Author: Sen Yang
Mar 23 2013, Arizona State University*/
void mexFunction (int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
	const mwSize *sDims;
	double *G, *P, *v, *fx, *fy;
    double *ind = NULL;
	double lam, rho;
	size_t Num, M, N, K, penalty,nodes=0;

	sDims = mxGetDimensions(prhs[0]);

	G = mxGetPr(prhs[0]);
	P = mxGetPr(prhs[1]);
	fx = mxGetPr(prhs[2]);
	fy = mxGetPr(prhs[3]);

	lam = mxGetScalar(prhs[4]);
    rho = mxGetScalar(prhs[5]);
    v = mxGetPr(prhs[6]);
	penalty = (size_t) mxGetScalar(prhs[7]);
    
    if(penalty==2){
        ind = mxGetPr(prhs[8]);
        nodes = (size_t) mxGetScalar(prhs[9]);
    }
    

	M =(size_t) sDims[0];
	N =(size_t) sDims[1];
	K = (size_t) sDims[2];

	Num = findFreeSet(G, P,  v, fx, fy, (int)M, (int)K, lam, rho,0.00001, (int)penalty,ind,(int)nodes);
	plhs[0] = mxCreateDoubleScalar((double)Num);
}

