#include "mgl.h"
/* proximal opertor of multiple graphical lasso
 * Author: Sen Yang
 * Sep. 11, 2013
 * Arizona State University
 */ 
void mexFunction (int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
	const mwSize *sDims;
	double *input, *output, *fx, *fy;
	double lam, rho;
    double *ind = NULL;
	size_t Num,M,N,K, penalty,nodes=0;
	sDims = mxGetDimensions(prhs[0]);

	output = mxGetPr(prhs[0]);
	input = mxGetPr(prhs[1]);
	fx = mxGetPr(prhs[2]);
	fy = mxGetPr(prhs[3]);

	lam = mxGetScalar(prhs[4]);
    rho = mxGetScalar(prhs[5]);
    Num = (size_t) mxGetScalar(prhs[6]);
    penalty = (size_t) mxGetScalar(prhs[7]);
    
    if(penalty==2){
        ind = mxGetPr(prhs[8]);
        nodes = (size_t) mxGetScalar(prhs[9]);
    }
    
	M =(size_t) sDims[0];
	N =(size_t) sDims[1];
	K = (size_t) sDims[2];
	mgl_prox(output, input,  fx, fy, lam, rho, sDims[0], sDims[2], Num, penalty,ind,(int)nodes);

}