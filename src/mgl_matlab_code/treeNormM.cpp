#include "mgl.h"
/* global screening for fused multiple graphical lasso
Author: Sen Yang
Mar 23 2013, Arizona State University*/
void mexFunction (int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
	const mwSize *sDims;
	double *G;
    double *ind = NULL;
	size_t M,N,K,nodes;
    double Num;
	sDims = mxGetDimensions(prhs[0]);

	G = mxGetPr(prhs[0]);

    ind = mxGetPr(prhs[1]);
    nodes = (size_t) mxGetScalar(prhs[2]);

	M =(size_t) sDims[0];
	N =(size_t) sDims[1];
	K = (size_t) sDims[2];
	Num = treeNormM(G,(int)M, (int)K,ind,(int)nodes);
	plhs[0] = mxCreateDoubleScalar((double)Num);
}
