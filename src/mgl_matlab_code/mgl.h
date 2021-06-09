#include <math.h>
#include <time.h>			
#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include "matrix.h"
#include <cstring>
#include "altra.h"
#include<iostream>
/*
This include sub functions for multiple graphical lasso and screening
 - fused multiple graphical lasso
 - group multiple graphical lasso
 - tree multiple graphical lasso 

 Author: Sen Yang
 Email: senyang@asu.edu
 Sep. 11, 2013
 Arizona State University
*/

/*screening scripts
  - sub function for fused multiple graphical lasso
*/
bool screeningFused(double* v, const int K, const double lam, const double rho, const double EPS){
	double sumT = 0;
	double threshold = rho;
	for (int k = 0; k < K - 1; k++){
		sumT += v[k];
		threshold += lam;
		if (fabs(sumT) + EPS > threshold)
			return false;
	}

	sumT += v[K - 1];
	threshold += lam - rho;
	if (fabs(sumT) + EPS > threshold)
		return false;

	for (int k = 1; k < K - 1; k++)
		for (int j = k; j < K - 1; j++){
			sumT = 0;
			threshold = (j - k + 1)*lam + 2 * rho;
			for (int i = k; i < j + 1; i++)
				sumT += v[i];
			if (fabs(sumT) + EPS > threshold)
				return false;
		}

	sumT = 0;
	threshold = rho;
	for (int k = K - 1; k > 0; k--){
		sumT += v[k];
		threshold += lam;
		if (fabs(sumT) + EPS > threshold)
			return false;
	}
	return true;
}

/*screening scripts
  - sub function for group multiple graphical lasso
*/
bool screeningGroup(double* v, const int K, const double lam, const double rho, const double EPS){
	double temp = 0;
	for (int k = 0; k < K; k++){
		if (v[k] < -lam)
			v[k] += lam;
		else if (v[k] > lam)
			v[k] -= lam;
		else
			v[k] = 0;
		temp += v[k] * v[k];
	}
	return temp + EPS <= rho * rho;
}

/*screening scripts
- sub function for tree multiple graphical lasso
*/
bool screeningTree(double* v, const int n, const double *ind, const int nodes, const double EPS){

	int i, j, m;
	double lambda, twoNorm, ratio;
	/*
	* test whether the first node is special
	*/
	if ((int) ind[0] == -1){

		/*
		*Recheck whether ind[1] equals to zero
		*/
		if ((int) ind[1] != -1){
			printf("\n Error! \n Check ind");
			exit(1);
		}

		lambda = ind[2]-EPS;

		for (j = 0; j<n; j++){
			if (v[j]>lambda)
				v[j] -= lambda;
			else
				if (v[j] < -lambda)
					v[j] += lambda;
				else
					v[j] = 0;
		}

		i = 1;
	}
	else{
		i = 0;
	}

	/*
	* sequentially process each node
	*
	*/
	for (; i < nodes; i++){
		/*
		* compute the L2 norm of this group
		*/
		twoNorm = 0;
		for (j = (int) ind[3 * i] - 1; j < (int) ind[3 * i + 1]; j++)
			twoNorm += v[j] * v[j];
		twoNorm = sqrt(twoNorm);

		lambda = ind[3 * i + 2] - EPS; 
		if (twoNorm > lambda){
			ratio = (twoNorm - lambda) / twoNorm;

			/*
			* shrinkage this group by ratio
			*/
			for (j = (int) ind[3 * i] - 1; j < (int) ind[3 * i + 1]; j++)
				v[j] *= ratio;
		}
		else{
			/*
			* threshold this group to zero
			*/
			for (j = (int) ind[3 * i] - 1; j < (int) ind[3 * i + 1]; j++)
				v[j] = 0;
		}
	}

	for (int i = 0; i < n; i++)
		if (v[i] != 0)
			return false;
	return true;
}

/*sub function
  - solving fused proximal operator
  This function is written by Dr.Laurent Condat

  Reference:
  L. Condat, A direct algorithm for 1D total variation denoising, preprint hal-00675043, 2011.

  Note:
  lambda is the fused regualrization paramter, and mu is the l1 regualrization parameter!!
*/
void fused_lasso(double* output, double* input, const size_t width, const double lambda, const double mu) {
	if (width > 0) {
		/*to avoid invalid memory access to input[0]*/
		int k = 0, k0 = 0;			/*k: current sample location, k0: beginning of current segment*/
		double umin = lambda, umax = -lambda;	/*u is the dual variable*/
		double vmin = input[0] - lambda, vmax = input[0] + lambda;	/*bounds for the segment's value*/
		int kplus = 0, kminus = 0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const double twolambda = 2.0*lambda;	/*auxiliary variable*/
		const double minlambda = -lambda;		/*auxiliary variable*/
		for (;;) {
			/*simple loop, the exit test is inside*/
			while (k == width - 1) {
				/*we use the right boundary condition*/
				if (umin<0.0) {
					/*vmin is too high -> negative jump necessary*/
					vmin = vmin > mu ? vmin - mu : vmin < -mu ? vmin + mu : 0.0;
					do output[k0++] = vmin; while (k0 <= kminus);
					umax = (vmin = input[kminus = k = k0]) + (umin = lambda) - vmax;
				}
				else if (umax > 0.0) {
					/*vmax is too low -> positive jump necessary*/
					vmax = vmax > mu ? vmax - mu : vmax < -mu ? vmax + mu : 0.0;
					do output[k0++] = vmax; while (k0 <= kplus);
					umin = (vmax = input[kplus = k = k0]) + (umax = minlambda) - vmin;
				}
				else {
					vmin += umin / (k - k0 + 1);
					vmin = vmin > mu ? vmin - mu : vmin < -mu ? vmin + mu : 0.0;
					do output[k0++] = vmin; while (k0 <= k);
					return;
				}
			}

			if (k > width - 1)
				break;

			if ((umin += input[k + 1] - vmin)<minlambda) {
				/*negative jump necessary*/
				vmin = vmin > mu ? vmin - mu : vmin < -mu ? vmin + mu : 0.0;
				do output[k0++] = vmin; while (k0 <= kminus);
				vmax = (vmin = input[kplus = kminus = k = k0]) + twolambda;
				umin = lambda; umax = minlambda;
			}
			else if ((umax += input[k + 1] - vmax) > lambda) {
				/*positive jump necessary*/
				vmax = vmax > mu ? vmax - mu : vmax < -mu ? vmax + mu : 0.0;
				do output[k0++] = vmax; while (k0 <= kplus);
				vmin = (vmax = input[kplus = kminus = k = k0]) - twolambda;
				umin = lambda; umax = minlambda;
			}
			else {
				/*no jump necessary, we continue*/
				k++;
				if (umin >= lambda) {
					/*update of vmin*/
					vmin += (umin - lambda) / ((kminus = k) - k0 + 1);
					umin = lambda;
				}
				if (umax <= minlambda) {
					/*update of vmax*/
					vmax += (umax + lambda) / ((kplus = k) - k0 + 1);
					umax = minlambda;
				}
			}
		}
	}
}

/*sub function
  - solving group proximal operator
*/
void group_lasso(double *out, double *in, int K, double lam, double rho){
	double temp = 0;
	for (int k = 0; k < K; k++){
		if (in[k] < -lam)
			in[k] += lam;
		else if (in[k] > lam)
			in[k] -= lam;
		else
			in[k] = 0;
		temp += in[k] * in[k];
	}
	temp = sqrt(temp);
	temp = temp > rho ? (1 - rho / temp) : 0;
	for (int k = 0; k < K; k++)
		out[k] = temp * in[k];
}

/*sub function
  - solving tree proximal operator
*/
void tree_lasso(double *out, double *in, int K, double lam, double rho){
	/*To be implemented.........*/
	/* use altra(double *x, double *v, int n, double *ind, int nodes)*/
}

/*
 Find free set for shrinking techique
 penalty:
   0: fused
   1: group
   2: tree
*/
size_t findFreeSet(double *G, double *P, double *v, double *fx, double *fy, const int N, const int K,
	const double lam, const double rho, const double EPS, int penalty, const double* ind = NULL, const int nodes = 0){
	size_t indxNum = 0;
	bool flag = true;
	double *ind_work = NULL;
	size_t indx, indxi, N2 = N*N;
	
	if (penalty == 2 && ind){
		ind_work = new double[3 * nodes];
		memcpy(ind_work, ind, 3 * nodes*sizeof(double) );
		for (int i = 0; i < nodes; i++)
			ind_work[3 * i + 2] *= rho;
	}

	for (int i = 0; i < N; i++){
		for (int j = i + 1; j < N; j++){
			indx = i*N + j;
			flag = true;
			for (int k = 0; k < K; k++){
				indxi = indx + k*N2;
				if (P[indxi] != 0.0){
					flag = false;
					break;	
				}else
			    	v[k] = G[indxi];
			}

			if (flag)
			{
				switch (penalty){
				case 0:
					flag = screeningFused(v, K, lam, rho, EPS);
					break;
				case 1:
					flag = screeningGroup(v, K, lam, rho, EPS);
					break;
				case 2:
					flag = screeningTree(v, K, ind_work, nodes, 0.0);
					break;
				default:
					flag = screeningFused(v, K, lam, rho, EPS);
					break;
				}
			}

			if (flag)
				continue;

			fx[indxNum] = (double) i;
			fy[indxNum] = (double) j;
			indxNum++;
		}
	}

	if (ind_work)
		delete [] ind_work;
    return indxNum;
}

/*
  sub function
  - solving proximal operator of multiple grpahical lasso
*/
void mgl_prox(double *output, double *input, const double *fx, const double *fy, double lam, double rho,
	const int N, const int K, const size_t Num, const int penalty,const double* ind = NULL, const int nodes=0){
	size_t N2 = N * N;
	double *sarrayi = new double[K];
	double *sarrayo = new double[K];
	double *ind_work = NULL;
	size_t idx1, idx2;
    
	if (penalty == 2 && ind){
		ind_work = new double[3 * nodes];
		memcpy(ind_work, ind, 3 * nodes*sizeof(double) );
		for (int i = 0; i < nodes; i++)
			ind_work[3 * i + 2] *= rho;
	}

	for (size_t i = 0; i < Num; i++){
		double temp = fx[i] * N + fy[i];
		idx1 = (size_t) temp;
		temp = fy[i] * N + fx[i];
		idx2 = (size_t) temp;
		for (int k = 0; k < K; k++)
			sarrayi[k] = input[idx1 + k*N2];
		switch (penalty){
		case 0:
			fused_lasso(sarrayo, sarrayi, K, rho, lam);
			break;
		case 1:
			group_lasso(sarrayo, sarrayi, K, lam, rho);
			break;
		case 2:
			altra(sarrayo, sarrayi, K, ind_work, nodes);
			break;
		default:
			fused_lasso(sarrayo, sarrayi, K, rho, lam);
			break;
		}
		
		for (int k = 0; k < K; k++)
			output[idx1 + k*N2] = output[idx2 + k*N2] = sarrayo[k];
	}

	for (size_t k = 0; k < K; k++){
		size_t kN2 = k*N2;
		for (size_t i = 0; i < N; i++)
			output[i*N + i + kN2] = input[i*N + i + kN2];
	}

	delete [] sarrayi;
	delete [] sarrayo;
	if (ind_work)
		delete [] ind_work;
}

/*
  global screening for multiple graphical lasso
*/
size_t globalScreening(double *G, double *adj, const int N, const int K, const double lam, const double rho, const int penalty,
	                   const double* ind = NULL, const int nodes = 0){
	size_t indxNum = 0;
	size_t indx, N2 = N * N;
	double *ind_work = NULL;

	if (penalty == 2 && ind){
		ind_work = new double[3 * nodes];
		memcpy(ind_work, ind, 3 * nodes*sizeof(double) );
		for (int i = 0; i < nodes; i++)
			ind_work[3 * i + 2] *= rho;
	}

	double *v = new double[K];
	for (int i = 0; i < N; i++){
		for (int j = i + 1; j < N; j++){
			indx = i * N + j;
			for (int k = 0; k < K; k++)
				v[k] = G[indx + k * N2];
			bool flag = false;
			switch (penalty){
			case 0:
				flag = screeningFused(v, K, lam, rho, 0.0);
				break;
			case 1:
				flag = screeningGroup(v, K, lam, rho, 0.0);
				break;
			case 2:
				flag = screeningTree(v, K, ind_work, nodes, 0.0);
				break;
			default:
				flag = screeningFused(v, K, lam, rho, 0.0);
				break;
			}

			if (flag)
				continue;
			adj[i*N + j] = adj[N*j + i] = 1;
			indxNum++;
		}
	}
	delete [] v;
	if (ind_work)
		delete ind_work;
	return indxNum;
}

double treeNormM(double *G, int N, int K, double *ind, int nodes){
	double tN = 0.0;
	double tt = 0.0;
	size_t indx, N2 = N * N;
	double *v = new double[K];
	for (int i = 0; i < N; i++)
		for (int j = i + 1; j < N; j++){
			indx = i * N + j;
			for (int k = 0; k < K; k++)
				v[k] = G[indx + k * N2];
			
			tN += treeNormR(v, K, ind, nodes);
		}
		delete [] v;
	return tN * 2;
}