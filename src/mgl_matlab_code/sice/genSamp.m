function [ X ] = genSamp( C, m )
%GENSAMP Generate samples according to a sparse inverse covariance matrix

% Input: C      -- n x n sparse inverse covariance matrix
%        m      -- number of samples
% Ouput: X      -- n x m matrix of samples, each column is sampled i.i.d
%                  from N(0, inv(C))

% Author:    Richard Y. Zhang (ryz@alum.mit.edu)
% Url:       http://alum.mit.edu/www/ryz
% Date:      Feb 2018
% Reference: R.Y. Zhang, S. Fattahi, S. Sojoudi, "Linear-Time Algorithm for 
%            Learning Large-Scale Sparse Graphical Models".

n = size(C,1);
assert(size(C,2) == n, 'Inverse covariance must be square');
assert(norm(C - C' ,'fro') == 0, 'Inverse covariance must be perfectly symmetric');
assert(isreal(C), 'Inverse covariance must be real');

% Attempt to factor
p = amd(C);
[L, fail] = chol(C(p,p),'lower');
assert(fail==0, 'Inverse covariance must be posdef');

% Generate samples
X = randn(n,m); % "stochastic gem"
X = L' \ X;
X(p,:) = X;

end

