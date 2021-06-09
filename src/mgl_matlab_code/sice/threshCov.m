function S = threshCov(X, tau)
%THRESHCOV Given the zero-mean sample matrix X, soft-threshold the 
% associated sample covariance matrix M = 1/m * X*X' by the threshold tau.
%       For each i,j = 1...n set
%          Sij = Mij         if i = j
%                Mij - tau   if Mij > tau
%                0           if |Mij| < tau
%                Mij + tau   if Mij < -tau
% To recover the graphical lasso estimator, run X = fastmdmc(S);

% This routine partitions the matrix M into 4000 x 4000 subblocks, and
% thresholds using dense linear algebra.

% Author:    Richard Y. Zhang (ryz@alum.mit.edu)
% Url:       http://alum.mit.edu/www/ryz
% Date:      Feb 2018
% Reference: R.Y. Zhang, S. Fattahi, S. Sojoudi, "Linear-Time Algorithm for 
%            Learning Large-Scale Sparse Graphical Models".

% If we allocate more than this number times n nonzeros, terminate
% prematurely.
MEMORY_LIMIT = 4000; % x n

% Get partition size
[n, m] = size(X);
skip = ceil(sqrt(n*MEMORY_LIMIT));
X = X/sqrt(m);
xd = sum(X.^2,2);
if n <= 10*MEMORY_LIMIT 
    % If problem size small, just do it explicitly
    M = triu(X*X');
    [ii_pos, jj_pos, kk_pos] = find(max(M-tau,0));
    [ii_neg, jj_neg, kk_neg] = find(min(M+tau,0));
    S = sparse([ii_pos; ii_neg],[jj_pos; jj_neg],[kk_pos; kk_neg],...
               n,n);
    S = S + triu(S,1)';
   return
end

% Partition x into cells
tmp = floor(n/skip);
rows = [repmat(skip,tmp,1); n - tmp*skip];
cumrows = [0;cumsum(rows)];
X = mat2cell(X, rows);
p = numel(X);

% Threshold one submatrix of M = x*x' at a time.
[ii,jj,kk] = deal(cell(1,p)); 
total_nz = 0;
for j = 1:p
    [this_ii, this_jj, this_kk] = deal(cell(1,j));
    for i = 1:j
        % Form matrix
        Mij = X{i}*X{j}';
        if i == j
            Mij = triu(Mij,1);
        end
        % Do soft-thresholding
        [ii_pos, jj_pos, kk_pos] = find(max(Mij-tau,0));
        [ii_neg, jj_neg, kk_neg] = find(min(Mij+tau,0));
        % Record nonzeros
        this_ii{i} = [ii_pos; ii_neg] + cumrows(i);
        this_jj{i} = [jj_pos; jj_neg] + cumrows(j);
        this_kk{i} = [kk_pos; kk_neg];
        % Sum nonzeros
        total_nz = total_nz + numel(this_ii{i});
        % Check for memory issues
        if total_nz > MEMORY_LIMIT * n
            error('REACHED MEMORY LIMIT. EXITING....');
        end
    end
    % Assemble this column
    ii{j} = cat(1, this_ii{:});
    jj{j} = cat(1, this_jj{:});
    kk{j} = cat(1, this_kk{:});
end
% Assemble all columns
ii = cat(1, ii{:});
jj = cat(1, jj{:});
kk = cat(1, kk{:});

% form sparse matrix
S = sparse(ii,jj,kk,n,n);
S = S + S';
S = S + diag(sparse(xd));
end