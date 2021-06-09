function [ X , S] = fastmdmc(C)
%FASTMDMC Use Newton-CG to solve the maximum determinant matrix 
%completion problem
%  minimize tr CX - logdet X s.t. X_{i,j}=0 if C_{i,j}=0
% via its dual
%  maximize logdet Z s.t. Z_{i,j} = C_{i,j} if C_{i,j}=0
% The maxdet matrix completion is recovered as Z = inv(X). 

% Output: primal-dual solutions X and S satisfying 
%    Xij = 0           for all Cij  = 0 (primal feasibility)
%    Sij = Cij         for all Cij ~= 0 (dual feasibility)
%    Sij = [inv(X)]_ij for all Sij ~= 0 (complementary slackness)

% Author:    Richard Y. Zhang (ryz@alum.mit.edu)
% Url:       http://alum.mit.edu/www/ryz
% Date:      Feb 2018
% Reference: R.Y. Zhang, S. Fattahi, S. Sojoudi, "Linear-Time Algorithm for 
%            Learning Large-Scale Sparse Graphical Models". 

tt = tic;

% Configuration
bt_alpha      = 0.01;  % BT linesearch minimum decrement
bt_beta       = 0.5;   % BT linesearch shrinkage
ipm_tol       = 1e-12; % IPM exit threshold (Newton decrement)
max_ipm_iter  = 100;   % Max number of outer IPM iterations
cg_tol        = 1e-12;  % CG exit threshold (suboptimality)
max_cg_iter   = 100;   % Max number of inner CG iterations
debug         = false; % debug mode is slower but performs checks

% Symmetricize C and obtain data
n = size(C,1);
assert(size(C,2) == n, 'C must be square');
assert(issparse(C), 'C must be sparse');
C = (C+C')/2;

% Chordal completion for the sparsity pattern of C.
% To obtain the explicit chordal embedding, run 
%   E = Lsym | Lsym'; E(p,p) = E; 
% Here, p is the associated perfect elimination ordering.
p = amd(C); C = C(p,p);
[~, ~, par, post, Lsym] = symbfact(C,'sym','lower'); 

% Preapply postordering to avoid random access
p = p(post); 
Lsym = Lsym(post,post); C = C(post,post); 
par = reorderTree(par, post);
ch = p2ch(par); % children of each node on the elimination tree

% Data structure for compressed column storage
[row_ind, row_cnt] = compressCol(Lsym); % row indices for every column
mask = getmask(row_ind, par); 

% Sparsity pattern
fill = find(Lsym); orig = find(tril(C));
add_mask = ~ismember(fill,orig);
add = fill(add_mask); % newly added edges
m = numel(add);

% Begin loop
y = zeros(m,1);
total_cg = 0; total_time = 0;
if m > 0 % Closed-form solution?
    flag = 0;
    fprintf('     | Potential  |           |           |     |           |          |\n');
    fprintf(' it# | function   | Gradient  | Newt Decr | PCG | Step Size | Time (s) |\n');
    fprintf('-------------------------------------------------------------------------\n');
    for iter = 0:max_ipm_iter
        % Evaluate objective function and gradient
        % R is used in the Hessian
        [Lccs,D,fail,R] = matComp(C - Aop(y));        
        if fail > 0, error('y=0 is not a valid starting point'); end
        [f,g] = obj_fun(Lccs,D);
        
        % Display information
        this_time = toc(tt); tt = tic;
        total_time = total_time + this_time;
        if iter == 0
            fprintf(' %3d | %+.3e | %.3e |                             | %.3e \n',...
            iter,f,norm(g),this_time);
        else
            fprintf(' %3d | %+.3e | %.3e | %.3e | %3d | %.3e | %.3e \n',...
            iter,f,norm(g),newt_decr, numel(resvec), t, this_time);
            total_cg = total_cg + numel(resvec);
            if flag == 1, break; end
        end
        
        
        % Compute Newton direction
        [dy,~,~,~,resvec] = my_cg(@(y)mvSchur(Lccs,D,R,y), g, cg_tol, max_cg_iter);
        semilogy(resvec); pause(0.01);
        % Gradient descent if CG does not provide a useful direction
        if all(dy==0), dy = g; end
        
        % Backtracking Line search
        % See pg. 464 of Boyd & Vandenberghe
        newt_decr = full(dy'*g); t = 1;
        if newt_decr > ipm_tol
            while obj_fun(y - t*dy) > f - bt_alpha*t*newt_decr
                t = bt_beta*t;
            end
        else % get ready to step out
            flag = 1;
        end
        y = y - t*dy;
    end
end

% Post process the output
S = C - Aop(y); 
[L,D] = matComp(S);
X = chol2mat(L,D);
fprintf('\n\n');
fprintf(' Optimality:   '); G = C - projInv(L,D);
fprintf(' |P(C - inv(X))|/|C| = %.3e \n', norm(G(orig))/norm(C(:)));
fprintf(' Feasibility:  ');
fprintf(' |P(X)-X|/|X|        = %.3e \n', norm(X(add))/norm(X(:)));
fprintf('\n');
fprintf(' Total cg iters  = %d \n', total_cg);
fprintf(' Total time      = %.3e sec \n', total_time);

% Rearrange
S(p,p) = S; X(p,p) = X;



%==========================================================================
% Utility functions (nested)
%==========================================================================
function [val] = sparse2ccs(Matrix)
% convert a sparse matrix into the compressed column storage format
    if debug
        assert(isempty(setdiff(find(Matrix), fill)),...
            'Input matrix does not have the required sparsity pattern');
    end
    val = cell(1,n);
    for j = 1:n
        val{j} = full(Matrix(row_ind{j},j));
    end
end

function [Matrix] = ccs2sparse(val)
    % Convert CCS back to sparse matrix
    col_ind = cell(1,n);
    for j = 1:n
        col_ind{j} = j*ones(size(row_ind{j}));
    end
    ii = cat(1,row_ind{:});
    jj = cat(1,col_ind{:});
    kk = cat(1,val{:});
    Matrix = sparse(ii,jj,kk,n,n,numel(kk));
end
%==========================================================================
% IPM functions (nested)
%==========================================================================
function [f,g] = obj_fun(L,D)
% Evaluate phi_*(C-Aop(y))
    if nargin < 2 || isempty(D)
    [L,D,fail_] = matComp(C - Aop(L));
        if fail_>0
            f = inf; 
            g = nan(m,1); 
            return; 
        end
    end
    f = sum(log(D));
    if nargout > 1
        g = Aadjop(chol2mat(L,D));
    end
end
function Y = Aop(y)
% Convert the degrees of freedom into the matrix
    Y = sparse(add,1,y,n^2,1);
    Y = reshape(Y,n,n);
    Y = Y + tril(Y,-1)';
end
function y = Aadjop(Y)
    y = Y(add);
end
function [L, D, fail, R] = matComp(S)
%MATCOMP Compute X satisfying S = P(inv(X)), where S must have 
% chordal sparsity pattern Lsym

% Compute tree decomposition
S = sparse2ccs(tril(S)); % cliques
fail = 0;

% Andersen, Dahl, Vandenberghe -- Algorithm 4.2
L = cell(1,n);
D = zeros(1,n);
V = cell(1,n); % Update matrices
R = cell(1,n); % Update matrices
for j = n:-1:1 % reverse postorder = preorder    
    % j-th column of S
    % NB: S is lower-triangular, so Sj(1) = j by construction.
    Sj = S{j}; 
    Sjj = Sj(1); SIj = Sj(2:end);
    Jj = row_ind{j}; 

    % Run formula
    if numel(Jj) > 1 % Not a root node
        LIj = -R{j}'\(R{j}\SIj);
        invDj = (Sjj + SIj'*LIj);
    else
        invDj = Sjj;
        LIj = [];
    end
    
    % Failure check
    if invDj > 0
        D(j) = 1/invDj;
    else
        fail = j; break;
    end

    % Store
    L{j} = [1; LIj];

    % Propagate to children
    tmp = [Sjj, SIj'; SIj, V{j}];
    if nargout < 4 
        % release update matrices if no longer needed
        V{j} = []; R{j} = []; 
    end
    for i = ch{j}
        V{i} = tmp(mask{i},mask{i});
        R{i} = chol(V{i},'lower');
    end
end
if nargout < 2 
    % In case of one output, we will form the explicit matrix
    L = chol2mat(L,D);
end
end
function [S, fail] = projInv(L, D)
%PROJINV Compute P(inv(X)), where X must have chordal sparsity pattern Lsym
if nargin < 2 || isempty(D) 
    % In case of one argument, we will have to do the factoring
    [L, D, fail] = mat2chol(L);
    if nargout < 2 && fail>0
        error('projInv only works with pos def matrices');
    end
end

% Andersen, Dahl, Vandenberghe -- Algorithm 4.1
S = cell(1,n); % off-diagonal elements of S
V = cell(1,n); % Update matrices
for j = n:-1:1 % reverse postorder = preorder    
    % j-th column of L
    % NB: L is lower-triangular, so Jj(1) = j by construction.
    Lj = L{j};

    % Run formula
    if numel(Lj) > 1 % Not a root node
        LIj = Lj(2:end);
        SIj = -V{j}*LIj;
        Sjj = 1/D(j) - SIj'*LIj;
    else
        Sjj = 1/D(j);
        SIj = [];
    end

    % Store
    S{j} = [Sjj; SIj];

    % Propagate to children
    tmp = [Sjj, SIj'; SIj, V{j}];
    V{j} = []; % Free this update matrix
    for i = ch{j}
        V{i} = tmp(mask{i},mask{i});
    end
end

% Form symmetric matrix from lower-triangular part
S = ccs2sparse(S);
S = S + tril(S,-1)';
end
function X = chol2mat(L,D)
    L = ccs2sparse(L);
    X = L*diag(sparse(D)/2)*L';
    X = X + X';
end
function [L,D,fail] = mat2chol(X)
    [L, fail] = chol(X,'lower');
    if fail == 0
        D = diag(L); L = L/diag(D); D = D.^2;
    end
    L = sparse2ccs(L);
end
%==========================================================================
% Custom fast matrix-vector product with the Schur complement matrix
%       data = factorSchur(L,D,S); z = mvSchur(data,y);
% is equivalent to
%       z = Aadjop(invHess(L,D,S,Aop(y)))
%==========================================================================
function z = mvSchur(L,D,R,y)
% Put matrices into CSS form
% Faster version of Yccs = sparse2ccs(Aop(y))    
T = zeros(numel(fill),1);
T(add_mask) = y;
T = mat2cell(T,row_cnt);

% Andersen, Dahl, Vandenberghe -- Algorithm 5.2
% Step 1 and 2 combined
K = cell(1,n); 
Vp = cell(1,n); % Update matrices
for j = n:-1:1 % reverse postorder = preorder
    
    % j-th column of M
    Tj = T{j}; Tjj = Tj(1); TIj = Tj(2:end);
    Dj = D(j); LIj = L{j}(2:end); % LIj = Lj(2:end)
    Rj = R{j}; % Rj * Rj' = Vj = S_{IjIj}.

    % Run formula
    % Here, we explicitly expand the Lslice'*T_front*Lslice formula
    nn = numel(TIj);
    if nn > 0 % Not a root node
        MIj = TIj + Vp{j}*LIj;
        Mjj = Tjj + LIj'*(MIj+TIj);
        Kjj = Dj^2*Mjj;
        KIj = Dj*(Rj'\(Rj\MIj));
    else
        Kjj = Dj^2*Tjj;
        KIj = [];
    end
    K{j} = [Kjj; KIj];

    % Propagate to children
    T_front = [Tjj, TIj'; TIj, Vp{j}];
    Vp{j} = []; % Free this update matrix
    for i = ch{j}
        Vp{i} = T_front(mask{i},mask{i});
    end
end

% Step 3
Y = cell(1,n); 
Up = cell(1,n); % Update matrices
for j = 1:n % postorder
    
    % j-th column of K
    Kj = K{j}; Lj = L{j};
    Kjj = Kj(1); KIj = Kj(2:end);

    % Run formula
    nj = numel(KIj);
    if nj > 0
        % We explicitly expand the Lslice*K_front*Lslice' formula
        Y_front = Lj*[0; KIj]'; 
        Y_front = Y_front + Y_front'; 
        Y_front = Y_front + Kjj*(Lj*Lj');
    else
        Y_front = Kjj;
    end
    for i = ch{j}
        Y_front(mask{i},mask{i}) = Y_front(mask{i},mask{i}) - Up{i};
        Up{i} = []; % Free the update matrix
    end
    Y{j} = Y_front(:,1);
    Up{j} = -Y_front(2:end,2:end);
end

% Form symmetric matrix from lower-triangular part
% Faster version of z = Aadjop(sparse2ccs(Y))
Y = cat(1,Y{:});
z = Y(add_mask);
end
end

%==========================================================================
% Utility functions
%==========================================================================
function [x,flag,f,ii,resvec,r] = my_cg(A,b,tol,maxit)
%MY_CG Conjugate Gradients Method.
%      Maximize f(x) = x'*A*x - 2*x'*b with x_0 = [0,...,0]'
%      Exit when f(x_{k+1}) - f(x_{k}) < tol * f(x_{k+1})
%      Modifications:
%      -- Removed stagnation detection, which gives up too early around the
%         numerical floor.
%      -- Always output the final iterate.
%      -- "iter" is the total number of iterations, not the best.
%      -- Removed support for matrix input, missing arguments
%      (Simplified and significantly modified from the MATLAB internal code)
n = numel(b);

% Check for all zero right hand side vector => all zero solution
if (norm(b) == 0)                  % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    f = 0;                         % the relative residual is actually 0/0
    ii = 0;                        % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    r = zeros(n,1);
    return
end

% Set up for the method
x = zeros(n,1);
f = 0;
flag = 1;
r = b;
resvec = zeros(maxit,1);         % Preallocate vector for norm of residuals
rho = 1;

% loop over maxit iterations (unless convergence or failure)
for ii = 1 : maxit
    rho1 = rho;
    rho = r' * r;
    if ((rho == 0) || isinf(rho))
        flag = 4;
        break
    end
    if (ii == 1)
        p = r;
    else
        beta = rho / rho1;
        if ((beta == 0) || isinf(beta))
            flag = 4;
            break
        end
        p = r + beta * p;
    end
    q = A(p);
    pq = p' * q;
    if ((pq <= 0) || isinf(pq))
        flag = 4;
        break
    else
        alpha = rho / pq;
    end
    if isinf(alpha)
        flag = 4;
        break
    end    
    x = x + alpha * p;             % form new iterate
    r = r - alpha * q;
    df  = 2*alpha*x'*q + alpha^2*pq - 2*alpha*b'*p;
    f = f + df;
    resvec(ii) = df/f;
    
    % check for convergence
    if resvec(ii) <= tol
        flag = 0;
        break
    else
        
    end
end                                % for ii = 1 : maxit

% returned solution is first with minimal residual
resvec = resvec(1:ii,:);

end

function [row_ind, row_cnt] = compressCol(Pattern)
% Generate definition of CCS format
row_ind = cell(1,size(Pattern,2));
row_cnt = zeros(1,size(Pattern,2));
for j = 1:size(Pattern,2)
    row_ind{j} = find(Pattern(:,j));
    row_cnt(j) = numel(row_ind{j});
end
end

function [ch, root] = p2ch(p)
% Convert the list of parent vectors into a cell array of children vectors.
n = length(p);
ch = cell(size(p));
root = [];
for i = 1:n
    par = p(i);
    if par > 0
        % Push child into parent
        ch{par} = [ch{par}, i];
    elseif par == 0
        root = [root, i];
    end
end
end

function mask = getmask(row_ind, par)
% Precompute the masks such that 
%   M(mask{i}, mask{i}) = E_{JjIi}' * M * E_{JjIi} where j = par(i)
% and
%   E_{IJ}(i,j) = 1 if I(i) = J(j) and 0 otherwise,
%   Jj = row_ind{j}
%   Ii = row_ind{i}(2:end)
% 
n = numel(row_ind);
mask = cell(1,n);
for i = 1:n
    j = par(i);
    if j>0
        Jj = row_ind{j}; Ji = row_ind{i};
        mask{i} = ismember(Jj, Ji(2:end));
    end
end
end

function [ parent2 ] = reorderTree( parent, v )
%RELABELTREE Reorder a given tree vector
% Let v : order -> vertex. Then define parent2 such that
%   if j = parent2(i), then v(j) = parent( v(i) )
%      0 = parent2(i), then   0  = parent( v(i) )

n = numel(parent); v = v(:)';
assert(all(sort(v) == 1:n), 'provided v is not an ordering!');
vinv(v) = 1:n; % inverse v : vertex -> order
parent2 = zeros(1,n);
for i = 1:n
    vj = parent(v(i)); 
    if vj > 0 % not a root
        parent2(i) = vinv(vj);
    end
end
end

