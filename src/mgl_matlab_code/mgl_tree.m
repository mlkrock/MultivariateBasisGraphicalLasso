function [P,funVal] = mgl_tree(S,params)

global marks

if nargin < 2
  error('params.lambda and params.rho should be specified!');
end

if isfield(params,'lambda')
    lambda = params.lambda;
else
   error('params.lambda should be specified!');
end

if ~isfield(params,'maxlineiter')
    params.maxlineiter = 50;
end
if ~isfield(params,'SPGreltol')
    params.SPGreltol = 1e-6;
end
if ~isfield(params,'SPGmaxiter')
    params.SPGmaxiter = 50;
end
if ~isfield(params,'Newtontol')
   params.Newtontol = 1e-6;
end

if ~isfield(params,'maxiter')
   params.maxiter = 1500;
end
if ~isfield(params,'sigma')
   params.sigma = 0.001;
end

[n,n,K] = size(S);
if isfield(params, 'P0')
    P = params.P0;
else
    P = zeros(n,n,K);
    for k = 1:K
        P(:,:,k) = diag(1./diag(S(:,:,k)));
    end
end

if(~(isfield(params, 'ind')))
    error('For tree, params.ind should be specified!');
end

marks = zeros(n,n,K);
for k = 1:K
    marks(:,:,k) = 1 - eye(n,n);
end

funVal = [];
R = zeros(n,n,K);
B = R;
for k = 1:K
    [R(:,:,k),info] = chol(P(:,:,k));
    if info > 0
        error('%d not positive definite', k);
    end
    invR = inv(R(:,:,k));
    B(:,:,k) = invR*invR';
end

[fobjp, l1normX] = computLogDet(P, S, R, lambda, params.ind);
funVal(1) = fobjp;
fx = zeros((n-1)*n/2,1);
fy = zeros((n-1)*n/2,1);
tev = zeros(K,1);
iter = 1;
while (1)
    iter = iter+1;
    G = S - B;
    numF = findFreeSet(G,P,fx,fy,0,lambda,tev,2,params.ind,size(params.ind,2));
    params.NumF = numF;
    D = mgl_spg(S,P*1.0,params,B,fx,fy,max(params.SPGreltol,0));

    if (norm(D(:),'inf') < 5e-5)
        break;
    end
    
   %linesearch
    alpha = 1;
    PW = (P+D).*marks;
    l1normXD = lambda*treeNormM(PW,params.ind,size(params.ind,2));
    trdg = 0;
    
    for k = 1:K
        t1 = G(:,:,k);
        t2 = D(:,:,k);
        trdg = trdg + t1(:)'*t2(:);
    end
    
    fobjp1 = 1e15;
    for liter = 1 : params.maxlineiter
        W = P+alpha*D;
        flag = 0;
        for k = 1:K
            [Tz,info] = chol(W(:,:,k));
            if info > 0
                flag = 1;
                break;
            end
            R(:,:,k) = Tz;
        end
        
        if flag == 1
            alpha = alpha/2;
            continue;
        end
        
        [fobj, l1normX1] = computLogDet(W, S, R, lambda, params.ind);
      % fprintf('Linear search decreases %f and norm of D is %f, number of free set %d\n', fobjp - fobj, norm(D(:),'inf'), numF);
        fprintf('outer:%d inner%d, Linear search decreases %f and norm of D is %f, number of free set %d\n', iter, liter,fobjp - fobj, norm(D(:),'inf'), numF);
        if fobj <= fobjp + alpha*params.sigma*(trdg + l1normXD - l1normX);
            l1normX = l1normX1;
            break;
        end
        
        if fobjp1 < fobj
            l1normX = l1normX1;
            break
        end
        fobjp1 = fobj;
        alpha = alpha/2;
    end
    P = W;
    
   % fprintf('fset %f, fist %f,  linet %f\n', fset,fist,linet);
    
    funVal(iter) = fobj;
    
    fobjp = fobj;
    if iter == 1
        for k = 1:K
            invR = inv(R(:,:,k));
            B(:,:,k) = invR*invR';
        end
        continue
    end
    
    if(abs(funVal(iter) - funVal(iter-1)) < params.Newtontol*abs(funVal(iter-1))) || norm(D(:),'inf') < 5e-5
        break;
    end
    
    if iter>params.maxiter
        break;
    end
    
    for k = 1:K
        invR = inv(R(:,:,k));
        B(:,:,k) = invR*invR';
    end
end

end

function [fobj, l1normX] = computLogDet(W, S, R, lambda, ind)
global marks
[n,n,K] = size(S);
logdet = 0;
trdg = 0;
for k = 1:K
    logdet = logdet + 2*sum(log(diag(R(:,:,k))));
    t1 = S(:,:,k);
    t2 = W(:,:,k);
    trdg = trdg + t1(:)'*t2(:);
end
W = W.*marks;
l1normX = lambda*treeNormM(W,ind,size(ind,2));

fobj = trdg - logdet + l1normX;
end


function opterr = objcompute(D,hatP, B, S, lambda, ind)
global marks
hatP = hatP.*marks;
opterr = 0;
[n,n,K] = size(hatP);
for k = 1:K
    Wz = B(:,:,k)*D(:,:,k);
    opterr = opterr + Wz(:)'*Wz(:)/2 + ...
        trace((S(:,:,k) - B(:,:,k))*D(:,:,k));
end
opterr = opterr + lambda*treeNormM(hatP,ind,size(ind,2));
end

function [W,funVal,lambdak] = mgl_spg(S,P,params,B,fx,fy,termtol)
lambdamin = 1e-30;
lambdamax = 1e5;
gamma = 1e-4;
eta = 2;
M = 10;
f_1 = -inf*ones(M,1);

K = size(S,3);
l1reg = params.lambda;

SB = S - 2*B;
grad = SB;
hatP = 1.0*P; P_start = 1.0*P;
iter = 1;
for k = 1:K
    grad(:,:,k) = B(:,:,k)*hatP(:,:,k)*B(:,:,k);
end
grad  = grad + SB;
z = hatP - grad;
mgl_prox(P,z,fx,fy,0, l1reg, params.NumF,2,params.ind,size(params.ind,2));
%%
f_1(iter) = objcompute(hatP-P_start,hatP, B, S, l1reg,params.ind);
funVal(iter) = f_1(iter);
pomga = P - hatP;
err = norm(pomga(:),'inf');
if err <= termtol
    W = P - P_start;
    return;
end

if isfield(params,'SPGlambdak')
    lambdak =  params.SPGlambdak;
else
    lambdak = 1/(err + 1e-15);
end
P = 1.0*hatP;
iter = iter+1;

while (iter < params.SPGmaxiter)
    
    cont_inner = 1;
    while cont_inner
        z = hatP - 1/lambdak*grad;
        mgl_prox(P,z,fx,fy,0,l1reg/lambdak, params.NumF,2,params.ind,size(params.ind,2));
        D = P - hatP;
        fnew = objcompute(P-P_start,P, B, S, l1reg,params.ind);
        ddr = D(:)'*D(:);
        f_threshold = max(f_1) - 0.5*gamma*lambdak*ddr;
        if fnew <= f_threshold
            cont_inner=0;
        else
            lambdak = eta*lambdak;
            %fprintf(1,' f=%10.6e, increasing lambdak to %6.2e\n', fnew, lambdak);
        end
    end
    f_1(mod(iter,M)+1) = fnew;
	iter = iter+1;
    funVal(iter) = fnew;
    hatgrad = grad;
    for k = 1:K
        grad(:,:,k) = B(:,:,k)*P(:,:,k)*B(:,:,k);
    end
    grad  = grad + SB;
    gk = grad - hatgrad;
    td = D(:)'*gk(:);
    if td <=0
        lambdak = lambdamax;
        continue
    end    
    td = td/ddr;
    lambdak = max(lambdamin,min(td,lambdamax));
    
    
    pomga = P - hatP;
    err = norm(pomga(:),'inf')/norm(P(:),'inf');
    
    if err <= termtol
        break;
    end
    
    hatP = 1.0*P;   
end
W = P - P_start;
end    

