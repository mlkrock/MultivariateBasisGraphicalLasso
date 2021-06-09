function [P,funVal,iterTime, subOpt] = mgl_adp(S,params)

global marks FX FY

if nargin < 2
    error('params.lambda and params.rho should be specified!');
end

% set up the params
params = setParms(params);

if isfield(params,'lambda')
    lambda = params.lambda;
else
    error('params.lambda should be specified!');
end

if isfield(params,'rho')
    rho = params.rho;
else
    error('params.rho should be specified!');
end


% 
[n,n,K] = size(S);
if isfield(params, 'P0')
    P = params.P0;
else
    P = zeros(n,n,K);
    for k = 1:K
        P(:,:,k) = diag(1./diag(S(:,:,k)));
    end
end

% set up global varaibles 
marks = zeros(n,n,K);
for k = 1:K
    marks(:,:,k) = 1 - eye(n,n);
end

FX = zeros((n-1)*n/2,1);
FY = zeros((n-1)*n/2,1);
iter = 1;
for i = 1:n
    for j = i+1:n
        FX(iter) =i-1;
        FY(iter) = j-1;
        iter = iter+1;
    end
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
tic
[fobjp, l1normX] = computLogDet(P, S, R, lambda, rho, params.penalty);
iterTime(1) = toc;
funVal(1) = fobjp;
fx = zeros((n-1)*n/2,1);
fy = zeros((n-1)*n/2,1);
tev = zeros(K,1);
iter = 1;
Gtf = zeros(size(P));
pgfm = -1;
tau = 1/params.AdaptivePar;
adpeta = params.NewtonEta;

% main loop
while (1)
    iter = iter+1;
    tic
    
    % compute B
    if iter ~= 2
        for k = 1:K
            invR = inv(R(:,:,k));
            B(:,:,k) = invR*invR';
        end
    end
    
    G = S - B;
    numF = findFreeSet(G,P,fx,fy,lambda,rho,tev,params.penalty);
    params.NumF = numF;
    
    dflag = 0;
    % compute subopt
    ztemp = P - tau*(S - B);
    mgl_prox(Gtf,ztemp,FX,FY,lambda*tau, rho*tau, length(FX), params.penalty);
    gfm = norm(P(:) - Gtf(:));
    subOpt(iter-1) = gfm;
    if params.Adaptive
        if iter == 2
            adpeta = params.NewtonEta;
        else
            adpeta = min(norm(pGtf(:) - Gtf(:))/pgfm, params.NewtonEta);
        end
    end
    
    [D,~,~,pGtf] = mgl_spg_adp(S,P*1.0,params,B,fx,fy, l1normX, adpeta*gfm);
    pgfm = gfm;
    
    if norm(D(:),'inf') < params.Dtol 
        fprintf('Step size is too small: ddr %5.2e\n', norm(D(:),'inf'));
        dflag = 1;
    end
    
    %linesearch
    alpha = 1;
    PW = (P+D).*marks;
    l1normXD = sum(abs(PW(:)));
    
    switch params.penalty
        case 0
            PW = PW(:,:,1:K-1) - PW(:,:,2:K);
            PW = sum(abs(PW),3);
        case 1
            PW = PW.*PW;
            PW = sqrt(sum(abs(PW),3));
        otherwise
            PW = PW(:,:,1:K-1) - PW(:,:,2:K);
            PW = sum(abs(PW),3);
    end
    l1normXD = lambda*l1normXD + rho*sum(PW(:));
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
        
        [fobj, l1normX1] = computLogDet(W, S, R, lambda, rho, params.penalty);
        if params.display
            fprintf('outer:%d inner%d, Linear search decreases %5.2e and norm of D is %5.2e, number of free set %d\n', iter, liter,fobjp - fobj, norm(D(:),'inf'), numF);
        end
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
    
    iterTime(iter)=toc;
    funVal(iter) = fobj;
    
    if dflag && params.NewtonStop~=4
        break;
    end
    
    fobjp = fobj;

    if params.NewtonStop == 1 && isfield(params,'NewtonObj')
        if funVal(iter) <= params.NewtonObj || abs(funVal(iter) - params.NewtonObj) < 1e-10
            break;
        end
    elseif params.NewtonStop == 2
        if subOpt(end) <= params.Newtontol
            break;
        end
    elseif params.NewtonStop == 3
        if(abs(funVal(iter) - funVal(iter-1)) <= params.Newtontol*abs(funVal(iter-1))) || norm(D(:),'inf') <  params.Dtol
            break;
        end
    end
    
    if iter>params.maxiter
        break;
    end
end

for k = 1:K
    invR = inv(R(:,:,k));
    B(:,:,k) = invR*invR';
end
mgl_prox(Gtf, P - tau*(S - B),FX,FY,lambda*tau, rho*tau, length(FX), params.penalty);
gfm = norm(P(:) - Gtf(:));
subOpt(iter) = gfm;

end

function [fobj, l1normX] = computLogDet(W, S, R, lambda, rho, penalty)
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
l1normX = sum(abs(W(:)));
PW = W;
switch penalty
    case 0
        PW = PW(:,:,1:K-1) - PW(:,:,2:K);
        PW = sum(abs(PW),3);
    case 1
        PW = PW.*PW;
        PW = sqrt(sum(abs(PW),3));
    otherwise
        PW = PW(:,:,1:K-1) - PW(:,:,2:K);
        PW = sum(abs(PW),3);
end
l1normX = lambda*l1normX + rho*sum(PW(:));

fobj = trdg - logdet + l1normX;
end


function opterr = objcompute(D,hatP, B, S, K, lambda, rho, penalty)
global marks
hatP = hatP.*marks;
opterr = 0;
for k = 1:K
    Wz = B(:,:,k)*D(:,:,k);
    opterr = opterr + Wz(:)'*Wz(:)/2 + ...
        trace((S(:,:,k) - B(:,:,k))*D(:,:,k));
end
opterr = opterr + lambda*sum(abs(hatP(:)));
PW = hatP;
switch penalty
    case 0
        PW = PW(:,:,1:K-1) - PW(:,:,2:K);
        PW = sum(abs(PW),3);
    case 1
        PW = PW.*PW;
        PW = sqrt(sum(abs(PW),3));
    otherwise
        PW = PW(:,:,1:K-1) - PW(:,:,2:K);
        PW = sum(abs(PW),3);
end
opterr = opterr + rho*sum(PW(:));
end

function [W,funVal,lambdak,Gtf] = mgl_spg_adp(S,P,params,B,fx,fy, pqobj,adgfm)

global FX FY

tau = 1/params.AdaptivePar;
lambdamin = 1e-30;
lambdamax = 1e30;
gamma = 1e-4;
eta = 2;
M = 10;
f_1 = -inf*ones(M,1);

K = size(S,3);
l1reg = params.lambda;
freg = params.rho;

SB = S - 2*B;
grad = SB;
hatP = 1.0*P; P_start = 1.0*P;
iter = 1;
for k = 1:K
    grad(:,:,k) = B(:,:,k)*hatP(:,:,k)*B(:,:,k);
end
grad  = grad + SB;
z = hatP - grad;
Gtf = zeros(size(z));
mgl_prox(P,z,fx,fy,l1reg, freg, params.NumF,params.penalty);
f_1(iter) = objcompute(hatP-P_start,hatP, B, S, K, l1reg, freg,params.penalty);
funVal(iter) = f_1(iter);

pomga = P - hatP;
err = norm(pomga(:),'inf');
% if ~params.Adaptive
%     if err <= params.SPGreltol*norm(hatP(:),'inf');
%         W = P - P_start;
%         return;
%     end
% end

if isfield(params,'SPGlambdak')
    lambdak =  params.SPGlambdak;
else
    lambdak = min(1/(err + 1e-15),10);
end
P = 1.0*hatP;
iter = iter+1;
count = 0;
while (iter < params.SPGmaxiter)
    
    cont_inner = 1;
    while cont_inner
        z = hatP - 1/lambdak*grad;
        %%% Feb 10, fix a bug?
        if l1reg/lambdak < 1e-16 && freg/lambdak <1e-16
            P = z;
        else
            mgl_prox(P,z,fx,fy,l1reg/lambdak, freg/lambdak, params.NumF, params.penalty);
        end
        D = P - hatP;
        fnew = objcompute(P-P_start,P, B, S, K, l1reg, freg,params.penalty);
        ddr = D(:)'*D(:);
        f_threshold = max(f_1) - 0.5*gamma*lambdak*ddr;
        if fnew <= f_threshold
            cont_inner=0;
        else
            lambdak = eta*lambdak;
        end
    end
    f_1(mod(iter,M)+1) = fnew;
    funVal(iter) = fnew;
    hatgrad = grad;
    for k = 1:K
        grad(:,:,k) = B(:,:,k)*P(:,:,k)*B(:,:,k);
    end
    grad  = grad + SB;
    gk = grad - hatgrad;
    td = D(:)'*gk(:);
    td = td/ddr;
    
    if td == 0
        fprintf('step size is too small...\n');
        break;
    end
    if td <=0
        lambdak = lambdamax;
        fprintf('Safe guard is activated...\n')
        continue
    end

    lambdak = max(lambdamin,min(td,lambdamax));
    
    iter = iter+1;
    
    
    if ~params.Adaptive
        pomga = P - hatP;
        err = norm(pomga(:),'inf')/norm(P(:),'inf');
        
        if err <= params.SPGreltol
            count = count + 1;
            if count > 1
            break;
            end
        end
    elseif params.Adaptive == 1
        
        if  fnew < pqobj
            ztemp = P - tau*grad;
            mgl_prox(Gtf,ztemp,FX,FY,l1reg*tau, freg*tau, length(FX), params.penalty);
            gtfm = norm(Gtf(:) - P(:));
            if gtfm <= adgfm
                fprintf('gtfm %5.2e ddr %5.2e\n',gtfm, ddr);
                break;
            end
        end
    else
        ztemp = P - tau*grad;
        mgl_prox(Gtf,ztemp,FX,FY,l1reg*tau, freg*tau, length(FX), params.penalty);
        gtfm = norm(Gtf(:) - P(:));
        if gtfm <= params.SPGreltol
            fprintf('Accurate: gtfm %5.2e ddr %5.2e\n', gtfm, ddr);
            break;
        end
    end
    
    hatP = 1.0*P;
end
fprintf('Number iterations %d\n',iter)
W = P - P_start;
end

