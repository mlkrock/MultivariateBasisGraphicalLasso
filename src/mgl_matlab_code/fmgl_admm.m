function [P, funVal,Z,iterTime] = fmgl_admm(S,params)


lambda = params.lambda;
gamma = params.rho;

[n,n,K] = size(S);
P = S;
marks = zeros(n,n,K);
for k = 1:K
    marks(:,:,k) = 1 - eye(n,n);
end

rho = 5;
for i = 1:K
    P(:,:,i) = eye(n);
end
U = zeros(n,n,K);
Z = U;
Zo = Z;
D = zeros(size(P));
%%% for C
F = zeros((n-1)*n/2,2);
iter = 1;
for i = 1:n
    for j = i+1:n
        F(iter,:) =[i,j];
        iter = iter+1;
    end
end


for iter = 1:params.ADMMmaxiter
    tic
    W = -S + rho*(Z - U);
    for k = 1:K
        [V,Dz] = eig(W(:,:,k));
        wd = diag(Dz);
        D(:,:,k) = diag(wd + sqrt(wd.^2 + 4*rho))./(2*rho);
        P(:,:,k) = V*D(:,:,k)*V';
    end
    
    W = P + U;
    mgl_prox(Z,W,F(:,1)-1,F(:,2)-1,lambda/rho, gamma/rho,(n-1)*n/2,params.penalty);
    U = U + (P - Z);
    
    
    funVal(1,iter) = computLogDet(D, P, S, lambda, gamma,marks);
    
    if funVal(1,iter) <= params.admmobj || abs(funVal(1,iter) - params.admmobj)<params.admmobj*1e-10
        break;
    end
    %fprintf('maxiter %f...\n', funVal(1,iter) - params.admmobj);
    
    W = P - Z;
    funVal(2,iter) = norm(W(:),'inf');
    W = Z - Zo;
    funVal(3,iter) = norm(W(:),'inf');
    if(iter ==1)
        continue
    end
    
    if (abs(funVal(2,iter) - funVal(2,iter-1)) < params.admmtol && ...
           abs(funVal(3,iter) - funVal(3,iter-1)) < params.admmtol )
       fprintf('step size is too small %f, %f...\n', abs(funVal(2,iter) - funVal(2,iter-1)),  abs(funVal(3,iter) - funVal(3,iter-1)));
       break;
    end
    
    
    if(iter > params.ADMMmaxiter)
    fprintf('maxiter %f...\n', funVal(1,iter) - params.admmobj);
    end
    Zo = Z;
    iterTime(iter) = toc;
    if mod(iter,10) == 0
        fprintf('The %d ith and Time %f: ADMM %f obj %f\n', iter, sum(iterTime(1:iter)),funVal(1,iter),params.admmobj)
    end
end

end

function td = computLogDet(D, P, S, lambda, rho,marks)

td = 0;
K = size(P,3);
for i = 1:K
%    td = td -log(det(P(:,:,i))) + trace(S(:,:,i)*P(:,:,i));
   td = td -log(prod(diag(D(:,:,i)))) + trace(S(:,:,i)*P(:,:,i));
end
P = P.*marks;
td = td + lambda*sum(abs(P(:)));
P = P(:,:,1:K-1) - P(:,:,2:K);
td = td + rho*sum(abs(P(:)));
end


function td = computLogDetP(D, P, S, lambda, rho,marks)

td = 0;
K = size(P,3);
for i = 1:K
   td = td -log(det(P(:,:,i))) + trace(S(:,:,i)*P(:,:,i));
   %td = td -log(prod(diag(D(:,:,i)))) + trace(S(:,:,i)*P(:,:,i));
end
P = P.*marks;
td = td + lambda*sum(abs(P(:)));
P = P(:,:,1:K-1) - P(:,:,2:K);
td = td + rho*sum(abs(P(:)));
end
    