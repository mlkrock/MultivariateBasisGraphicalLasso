function P = genPric(n,K,num,L)
TT = generateTree(K);
P = zeros(n,n,K);
for i = 1:K
    P(:,:,i) = 0.25*eye(n);
end

ns = fix(n/L);
Z = zeros(n,n);
for i =1:L-1
    Z((i-1)*ns+1:i*ns,(i-1)*ns+1:i*ns) = ones(ns);
end

if L == 1
  Z = ones(n);
else
  Z((L-1)*ns+1:n,(L-1)*ns+1:n) = ones(n-ns*(L-1));
end


Z = tril(Z,-1);
idx = find(Z~=0);
ridx = randsample(length(idx),num);
idx = idx(ridx);
[I,J] = ind2sub([n,n],idx);
for i = 1:length(I)
    ri = I(i);
    ci = J(i);
    T = sampleTree(K,TT);
    for it = 1:size(T,2)
        for sit = T(1,it):T(2,it)
           sigma = rand(1)/4+0.1;
           P(ri,ci,sit) = P(ri,ci,sit) - sigma;
           P(ci,ri,sit) = P(ri,ci,sit);
           P(ri,ri,sit) = P(ri,ri,sit) + sigma;
           P(ci,ci,sit) = P(ci,ci,sit) + sigma;
        end
    end
end