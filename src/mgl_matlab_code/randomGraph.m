function Tran = randomGraph(p,n)
Tran  = tril(ones(p,p),-1)';
idx = find(Tran>0);
len = numel(idx);
eidx = randsample(len,n);
Tran=zeros(p,p);
Tran(idx(eidx))=1;
Tran = Tran+Tran';
