function T = sampleTree(n,TT)

T = TT;
nT = size(T,2);
nt = randsample(nT,1);

nt = min(nT/2,nt);

idx = randsample(nT,nt);

idx = sort(idx);
Tz = T(:,idx);

[zd,zidx] = sort(Tz(1,:));
Tz = Tz(:,zidx);
count = size(Tz,2);
i = 1;
T = Tz;
while(i <= count)
j = i+1;
while (j <= count)
    if ~(T(1,i) >T(2,j) || T(2,i) < T(1,j))
        T(:,i) = [min(T(1,i),T(1,j));max(T(2,i),T(2,j))];
        T(:,j) = [];
        count = count - 1;
    else
        j = j + 1;
    end
end
i = i + 1;
end
end


