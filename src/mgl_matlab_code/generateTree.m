function T = generateTree(n)
T =  generateTreeIT(n);
n = size(T,2);
idx = [];
for i = 1:n
    if T(1,i) == T(2,i)
       idx =[idx i];
    end
end
T(:,idx) = [];
end

function T = generateTreeIT(n)

if(n==1)
    T = [ ];
    return
end

T2 = generateTreeIT(n/2);    
if length(T2)==0
   T = [T2 [1;n/2] [n/2+1;n]];
else
    T = [T2 n/2+T2 [1;n/2] [n/2+1;n]];
end
end