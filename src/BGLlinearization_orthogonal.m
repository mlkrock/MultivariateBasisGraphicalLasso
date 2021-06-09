function [S_star_array] = BGLlinearization_orthogonal(Qarray,inv_error_variances,projected_data)

p = length(inv_error_variances);
l = size(Qarray,3);
m = size(projected_data,2);

idiag=diag(inv_error_variances);
S_star_array = Qarray;

for i=1:l
    blockindices = (i-1)*p + (1:p);
    
    ll = inv(Qarray(:,:,i)+idiag);
    tt = ll*projected_data(blockindices,:);
    S_star_array(:,:,i) = ll + tt*tt'/m;
    
end

end

%
% for i=1:l
%     Qcell{i}  = Qcell{i}+ldiag;
%     Qcell{i} = inv(Qcell{i});
%
%     blockindices = (i-1)*p + (1:p);
%     tracecross = Qcell{i} * projected_data(blockindices,:);
%
%     S_star_array(:,:,i) = Qcell{i} + (tracecross * tracecross')/m;
% end