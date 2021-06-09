function [Qarray, counter] = BGLblockdiag_nopenalty(inv_error_variances,projected_data)


p = length(inv_error_variances);
l = size(projected_data,1)/p;


Qarray = zeros([p p l]);
for i = 1:l
    Qarray(:,:,i) = eye(p);
end

tol=0.05;
my_norm=1;
counter=1;
max_iterations=1000;


while my_norm >= tol && counter <= max_iterations
    
    
    fprintf('Starting Iter: %d', counter);
    fprintf('...linearization...');
    
    S_star_array = BGLlinearization_orthogonal(Qarray,inv_error_variances,projected_data);

    sum_sq = zeros([l 1]);
    sum_sq_diff = zeros([l 1]);
    
    fprintf('...BGL...\n')
    
    for i=1:l
        
        Qguess = inv(S_star_array(:,:,i));
        sum_sq_diff(i) =  norm(Qguess - Qarray(:,:,i),'fro')^2;
        sum_sq(i) =  norm(Qarray(:,:,i),'fro')^2;
        Qarray(:,:,i) = Qguess;
        
        
    end
    
    my_norm = sqrt(sum(sum_sq_diff))/sqrt(sum(sum_sq));
    
    fprintf('\n');
    fprintf('Ending Iter: %d, norm: %.12e\n', counter, my_norm );
    fprintf('\n');
    
    counter=counter+1;
    
end

end
