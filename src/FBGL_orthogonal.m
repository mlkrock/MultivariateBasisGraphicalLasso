function [Qarray, counter] = FBGL_orthogonal(lambda,rho,inv_error_variances,projected_data,Qarray)

p = length(inv_error_variances);
l = size(projected_data,1)/p;

%guess should be sparse, here from speye

params.sigma = 1e-3;
params.maxlineiter = 50;
params.maxiter = 50;
params.linesearch = true;
params.SPGmaxiter = 50;
params.SPGreltol = 1e-8;
params.Newtontol = 1e-5;
% NewtonStop 1: only works if you set an objective value, the
% algorithm won't stop util the objective value is smaller than the
% specified one
% NewtonStop 2: stop when subgradient optimial condition, which can be
% checked from output
% NewtonStop 3: stop when relative error is less than Newtontol
params.NewtonStop = 2;
% use adaptive stop condition in NSPG
params.Adaptive = 1;

% fused penalty
params.penalty = 0;
% group penalty
%params.penalty = 1;

params.display=0;



tol=0.05;
my_norm=1;
counter=1;
max_iterations=100;
params.lambda = lambda;
params.rho = rho;




while my_norm >= tol && counter <= max_iterations
    
    
    fprintf('Starting Iter: %d', counter);
    fprintf('...linearization...');
    
    S_star_array = BGLlinearization_orthogonal(Qarray,inv_error_variances,projected_data);
    
    fprintf('...FMGL...\n')
    
    %[new_guess,funVal,iterTime,subopt] = mgl_adp(S_star_array,params);
    
    params.P0 = Qarray;
    [new_guess,~] = mgl(S_star_array,params);
    
    
    sum_sq = zeros([l 1]);
    sum_sq_diff = zeros([l 1]);
    
    
    for i=1:l
        
        Qguess =  new_guess(:,:,i);
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
