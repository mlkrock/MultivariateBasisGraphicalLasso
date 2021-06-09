function params = setParms(params)
% function params = setParms(params)
% set up the default params


if ~isfield(params,'maxlineiter')
    params.maxlineiter = 50; % the max number of iteration in line search
end
if ~isfield(params,'SPGreltol')
    params.SPGreltol = 1e-6;  % the tolerance of SPG
end
if ~isfield(params,'Newtontol')
    params.Newtontol = 1e-6;   % the tolerance of newton method
end

if ~isfield(params,'maxiter')
    params.maxiter = 1500;     % the max number of iteration
end

if ~isfield(params,'sigma')
    params.sigma = 1e-3;
end

if ~isfield(params,'SPGmaxiter')
    params.SPGmaxiter = 150;
end

if ~isfield(params,'Adaptive')
    params.Adaptive = 1;
end

if ~isfield(params,'AdaptivePar')
    params.AdaptivePar = 100;
end

if ~isfield(params, 'NewtonStop')
    params.NewtonStop = 3;
end

if ~isfield(params, 'NewtonEta')
    params.NewtonEta = 0.1;
end

if ~isfield(params, 'penalty')
    params.penalty = 0;
end

if ~isfield(params, 'display')
    params.display = 0;
end

if ~isfield(params, 'Dtol')
    params.Dtol = 1e-8;
end
