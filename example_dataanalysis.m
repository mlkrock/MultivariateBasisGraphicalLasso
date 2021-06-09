%% Introduction and Setup Data

clc
clear all

addpath(genpath('src'));

% read in your data as a p x n x m 'data_array'
% number of variables x number of locations x number of realizations

p = size(data_array,1);
n = size(data_array,2);
m = size(data_array,3);

%% Basis Setup

% construct your orthogonal basis matrix 'smallPhi' here
% we use EOFs 
% be sure to save it if your basis setup is costly like this SVD
% smallPhi dimension is n x l where l = 2000 is the number of basis functions
% need smallPhi' * smallPhi = eye(l)

reshaped_data_array = zeros([n p*m]);
for i = 1:n
    reshaped_data_array(i,:) = vec(data_array(:,i,:));
end

[U,~,~] = svd(reshaped_data_array,'econ');
clear reshaped_data_array
smallPhi= U(:,1:2000);
l = size(smallPhi,2);

%% nugget estimate

% constant diagonal Q

x0 = [1,1];
lb = [0 0];
ub = [50 50];
nuggetestimates = zeros([p 1]);
smoothnessestimates = zeros([p 1]);

for i=1:p
    tt = smallPhi'* squeeze(data_array(i,:,:));
    smallPhi_S_Phi = tt*tt'/m;
    trS = norm(squeeze(data_array(i,:,:)),"fro")^2/m;
    [xval, fval]= fmincon(@(x) nuggetlikelihood_orthog(x,smallPhi_S_Phi,trS,n),x0,[], [], [], [], lb, ub);
    nuggetestimates(i) = xval(1);
    smoothnessestimates(i) = xval(2);
end

% exponentially varying diagonal Q for each level
% 
% x0 = [1,0.05,1];
% lb = [0 0 0];
% ub = [50 0.1 50];
% 
% nuggetestimates_exp = zeros([p 1]);
% smoothnessestimates_exp = zeros([p 1]);
% marginalvarianceestimates_exp = zeros([p 1]);
% 
% for i=1:p
%     tt = smallPhi'* squeeze(data_array(i,:,:));
%     smallPhi_S_Phi = tt*tt'/m;
%     trS = norm(squeeze(data_array(i,:,:)),"fro")^2/m;
%     [xval, fval]= fmincon(@(x) nuggetlikelihood_expdiag_orthog(x,smallPhi_S_Phi,trS,n),x0,[], [], [], [], lb, ub);
%     nuggetestimates_exp(i) = xval(1);
%     smoothnessestimates_exp(i) = xval(2);
%     marginalvarianceestimates_exp(i) = xval(3);
% end

% see supplementary material for more explanation here
% nugget estimates were the same out to several decimal places
% for both constant diagonal and exponentially varying Q
% in our analysis

inv_error_variances = 1./nuggetestimates;

%% projecting data and creating matrix

% this is creating the matrix
% Phi' * Dinv * Y_datamatrix
% with multiplications dimensions (pl x pn) * (pn x pn) * (pn x m)
% efficiently using a Kronecker product trick

projected_data = zeros(p*l,m);

for i=1:m
    projected_data(:,i)= vec((inv_error_variances .* data_array(:,:,i)) * smallPhi);
end

%% main algorithms

% no fusion
lambda=20;
graphguess = BGLblockdiag_orthogonal(lambda,inv_error_variances,projected_data);

% fusion
rho=10;
fusedguess = FBGL_orthogonal(lambda,rho,inv_error_variances,projected_data,graphguess);

% mle, no penalty
mleguess = BGLblockdiag_nopenalty(inv_error_variances,projected_data);
