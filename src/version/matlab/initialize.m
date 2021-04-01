% v 1.0.2
function [sigmaPoints,weights,lambda,covariance,covarianceZ] = initialize(n,alpha,kappa,nmo)
%==============================================================
% Initialize variables to Weights as well as Sigma points
%==============================================================

sigmaPoints = zeros(2*n+1 , n);
weights = zeros(2*n + 1,1);
lambda = (alpha^2) * (n + kappa) - n;
covariance = zeros(n);
covarianceZ = zeros(nmo);

end