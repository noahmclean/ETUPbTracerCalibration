function negLogLik = likelihoodFunction_meanOD(maxLikSoln, data, covmats)
%LIKELIHOODFUNCTION_MEANOD Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    maxLikSoln (9,1)   double
    data       (3,:)   double
    covmats    (3,3,:) double
end

arguments (Output)
    negLogLik  (1,1)   double
end

n = size(data,2);
mu = maxLikSoln(1:3);

L = [maxLikSoln(4), 0,             0; 
     maxLikSoln(5), maxLikSoln(6), 0; 
     maxLikSoln(7), maxLikSoln(8), maxLikSoln(9)];
Xi = L*L';

% Initialize the negative log-likelihood
negLogLik = 0;

% Compute the likelihood based on the provided data and covariance matrices
for i = 1:n
    
    residual = data(:, i) - mu; % Calculate the difference
    covOD = covmats(:,:,i) + Xi;
    negLogLik = negLogLik + residual'*(covOD\residual) + log(det(covOD));

end % function likelihoodFunction_meanOD()