function [negLP,grad,H] = neglogposterior(prs,negloglifun,Cinv, linearFilterLength,lambda)
% [negLP,grad,H] = neglogposterior(prs,negloglifun,Cinv)
%
% Compute negative log-posterior given a negative log-likelihood function
% and zero-mean Gaussian prior with inverse covariance 'Cinv'.
%
% Inputs:
 %   prs [d x 1] - parameter vector
%    negloglifun - handle for negative log-likelihood function
%   Cinv [d x d] - response (spike count per time bin)
%
% Outputs:
%          negLP - negative log posterior
%   grad [d x 1] - gradient 
%      H [d x d] - Hessian (second deriv matrix)

% Compute negative log-posterior by adding quadratic penalty to log-likelihood
linearFilter = prs(1:linearFilterLength);
L2Smooth = [0 diff(linearFilter)];
switch nargout
        
    case 1  % evaluate function
        negLP = negloglifun(prs) + .5*sum(L2Smooth.^2) * lambda;
    
    case 2  % evaluate function and gradient
        [negLP,grad] = negloglifun(prs);
        negLP = negloglifun(prs) + .5*sum(L2Smooth.^2) * lambda;
        grad(1:linearFilterLength) = grad(1:linearFilterLength) + L2Smooth * lambda;

    case 3  % evaluate function and gradient
        [negLP,grad,H] = negloglifun(prs);
        negLP = negloglifun(prs) + .5*sum(L2Smooth.^2) * lambda;
        grad(1:linearFilterLength) = grad(1:linearFilterLength) + L2Smooth * lambda;
        H(1:linearFilterLength,1:linearFilterLength)  = H(1:linearFilterLength,1:linearFilterLength) + Cinv;
end

