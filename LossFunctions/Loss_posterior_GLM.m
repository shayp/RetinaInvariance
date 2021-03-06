function [negLP,grad,H] = Loss_posterior_GLM(postGlmParams,negloglifun,Cinv, linearFilterLength)
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
linearFilter = postGlmParams(1:linearFilterLength);
switch nargout
        
    case 1  % evaluate function
        negLP = negloglifun(postGlmParams)  + .5 * linearFilter' * Cinv * linearFilter;
    
    case 2  % evaluate function and gradient
        [negLP,grad] = negloglifun(postGlmParams);
        negLP = negLP  + .5 * linearFilter' * Cinv * linearFilter;
        grad(1:linearFilterLength) = grad(1:linearFilterLength) + Cinv * linearFilter;

    case 3  % evaluate function and gradient
        [negLP,grad,H] = negloglifun(postGlmParams);
        negLP = negLP  + .5 * linearFilter' * Cinv * linearFilter;
        grad(1:linearFilterLength) = grad(1:linearFilterLength)  + Cinv * linearFilter;
        H(1:linearFilterLength, 1:linearFilterLength) = H(1:linearFilterLength, 1:linearFilterLength) + Cinv;

end