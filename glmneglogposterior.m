function [negLP,grad,H] = glmneglogposterior(prs,negloglifun,Cinv, linearFilterLength)
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
switch nargout
        
    case 1  % evaluate function
        negLP = negloglifun(prs)  + .5 * linearFilter * Cinv * linearFilter';
    
    case 2  % evaluate function and gradient
        [negLP,grad] = negloglifun(prs);
        negLP = negloglifun(prs)  + .5 * linearFilter * Cinv * linearFilter';
        grad(1:linearFilterLength) = grad(1:linearFilterLength) + (Cinv * linearFilter')';
        grad(1) = 0;

    case 3  % evaluate function and gradient
        [negLP,grad,H, Hk, Hkb, Hkh, Hb, Hhb, Hh] = negloglifun(prs);
        negLP = negloglifun(prs)  + .5 * linearFilter * Cinv * linearFilter';
        grad(1:linearFilterLength) = grad(1:linearFilterLength)  + (Cinv * linearFilter')';
        Hk = Hk + Cinv;
        H = [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]];
end

