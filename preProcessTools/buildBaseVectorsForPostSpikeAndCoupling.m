% Make basis for spike history and coupling 
% (based on pillowLab work)
% -------
% Inputs: 
%            numOfVectors = numberOfVectorsTooBuild
%            hpeaks = 2-vector containg [1st_peak  last_peak], the peak 
%                     location of first and last raised cosine basis vectors
%            b = offset for nonlinear stretching of x axis:  y = log(x+b) 
%                     (larger b -> more nearly linear stretching)
%
%     dt = grid of time points for representing basis
%  --------
%  Outputs:  iht = time lattice on which basis is defined
%            ihbas = orthogonalized basis
%            ihbasis = original (non-orthogonal) basis 
%
function [iht, ihbas, ihbasis] = buildBaseVectorsForPostSpikeAndCoupling(numOfVectors,dt,hpeaks, b)

% Check input values
if (hpeaks(1)+b) < 0, 
    error('b + first peak location: must be greater than 0'); 
end

% nonlinearity for stretching x axis (and its inverse)
nlin = @(x)log(x+1e-20);
invnl = @(x)exp(x)-1e-20; % inverse nonlinearity

% Generate basis of raised cosines
yrnge = nlin(hpeaks+b);        % nonlinearly transformed first & last bumps
db = diff(yrnge)/(numOfVectors-1);    % spacing between cosine bump peaks
ctrs = yrnge(1):db:yrnge(2);   % centers (peak locations) for basis vectors
mxt = invnl(yrnge(2)+2*db)-b;  % maximum time bin
iht = (dt:dt:mxt)';
nt = length(iht);        % number of points in iht
ff = @(x,c,dc)(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2; % raised cosine basis vector
ihbasis = ff(repmat(nlin(iht+b), 1, numOfVectors), repmat(ctrs, nt, 1), db);

% compute orthogonalized basis
ihbas = orth(ihbasis);


%ihbas = zeros(size(ihbas,1),size(ihbas,2));