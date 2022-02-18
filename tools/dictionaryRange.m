function [H,uk] = dictionaryRange(f,r0,uk,rMinMax,N)
%[H,uk] = dictionary(f,N) Obtains the H matrix for a sound field at a
% position r0 comprised of N spherical waves with propagation vector fixed
% from the given uk and distances to the source between rMinMax.
%   Input:
%       - f         : frequency (Hz). Nf x 1
%       - r0        : point where the dictionary is created. 3 x 1
%       - uk        : propagation unit vector. 1 x 3
%       - rMinMax   : minimum and maximum distances to the source. 2 x 1
%       - N         : number of plane waves. Scalar
%   Output:
%       - H         : dictionary. N x Nf
%       - uk        : propagation unit vector. 3 x N
%
% Author: Antonio Figueroa Durán
% Date: February 2022

%% ERROR HANDLING
if nargin < 5, error('dictionaryRange Error: Not enough input parameters.'), end

%% MAIN CODE
c = 343;

% Propagation vector
k = c./(2*pi*f);
kv = uk'*k;

% Distance to the source
r = uk'*linspace(rMinMax(1),rMax(2),N);
d = r0-r;

% Dictionary
H = (1./vecnorm(d).').*exp(-1i*d'*kv);

% Normalisation
H = H./vecnorm(H,2,2);

end