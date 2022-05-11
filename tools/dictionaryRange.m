function [H,rs,N] = dictionaryRange(f,r,r0,uk,rMinMax,Res)
%[H,rs,N] = dictionaryRange(f,r,r0,uk,rMinMax,Res) Obtains the H matrix for
% a sound field at positions given by r, comprised of spherical waves
% generated by point sources placed along the line at the direction given
% by the unit vector uk and the limits rMinMax. The monopoles are located
% at a distance Res apart from each other.
%   Input:
%       - f         : frequency (Hz). Nf x 1
%       - r         : point where the dictionary is created. 3 x M
%       - r0        : center of the spherical array. 3 x 1
%       - uk        : propagation unit vector. 1 x 3
%       - rMinMax   : minimum and maximum distances to the source. 2 x 1
%       - Res       : spatial resolution in meters. Scalar
%   Output:
%       - H         : dictionary. N x Nf
%       - rs        : position of the point sources. 3 x N
%       - N         : number of point sources. Scalar
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
if nargin < 5, error('dictionaryRange Error: Not enough input parameters.'), end

%% MAIN CODE
c = 343;

% Position of the candidate sources
rs = r0+uk'*(rMinMax(1):Res:rMinMax(2));

N = size(rs,2);
M = size(r,2);
Nf = length(f);

% Propagation vector
k = c./(2*pi*f);

% Distance to the candidate sources
d = nan(M,N);
for nn = 1:N
    for mm = 1:M
        d(mm,nn) = vecnorm(rs(:,nn)-r(:,mm));
    end
end

% Dictionary
H = nan(M,N,Nf);
for ff = 1:Nf
    H(:,:,ff) = (1./d).*exp(1i*d*k(ff));
%     H(:,:,ff) = exp(1i*d*k(ff));
end

% Normalisation
H = H./vecnorm(H,2,2);

end