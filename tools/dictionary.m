function [H,uk] = dictionary(c,f,r0,N)
%[H,uk] = dictionary(c,f,r0,N) Obtains the H matrix for a sound field at a
% position r0 comprised of N plane waves with propagation vector sampled
% along the surface of a unit sphere.
%   Input:
%       - c         : speed of sound (m/s). Scalar
%       - f         : frequency (Hz). Nf x 1
%       - r0        : measurement point. 3 x M
%       - N         : number of plane waves. Scalar
%   Output:
%       - H         : dictionary. M x N x Nf
%       - uk        : propagation unit vector. 3 x N
%
% Author: Antonio Figueroa Durán
% Date: February 2022

%% ERROR HANDLING
if nargin < 4, error('dictionary Error: Not enough input parameters.'), end

%% MAIN CODE
Nf = length(f);
[~,M] = size(r0);

% Spatial sampling - Fibonacci Algorithm
uk = -fibonacciSampling(N);      % uk point inwards the sphere

% If f contains 0 Hz, it is substituted by the next sample
if f(1) == 0, f(1) = f(2); end

% Propagation vector
k = c./(2*pi*f);
kv = nan(3,N,Nf);
for ii = 1:Nf, kv(:,:,ii) = k(ii)*uk; end

% Dictionary
H = nan(M,N,Nf);
d = waitbar(0,'Loading...0\%','Name','Building dictionary...');
for ii = 1:Nf
    H(:,:,ii) = exp(-1i*r0'*squeeze(kv(:,:,ii)));
    waitbar(ii/Nf,d,strcat("Loading... ",string(round(100*ii/Nf,2)),"\,\%"));
end
delete(d)

% Normalisation
H = H./vecnorm(H,2,2);

disp('Plane Wave Dictionary... OK')
end