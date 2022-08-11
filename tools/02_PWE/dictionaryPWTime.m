function [H,uk] = dictionaryPWTime(c,t,Fs,r,r0,N)
%[H,uk] = dictionaryPWTime(c,t,Fs,r,r0,N) Generates a dictionary of N plane
%waves in time domain via convolution matrix.
%   Input:
%       - c         : speed of sound (m/s). Scalar
%       - t         : time vector (s). Nt x 1
%       - Fs        : sampling frequency (Hz). Scalar
%       - r         : measurement point. 3 x M
%       - r0        : center of the spherical array. 3 x 1
%       - N         : number of plane waves. Scalar
%   Output:
%       - H         : dictionary. M x N x Nf
%       - uk        : propagation unit vector. 3 x N
%
% Author: Antonio Figueroa Durán
% Date: June 2022

%% ERROR HANDLING
if nargin < 5, error('dictionaryPWTime Error: Not enough input parameters.'), end

%% MAIN CODE
Nt = length(t);
[~,M] = size(r);

% Spatial sampling - Fibonacci Algorithm
uk = fibonacciSampling(N);      % uk point inwards the sphere

% Position of the virtual source
d = t(1)*c+0.1;          % meters
kv = d*uk;
rs = r0+kv;

% Dictionary
H = nan(Nt*M,Nt*N);
for ii = 1:M    % Mics
    r2s = r(:,ii)-rs;               % Vector source -> mic
    projection = dot(r2s,-uk);      % Projection onto propagation plane
%     aux = sinc(Fs*(t'-projection/c));
    aux = sinc(Fs*(t'-projection/c)/pi);
    for jj = 1:N
        H((ii-1)*Nt+(1:Nt),(jj-1)*Nt+(1:Nt)) = toeplitz(aux(:,jj),[aux(1,jj) zeros(1,Nt-1)]);
%         H((ii-1)*Nt+(1:Nt),(jj-1)*Nt+(1:Nt)) = toeplitz(aux(:,jj));
    end
end

end

