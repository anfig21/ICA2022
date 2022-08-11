function h = reconstructionPWTime(c,t,Fs,r,r0,w)
%[H,uk] = dictionary(c,f,r0,N) Description
%   Input:
%       - c         : speed of sound (m/s). Scalar
%       - f         : frequency (Hz). Nf x 1
%       - r         : reconstruction points. 3 x M
%       - r0        : center of the spherical array. 3 x 1
%       - w         : plane wave coefficients. Nt x N
%   Output:
%       - h         : RIR. Nt x Nr
%
% Author: Antonio Figueroa Durán
% Date: June 2022

%% ERROR HANDLING
if nargin < 5, error('reconstructionPWTime Error: Not enough input parameters.'), end

%% MAIN CODE
Nt = length(t);
[norm,M] = size(r);
N = size(w,2);

% Spatial sampling - Fibonacci Algorithm
uk = fibonacciSampling(N);      % uk point inwards the sphere

% Position of the virtual source (far from the array)
d = 7;          % meters
kv = d*uk;
rs = r0+kv;

% Dictionary
h = zeros(Nt,M);
if M > 1, wbar = waitbar(0,'Loading...0\%','Name','Building dictionary...'); end
for ii = 1:M    % Mics
    r2s = r(:,ii)-rs;         % Vector source -> mic
    projection = sqrt(sum(r2s.*(-uk)));    % Projection onto propagation plane
    aux = (1/(4*pi*norm))*sinc(Fs*(t'-projection/c)/pi);      % sinc(x) = sin(pi*x)/(pi*x)
    
    for jj = 1:N
        h(:,ii) = h(:,ii)+toeplitz(aux(:,jj),[aux(1,jj) zeros(1,Nt-1)])*w(:,jj);
    end
    
    if M > 1, waitbar(ii/M,wbar,strcat("Loading... ",string(round(100*ii/M,2)),"\,\%")); end
end
if M > 1, delete(wbar), end

end

