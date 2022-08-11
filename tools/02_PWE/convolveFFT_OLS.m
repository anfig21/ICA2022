function y = convolveFFT_OLS(x,h,N,bZeroPad)

if nargin < 2, error('Incorrect number of inmputs arguments.'), end
if isempty(x), error('x must be a vector'), end
if isempty(h), error('h must be a vector'), end
if nargin < 4 || isempty(bZeroPad), bZeroPad = 'false'; end

%% INITIAL PARAMETERS
x = x(:);
h = h(:);

M = length(h);
Ny = length(x);

%% ZERO-PADDING
% Pre zero-padding
x = vertcat(zeros(M-1,1),x);

% Post zero-padding
if bZeroPad, x = vertcat(x,zeros(M-1,1)); end

Nx = length(x);

if nargin < 3 || isempty(N), [N,~] = optimalN(Nx,M); end
L = N-M+1;

% Number of frames
nFrames = ceil(Nx/L);

% Zero-pad impulse response
h = vertcat(h,zeros(L-1,1));

% Zero-pad at the end of the signal
x = vertcat(x,zeros(nFrames*L-Ny,1));

%% DFT
H = fft(h,N);

%% SEGMENTATION
for ii = 0:nFrames-1
    xm = x((1:L+M-1) + ii*L);
    Xs = fft(xm,N);
    Ys = Xs.*H;
    yS = ifft(Ys,N);
    y(ii*L+1:(ii+1)*L) = yS(M:M+L-1);
end

if bZeroPad, y = y(1:Ny+M-1); else, y = y(1:Nx); end


end

