function [X1,X2] = fftUniBi(x)
%[X1,X2] = fftUniBi(x) Obtains the uni- and bilateral spectra of the signal
%x. The unilateral spectrum has half the number of samples of the original
%signal.
%   Input:
%       - x         : time-domain signal. Nsamples x M
%   Output:
%       - X1        : unilateral spectrum. Nsamples/2 x M
%       - X2        : bilateral spectrum. Nsamples x M
%
% Author: Antonio Figueroa Durán
% Date: March 2022

Nsamples = size(x,1);
X2 = fft(x)/Nsamples;
X1 = [X2(1,:); 2*X2(2:end/2,:)];

end

