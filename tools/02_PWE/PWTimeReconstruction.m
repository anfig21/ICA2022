function h = PWTimeReconstruction(Data,x,t,N,R)
%h = PWTimeReconstruction(Data,x,t,N,R) Reconstructs the sound field at
%points R using the coefficients x.
%   Input:
%       - Data      : raw data. Structure
%       - x         : PWE coefficients. Nt*N x 1
%       - t         : time vector (s). Nt x 1
%       - N         : number of plane waves. Scalar
%       - R         : reconstruction points. 3 x Nr

%   Output:
%       - h         : reconstructed RIR. Nt x Nr
%
% Author: Antonio Figueroa Durán
% Date: July 2022


%% ERROR HANDLING
if nargin < 5, error('PWTimeReconstruction Error: Not enough input parameters.'), end

%% MAIN CODE
Nr = size(R,2);
Nt = length(t);
h = zeros(Nt,Nr);

for ii = 1:Nr
    [HR,~] = dictionaryPWTime(Data.c,t,Data.Fs,R(:,ii),Data.Sph.R0.',N);
    h(:,ii) = HR*x;
end

end

