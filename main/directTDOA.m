function TDOA = directTDOA(Data,Direct,plotFlag)
%TDOA = directTDOA(Data,Direct) Applies Time Difference Of Arrival (TDOA)
%to estimate the DOA and Range to the source.
%Ref:
%   [1] Tervo, Sakari, et al. “Spatial Decomposition Method for Room
%       Impulse Responses.” Aes: Journal of the Audio Engineering Society,
%       vol. 61, no. 1-2, AUDIO ENGINEERING SOC, 2013, pp. 17–28
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - plotFlag  : 'true' to plot setup & TDOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - TDOA      : DOA & Range estimation and dictionary. Structure
%
% Author: Antonio Figueroa Durán
% Date: February 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 3, plotFlag = false;
elseif nargin < 2, error('directTDOA Error: Not enough input parameters.'), end

%% MAIN CODE






end

