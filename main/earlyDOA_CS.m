function DOA = earlyDOA_CS(Data,Early,Dict,plotFlag)
%DOA = earlyDOA_CS(Data,Direct,Dict,plotFlag) Applies Compressive
%Sensing to the early reflections to obtain a reflection map.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early reflections. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - DOA        : DOA estimation and dictionary via CS. Structure
%
% Author: Antonio Figueroa Durán
% Date: April 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('earlyDOA_CS Error: Not enough input parameters.'), end

%% MAIN CODE
DOA.x = nan(Dict.Plane.N,length(Dict.f));
Nnorm = 1.1*Early.InnSph.NnormLcurve;

c = waitbar(0,'Loading...0\%','Name','CVX across frequencies...');
for ii = 1:length(Dict.f)
    Hii = squeeze(Dict.Plane.H(:,:,ii));
    pii = Early.InnSph.H(Data.f==Dict.f(ii),:).';
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
    variable x(Dict.Plane.N) complex;
    minimize norm(x,1);
    subject to
    norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    DOA.x(:,ii) = x;
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

%% PLOT
if plotFlag
    
    % Coefficients as a function of theta & phi
    if length(Dict.f) == 1
        [TH,PHI,~] = cart2sph(Dict.Plane.uk(1,:),Dict.Plane.uk(2,:),Dict.Plane.uk(3,:));
        
        figure
        plot(abs(Early.DOA.CS.x)), grid on
        xlabel('k Index'), ylabel('Coefficients Amplitude')
        applyAxisProperties(gca)
    end
end

end

