function DOA = earlyDOA_WL(Data,Early,Dict)
%DOA = earlyDOA_WL(Data,Early,Dict) Applies Weighted LASSO to the early
%reflections to obtain a reflection map.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early reflections. Structure
%       - Dict      : dictionary of plane waves. Structure
%   Output:
%       - DOA       : DOA estimation and dictionary via WL. Structure
%
% Author: Antonio Figueroa Durán
% Date: April 2022

%% ERROR HANDLING
if nargin < 3, error('earlyDOA_WL Error: Not enough input parameters.'), end

%% MAIN CODE
est = nan(Dict.N,length(Dict.f));
NoiseMargin = 10;           % dB
mu = 0.01;
weights = ones(Dict.N,1);

c = waitbar(0,'Loading...0\%','Name','earlyDOA_WL: CVX across frequencies...');
for ii = 1:length(Dict.f)
    Nnorm = 10^(NoiseMargin/20)*Early.InnSph.Nnorm(Data.f==Dict.f(ii));
    Hii = squeeze(Dict.H(:,:,ii));
    pii = Early.InnSph.H(Data.f==Dict.f(ii),:).';
    
    for jj = 1:4
        % CVX Formulation
        cvx_begin quiet
        cvx_precision high
            variable x(Dict.N) complex
            minimize norm(weights.*x,1)
        subject to
            norm(Hii*x-pii,2) <= Nnorm;
        cvx_end
        
        % Weights
        weights = 1./(abs(x)+mu);
    end
    
    est(:,ii) = x;
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

DOA.x = sum(abs(est),2);

disp('Early reflections: DOA - WL... OK')

end

