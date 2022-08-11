function DOA = dirDOA_CS(Data,Direct,Dict)
%DOA = dirDOA_CS(Data,Direct,Dict,plotFlag) Applies Compressive
%Sensing to the Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%   Output:
%       - DOA        : DOA estimation and dictionary via CS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
if nargin < 3, error('dirDOA_CS Error: Not enough input parameters.'), end

%% MAIN CODE
DOA.Est = nan(3,length(Dict.f));
NoiseMargin = 10;           % dB

c = waitbar(0,'Loading...0\%','Name','dirDOA_CS: CVX across frequencies...');
for ii = 1:length(Dict.f)
    Nnorm = 10^(NoiseMargin/20)*Direct.InnSph.Nnorm(Data.f==Dict.f(ii));
    Hii = squeeze(Dict.H(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dict.f(ii),:).';
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
        variable x(Dict.N) complex;
        minimize norm(x,1);
        subject to
            norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    [~,Idx] = max(abs(x));
    
    DOA.Est(:,ii) = -Dict.uk(:,Idx);
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

disp('Direct sound: DOA - CS... OK')

end

