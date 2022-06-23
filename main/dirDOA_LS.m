function DOA = dirDOA_LS(Data,Direct,Dict)
%DOA = dirDOA_LS(Data,Direct,Dict,plotFlag) Applies Least Squares
%solution to the Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%   Output:
%       - DOA        : DOA estimation and dictionary via LS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
if nargin < 3, error('dirDOA_LS Error: Not enough input parameters.'), end

%% MAIN CODE
DOA.Est = nan(3,length(Dict.f));
for ii = 1:length(Dict.f)
    x = pinv(Dict.H(:,:,ii))*Direct.InnSph.H(Data.f==Dict.f(ii),:).';
    
    [~,Idx] = max(abs(x));
    DOA.Est(:,ii) = -Dict.uk(:,Idx);
end

disp('Direct sound: DOA - LS... OK')

end

