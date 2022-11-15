function DOA = earlyDOA_RLS(Data,Early,Dict)
%DOA = earlyDOA_RLS(Data,Early,Dict) Applies Regularised Least Squares
%solution with Tikhonov regularisation and L-Curve method to the Early
%Reflections to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early reflections. Structure
%       - Dict      : dictionary of plane waves. Structure
%   Output:
%       - DOA       : DOA estimation and dictionary via RLS. Structure
%
% Author: Antonio Figueroa Durán
% Date: June 2022

%% ERROR HANDLING
if nargin < 3, error('earlyDOA_RLS Error: Not enough input parameters.'), end

%% MAIN CODE
x = nan(Dict.N,length(Dict.f));
for ii = 1:length(Dict.f)
    [x(:,ii),~] = reguLeastSquares(squeeze(Dict.H(:,:,ii)),Early.InnSph.H(Data.f==Dict.f(ii),:).');
end
DOA.x = sum(abs(x),2);

disp('Early reflections: DOA - RLS... OK')

end

