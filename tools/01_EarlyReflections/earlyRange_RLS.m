function Range = earlyRange_RLS(Data,Early,Dict)
%Range = earlyRange_RLS(Data,Early,Dict) Applies Regularised Least Squares
%solution with Tikhonov regularisation and L-Curve method to the early
%reflections to obtain the position of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early reflections. Structure
%       - Dict      : dictionary of plane waves. Structure
%   Output:
%       - Range     : Range estimation via RLS. Structure
%
% Author: Antonio Figueroa Durán
% Date: April 2022

%% ERROR HANDLING
if nargin < 3, error('earlyRange_RLS Error: Not enough input parameters.'), end

%% MAIN CODE
Range.Est = nan(3,length(Dict.f));
for ii = 1:length(Dict.f)
    [x,~] = reguLeastSquares(squeeze(Dict.H(:,:,ii)),Early.InnSph.H(Data.f==Dict.f(ii),:).');
    
    [~,Idx] = max(abs(x));
    Range.Est(:,ii) = Dict.r(:,Idx);
end

disp('Early reflections: RANGE - RLS... OK')

end

