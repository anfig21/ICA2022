function Range = dirRange_RLS(Data,Direct,Dict)
%Range = dirRange_RLS(Data,Direct,Dict) Applies Regularised Least Squares
%solution with Tikhonov regularisation and L-Curve method to the Direct
%Sound Field to obtain the position of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of spherical waves. Structure
%   Output:
%       - Range     : Range estimation via RLS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
if nargin < 3, error('dirRange_RLS Error: Not enough input parameters.'), end

%% MAIN CODE
Range.Est = nan(3,length(Dict.f));
for ii = 1:length(Dict.f)
    [x,~] = reguLeastSquares(squeeze(Dict.H(:,:,ii)),Direct.InnSph.H(Data.f==Dict.f(ii),:).');
    
    [~,Idx] = max(abs(x));
    Range.Est(:,ii) = Dict.r(:,Idx);
end

disp('Direct sound: RANGE - RLS... OK')

end

