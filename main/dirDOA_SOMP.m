function DOA = dirDOA_SOMP(Direct,Dict)
%DOA = dirDOA_SOMP(Direct,Dict) Applies Simultaneous Orthogonal Matching
%Pursuit over the Direct Sound Field to obtain the position of the source.
%Disclaimer: DOA estimation uses plane waves. Therefore the model is not
%accurate enough to obtain reliable results.
%   Input:
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%   Output:
%       - DOA       : DOA estimation and dictionary. Structure
%
% Author: Antonio Figueroa Durán
% Date: February 2022

%% ERROR HANDLING
if nargin < 2, error('dirDOA_SOMP Error: Not enough input parameters.'), end

%% MAIN CODE
x = somp(Dict.H, Direct.InnSph.H.', 1);
Idx = find(x,1);
DOA.Est = -Dict.uk(:,Idx);

disp('Direct sound: DOA - SOMP... OK')

end

