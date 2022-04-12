function DOA = dirDOA_SOMP(Data,Direct,Dict,plotFlag)
%DOA = dirDOA_SOMP(Data,Direct) Applies Simultaneous Orthogonal Matching
%Pursuit over the Direct Sound Field to obtain the position of the source.
%Disclaimer: DOA estimation uses plane waves. Therefore the model is not
%accurate enough to obtain reliable results.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - DOA       : DOA estimation and dictionary. Structure
%
% Author: Antonio Figueroa Durán
% Date: February 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('dirDOA_SOMP Error: Not enough input parameters.'), end

%% MAIN CODE
% SOMP
DOA.x = somp(Dict.Plane.H, Direct.InnSph.H.', Dict.Plane.K);
DOA.Idx = find(DOA.x,Dict.Plane.K);
DOA.Est = Dict.Plane.uk(:,DOA.Idx);

if plotFlag
    uk = Dict.Plane.uk+Data.Sph.R0.';
    
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(uk(1,:),uk(2,:),uk(3,:))
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),-DOA.Est(1),-DOA.Est(2),-DOA.Est(3),2,'Linewidth',4)
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

DOA.Error = vecnorm(Direct.TrueDOA.'-DOA.Est);

end

