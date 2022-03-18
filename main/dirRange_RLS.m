function RLS = dirRange_RLS(Data,Direct,RLS,Dict,plotFlag)
%RLS = dirRange_RLS(Data,Direct,RLS,Dict,plotFlag) Applies Regularised
%Least Squares solution with Tikhonov regularisation and L-Curve method to
%the Direct Sound Field to obtain the position of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - RLS       : DOA data. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - RLS       : Range estimation via RLS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 5, plotFlag = false;
elseif nargin < 4, error('dirRange_RLS Error: Not enough input parameters.'), end

%% MAIN CODE
RLS.Range.Est = nan(3,length(Dict.f));
for ii = 1:length(Dict.f)
    [x,~] = reguLeastSquares(squeeze(Dict.Sph.H(:,:,ii)),Direct.InnSph.H(Data.f==Dict.f(ii),:).');
    
    [~,Idx] = max(abs(x));
    RLS.Range.Est(:,ii) = Dict.Sph.r(:,Idx);
end

RLS.Range.Avg = mean(RLS.Range.Est,2);

if plotFlag
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(Dict.Sph.r(1,:),Dict.Sph.r(2,:),Dict.Sph.r(3,:))
    scatter3(RLS.Range.Est(1,:),RLS.Range.Est(2,:),RLS.Range.Est(3,:),150,'filled')
    scatter3(RLS.Range.Avg(1),RLS.Range.Avg(2),RLS.Range.Avg(3),150,'filled')
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Atoms','Estimation','Mean')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

