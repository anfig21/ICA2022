function Range = dirRange_RLS(Data,Direct,Dict,plotFlag)
%Range = dirRange_RLS(Data,Direct,Dict,plotFlag) Applies Regularised
%Least Squares solution with Tikhonov regularisation and L-Curve method to
%the Direct Sound Field to obtain the position of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - Range     : Range estimation via RLS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('dirRange_RLS Error: Not enough input parameters.'), end

%% MAIN CODE
Range.Est = nan(3,length(Dict.f));
for ii = 1:length(Dict.f)
    [x,~] = reguLeastSquares(squeeze(Dict.Sph.H(:,:,ii)),Direct.InnSph.H(Data.f==Dict.f(ii),:).');
    
    [~,Idx] = max(abs(x));
    Range.Est(:,ii) = Dict.Sph.r(:,Idx);
end

Range.Avg = mean(Range.Est,2);

if plotFlag
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(Dict.Sph.r(1,:),Dict.Sph.r(2,:),Dict.Sph.r(3,:))
    scatter3(Range.Est(1,:),Range.Est(2,:),Range.Est(3,:),150,'filled')
    scatter3(Range.Avg(1),Range.Avg(2),Range.Avg(3),150,'filled')
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Atoms','Estimation','Mean')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

disp('Direct sound: RANGE - RLS... OK')

end

