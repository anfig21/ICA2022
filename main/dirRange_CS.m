function CS = dirRange_CS(Data,Direct,CS,Dict,plotFlag)
%CS = dirRange_CS(Data,Direct,Dict,plotFlag) Applies Compressive
%Sensing to the Direct Sound Field to obtain the position of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - CS        : DOA data. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - CS        : Range estimation via via CS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 5, plotFlag = false;
elseif nargin < 4, error('dirRange_CS Error: Not enough input parameters.'), end

%% MAIN CODE
CS.Range.Est = nan(3,length(Dict.f));

c = waitbar(0,'Loading...0\%','Name','CVX across frequencies...');
for ii = 1:length(Dict.f)
    Hii = squeeze(Dict.Sph.H(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dict.f(ii),:).';
    Nnorm = 1.1*Direct.InnSph.Nnorm(Data.f==Dict.f(ii));
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
    variable x(Dict.Sph.N) complex;
    minimize norm(x,1);
    subject to
    norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    [~,Idx] = max(abs(x));
    CS.Range.Est(:,ii) = Dict.Sph.r(:,Idx);
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

CS.Range.Error = vecnorm(Direct.DOA.'-CS.Range.Est);

CS.Range.Avg = mean(CS.Range.Est,2);

%% PLOT
if plotFlag
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(Dict.Sph.r(1,:),Dict.Sph.r(2,:),Dict.Sph.r(3,:))
    scatter3(CS.Range.Est(1,:),CS.Range.Est(2,:),CS.Range.Est(3,:),100,'filled')
    scatter3(CS.Range.Avg(1),CS.Range.Avg(2),CS.Range.Avg(3),170,'filled')
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Atoms','Estimation','Mean')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

