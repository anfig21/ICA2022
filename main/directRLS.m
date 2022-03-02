function RLS = directRLS(Data,Direct,Dictionary,plotFlag)
%RLS = directRLS(Data,Direct,Dictionary,plotFlag) Applies Regularised Least
%Squares solution with Tikhonov regularisation and L-Curve method to the
%Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dictionary: dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - RLS       : DOA estimation and dictionary via LS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('directRLS Error: Not enough input parameters.'), end

%% MAIN CODE
% Regularised Least-Squares solution across frequency - 'l-curve' method
RLS.Est = nan(3,length(Dictionary.f));
for ii = 1:length(Dictionary.f)
    [RLS.U,RLS.s,RLS.V] = csvd(squeeze(Dictionary.CA(:,:,ii)));
    [RLS.lambda,~,~,~] = l_curve(RLS.U,RLS.s,Data.H(Data.f==Dictionary.f(ii),:).','Tikh',[],[],false);
    [RLS.x,~,~] = tikhonov(RLS.U,RLS.s,RLS.V,Data.H(Data.f==Dictionary.f(ii),:).',RLS.lambda);
    
    [~,RLS.Idx] = max(abs(RLS.x));
    RLS.Est(:,ii) = Dictionary.uk(:,RLS.Idx);
end

RLS.Error = vecnorm(Direct.DOA.'-RLS.Est);

RLS.Avg = -mean(RLS.Est,2);
RLS.Avg = RLS.Avg/vecnorm(RLS.Avg);

if plotFlag
    % Mean Squared Error
    figure, plot(Dictionary.f,RLS.Error), grid on
    xlabel('Frequency in Hz'), ylabel('DOA - Mean Squared Error')
    applyAxisProperties(gca)
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3))
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),RLS.Avg(1),RLS.Avg(2),RLS.Avg(3),2,'Linewidth',4)
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

