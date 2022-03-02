function LS = directLS(Data,Direct,Dictionary,plotFlag)
%LS = directLS(Data,Direct,Dictionary,plotFlag) Applies Least Squares
%solution to the Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dictionary: dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - LS        : DOA estimation and dictionary via LS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('directLS Error: Not enough input parameters.'), end

%% MAIN CODE
% Least-Squares solution across frequency
LS.Est = nan(3,length(Dictionary.f));
for ii = 1:length(Dictionary.f)
    LS.x = pinv(Dictionary.CA(:,:,ii))*Data.H(Data.f==Dictionary.f(ii),:).';
    
    [~,LS.Idx] = max(abs(LS.x));
    LS.Est(:,ii) = Dictionary.uk(:,LS.Idx);
end

LS.Error = vecnorm(Direct.DOA.'-LS.Est);

LS.Avg = mean(LS.Est,2);
LS.Avg = LS.Avg/vecnorm(LS.Avg);

if plotFlag
    % Mean Squared Error
    figure, plot(Dictionary.f,LS.Error), grid on
    xlabel('Frequency in Hz'), ylabel('DOA - Mean Squared Error')
    applyAxisProperties(gca)
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3))
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),LS.Avg(1),LS.Avg(2),LS.Avg(3),2,'Linewidth',4)
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

