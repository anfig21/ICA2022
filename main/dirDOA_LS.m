function DOA = dirDOA_LS(Data,Direct,Dict,plotFlag)
%DOA = dirDOA_LS(Data,Direct,Dict,plotFlag) Applies Least Squares
%solution to the Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - DOA        : DOA estimation and dictionary via LS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('dirDOA_LS Error: Not enough input parameters.'), end

%% MAIN CODE
% Least-Squares solution across frequency
DOA.Est = nan(3,length(Dict.f));
for ii = 1:length(Dict.f)
    x = pinv(Dict.Plane.H(:,:,ii))*Direct.InnSph.H(Data.f==Dict.f(ii),:).';
    
    [~,Idx] = max(abs(x));
    DOA.Est(:,ii) = -Dict.Plane.uk(:,Idx);
end

DOA.Error = vecnorm(Direct.TrueDOA.'-DOA.Est);

DOA.Avg = mean(DOA.Est,2);
DOA.Avg = DOA.Avg/vecnorm(DOA.Avg);

if plotFlag
    % Mean Squared Error
    figure, plot(Dict.f,DOA.Error), grid on
    xlabel('Frequency in Hz'), ylabel('DOA - Mean Squared Error'), ylim([0 2])
    applyAxisProperties(gca)
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),DOA.Avg(1),DOA.Avg(2),DOA.Avg(3),2,'Linewidth',4)
    axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

disp('Direct sound: DOA - LS... OK')

end

