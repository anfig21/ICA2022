function CS = directCS(Data,Direct,Dictionary,plotFlag)
%CS = directCS(Data,Direct,Dictionary,plotFlag) Applies Compressive Sensing
%to the Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dictionary: dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - CS        : DOA estimation and dictionary via CS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('directCS Error: Not enough input parameters.'), end

%% MAIN CODE
CS.Nnorm = mean(Direct.InnSph.Nnorm);

CS.Est = nan(3,length(Dictionary.f));

c = waitbar(0,'Loading...0\%','Name','CVX across frequencies...');
for ii = 1:length(Dictionary.f)
    Hii = squeeze(Dictionary.Plane.H(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dictionary.f(ii),:).';
    Nnorm = 1.1*CS.Nnorm;
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
        variable x(Dictionary.Plane.N) complex;
        minimize norm(x,1);
        subject to
            norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    [~,Idx] = max(abs(x));
    CS.Est(:,ii) = Dictionary.Plane.uk(:,Idx);
    
    waitbar(ii/length(Dictionary.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dictionary.f),2)),"\,\%"));
end
delete(c)

CS.Error = vecnorm(Direct.DOA.'-CS.Est);

CS.Avg = mean(CS.Est,2);
CS.Avg = CS.Avg/vecnorm(CS.Avg);

if plotFlag
    % Mean Squared Error
    figure, plot(Dictionary.f,CS.Error), grid on
    xlabel('Frequency in Hz'), ylabel('DOA - Mean Squared Error'), ylim([0 2])
    applyAxisProperties(gca)
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3))
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),-CS.Avg(1),-CS.Avg(2),-CS.Avg(3),2,'Linewidth',4)
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

