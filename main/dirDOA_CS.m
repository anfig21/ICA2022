function CS = dirDOA_CS(Data,Direct,Dict,plotFlag)
%CS = dirDOA_CS(Data,Direct,Dict,plotFlag) Applies Compressive
%Sensing to the Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
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
elseif nargin < 3, error('dirDOA_CS Error: Not enough input parameters.'), end

%% MAIN CODE
CS.DOA.Est = nan(3,length(Dict.f));

c = waitbar(0,'Loading...0\%','Name','CVX across frequencies...');
for ii = 1:length(Dict.f)
    Hii = squeeze(Dict.Plane.H(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dict.f(ii),:).';
    Nnorm = 1.1*Direct.InnSph.Nnorm(Data.f==Dict.f(ii));
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
        variable x(Dict.Plane.N) complex;
        minimize norm(x,1);
        subject to
            norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    [~,Idx] = max(abs(x));
    CS.DOA.Est(:,ii) = Dict.Plane.uk(:,Idx);
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

CS.DOA.Error = vecnorm(Direct.DOA.'-CS.DOA.Est);

CS.DOA.Avg = mean(CS.DOA.Est,2);
CS.DOA.Avg = CS.DOA.Avg/vecnorm(CS.DOA.Avg);

%% PLOT
if plotFlag
    % Mean Squared Error
    figure, plot(Dict.f,CS.DOA.Error), grid on
    xlabel('Frequency in Hz'), ylabel('DOA - Mean Squared Error'), ylim([0 2])
    applyAxisProperties(gca)
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),-CS.DOA.Avg(1),-CS.DOA.Avg(2),-CS.DOA.Avg(3),2,'Linewidth',4)
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end
