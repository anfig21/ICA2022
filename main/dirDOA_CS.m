function DOA = dirDOA_CS(Data,Direct,Dict,plotFlag)
%DOA = dirDOA_CS(Data,Direct,Dict,plotFlag) Applies Compressive
%Sensing to the Direct Sound Field to obtain the DOA of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - DOA        : DOA estimation and dictionary via CS. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('dirDOA_CS Error: Not enough input parameters.'), end

%% MAIN CODE
DOA.Est = nan(3,length(Dict.f));
Nnorm = 1.1*Direct.InnSph.NnormLcurve;

c = waitbar(0,'Loading...0\%','Name','dirDOA_CS: CVX across frequencies...');
for ii = 1:length(Dict.f)
    Hii = squeeze(Dict.Plane.H(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dict.f(ii),:).';
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
        variable x(Dict.Plane.N) complex;
        minimize norm(x,1);
        subject to
            norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    [~,Idx] = max(abs(x));
    
    DOA.Est(:,ii) = -Dict.Plane.uk(:,Idx);
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

DOA.Error = rad2deg(vecnorm(Direct.TrueDOA.'-DOA.Est));

DOA.Avg = mode(DOA.Est,2);
DOA.Avg = DOA.Avg/vecnorm(DOA.Avg);

%% PLOT
if plotFlag
    % Mean Squared Error
    figure, plot(Dict.f,DOA.Error), grid on
    xlabel('Frequency in Hz'), ylabel('DOA Error in Degrees'), ylim([0 2])
    applyAxisProperties(gca)
    
    figure
    stem(abs(x)), grid on
    xlabel('Wave Index'), ylabel('Coefficients Amplitude')
    applyAxisProperties(gca)
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    %scatter3(DOA.Est(1,:)+Data.Sph.R0(1),DOA.Est(2,:)+Data.Sph.R0(2),DOA.Est(3,:)+Data.Sph.R0(3))
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),DOA.Avg(1),DOA.Avg(2),DOA.Avg(3),2,'Linewidth',4)
    axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

disp('Direct sound: DOA - CS... OK')

end

