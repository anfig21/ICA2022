function DOA = earlyDOA_CS(Data,Early,Dict,plotFlag)
%DOA = earlyDOA_CS(Data,Early,Dict,plotFlag) Applies Compressive
%Sensing to the early reflections to obtain a reflection map.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early reflections. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - DOA        : DOA estimation and dictionary via CS. Structure
%
% Author: Antonio Figueroa Durán
% Date: April 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('earlyDOA_CS Error: Not enough input parameters.'), end

%% MAIN CODE
DOA.x = nan(Dict.Plane.N,length(Dict.f));
NoiseMargin = 10;           % dB

c = waitbar(0,'Loading...0\%','Name','earlyDOA_CS: CVX across frequencies...');
for ii = 1:length(Dict.f)
    Nnorm = 10^(NoiseMargin/20)*Early.InnSph.Nnorm(Data.f==Dict.f(ii));
    Hii = squeeze(Dict.Plane.H(:,:,ii));
    pii = Early.InnSph.H(Data.f==Dict.f(ii),:).';
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
    variable x(Dict.Plane.N) complex;
    minimize norm(x,1);
    subject to
    norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    DOA.x(:,ii) = x;
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

x = sum(abs(DOA.x),2);
uk = Dict.Plane.uk+Data.Sph.R0.';
[~,Idx] = maxk(abs(x),Early.R);
DOA.Est = -Dict.Plane.uk(:,Idx);
ukEst = -Dict.Plane.uk(:,Idx)+Data.Sph.R0.';

%% PLOT
if plotFlag
    if length(Dict.f) == 1
        % Coefficients
        figure
        stem(abs(DOA.x)), grid on
        xlabel('Wave Index'), ylabel('Coefficients Amplitude')
        applyAxisProperties(gca)
    end
    
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(uk(1,:),uk(2,:),uk(3,:),'MarkerEdgeColor', uint8([200 200 200]))
    scatter3(ukEst(1,:),ukEst(2,:),ukEst(3,:),100,'filled','MarkerEdgeColor','k')
    axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Dictionary Atoms','DOA Estimations')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

disp('Early Reflections: DOA - CS... OK')

end

