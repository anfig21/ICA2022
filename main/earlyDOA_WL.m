function DOA = earlyDOA_WL(Data,Early,Dict,plotFlag)
%DOA = earlyDOA_WL(Data,Early,Dict,plotFlag) Applies Weighted LASSO to the
%early reflections to obtain a reflection map.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early reflections. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - DOA       : DOA estimation and dictionary via WL. Structure
%
% Author: Antonio Figueroa Durán
% Date: April 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('earlyDOA_WL Error: Not enough input parameters.'), end

%% MAIN CODE
DOA.x = nan(Dict.Plane.N,length(Dict.f));
NoiseMargin = 10;           % dB
mu = 0.01;
weights = ones(Dict.Plane.N,1);

c = waitbar(0,'Loading...0\%','Name','earlyDOA_WL: CVX across frequencies...');
for ii = 1:length(Dict.f)
    Nnorm = 10^(NoiseMargin/20)*Early.InnSph.Nnorm(Data.f==Dict.f(ii));
    Hii = squeeze(Dict.Plane.H(:,:,ii));
    pii = Early.InnSph.H(Data.f==Dict.f(ii),:).';
    
    for jj = 1:4
        % CVX Formulation
        cvx_begin quiet
        cvx_precision high
        variable x(Dict.Plane.N) complex
        minimize norm(weights.*x,1)
        subject to
        norm(Hii*x-pii,2) <= Nnorm;
        cvx_end
        
        % Weights
        weights = 1./(abs(x)+mu);
    end
    
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
    if length(Dict.f) > 1
        figure
        % s = surf(Dict.f,1:Dict.Plane.N,abs(DOA.x)./max(abs(DOA.x)));
        s = surf(Dict.f,1:Dict.Plane.N,abs(DOA.x));
        set(s,'edgecolor','none'), view(2)
        colormap hot
        xlabel('Frequency in Hz'), ylabel('Wave Index')
        applyAxisProperties(gca)
    end
    
    figure
    stem(abs(x)), grid on
    xlabel('Wave Index'), ylabel('Coefficients Amplitude')
    applyAxisProperties(gca)
    
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

disp('Early reflections: DOA - WL... OK')

end

