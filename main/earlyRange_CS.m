function Range = earlyRange_CS(Data,Early,Dict,plotFlag)
%Range = earlyRange_CS(Data,Early,Dict,plotFlag) Applies Compressive
%Sensing to the early reflections to obtain the position of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early reflections. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - Range     : Range estimation via via CS. Structure
%
% Author: Antonio Figueroa Durán
% Date: April 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('earlyRange_CS Error: Not enough input parameters.'), end

%% MAIN CODE
%Range.Est = nan(3,length(Dict.f));
Nnorm = 1.1*Early.InnSph.NnormLcurve;

c = waitbar(0,'Loading...0\%','Name','earlyRange_CS: CVX across frequencies...');
cc = 0;
for rr = 1:Early.R
    for ii = 1:length(Dict.f)
        Hii = squeeze(Dict.SphEarly.H{rr}(:,:,ii));
        pii = Early.InnSph.H(Data.f==Dict.f(ii),:).';
        
        % CVX Formulation
        cvx_begin quiet
        cvx_precision high
        variable x(Dict.SphEarly.N{rr}) complex;
        minimize norm(x,1);
        subject to
        norm((Hii*x-pii),2) <= Nnorm;
        cvx_end
        
        % Coefficients plot
        figure
        stem(vecnorm(Dict.SphEarly.r{rr}-Data.Sph.R0.'),abs(x)), grid on
        xlabel('r in m'), ylabel('Coefficients')
        applyAxisProperties(gca)
        
        [~,Idx] = max(abs(x));
        Range.Est{rr}(:,ii) = Dict.SphEarly.r{rr}(:,Idx);
        
        cc = cc + 1;
        waitbar(cc/(length(Dict.f)*Early.R),c,strcat("Loading... ",...
            string(round(100*cc/(length(Dict.f)*Early.R),2)),"\,\%"));
    end
    
    Range.Avg{rr} = mode(Range.Est{rr},2);
end
delete(c)

%% PLOT
if plotFlag
    aux = Early.DOA.WL.Est+Data.Sph.R0.';
    
    figure
    %scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3))
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3)), hold on
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(aux(1,:),aux(2,:),aux(3,:),100,'filled','MarkerEdgeColor','k')
    for rr = 1:Early.R, scatter3(Range.Avg{rr}(1),Range.Avg{rr}(2),Range.Avg{rr}(3),170,'filled'), end
    for rr = 1:Early.R, scatter3(Dict.SphEarly.r{rr}(1,:),Dict.SphEarly.r{rr}(2,:),Dict.SphEarly.r{rr}(3,:)), end
    drawRoom(Data.D(1),Data.D(2),Data.D(3))
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Spherical Array','Source','DOAs','Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

disp('Early reflections: RANGE - CS... OK')

end

