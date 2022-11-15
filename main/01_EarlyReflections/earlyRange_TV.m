function Range = earlyRange_TV(Data,Early,Dict)
%Range = earlyRange_TV(Data,Direct,Dict) Applies Total Variation to
%the early reflections to obtain the position of the source
%   Input:
%       - Data          : raw data. Structure
%       - Early         : early reflections. Structure
%       - Dict          : dictionary of plane waves. Structure
%   Output:
%       - Range        : Range estimation via TV. Structure
%
% Author: Antonio Figueroa Durán
% Date: May 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('earlyRange_TV Error: Not enough input parameters.'), end

%% MAIN CODE
NoiseMargin = 10;           % dB

% Stencil
% st = 3; mask = [1; -1; 0];        % First order
st = 3; mask = [1; -2; 1];        % Second order
% st = 3; mask = [0.5; -1; 0.5];    % Weighted Second order

% st = 5; mask = [0.5; 1; -3; 1; 0.5];  % Fourth order

c = waitbar(0,'Loading...0\%','Name','earlyRange_TV: CVX across frequencies...');
cc = 0;
for rr = 1:Early.R
    Mask = padarray(mask, Dict.N{rr}-st,'post');
    M = circshift(Mask(:),Dict.N{rr}-(st-1)/2);
    D = zeros(Dict.N{rr});
    for ii=1:Dict.N{rr}
        D(ii,:)=circshift(M.',[0 ii-1]);
    end
    
    for ii = 1:length(Dict.f)
        Nnorm = 10^(NoiseMargin/20)*Early.InnSph.Nnorm(Data.f==Dict.f(ii));
        Hii = squeeze(Dict.H{rr}(:,:,ii));
        pii = Early.InnSph.H(Data.f==Dict.f(ii),:).';
        %     pii = unwrap(angle(pii));
        
        % CVX Formulation
        cvx_begin quiet
        cvx_precision high
        variable x(Dict.N{rr}) complex;
        minimize norm(D*x,1);
        subject to
        norm((Hii*x-pii),2) <= Nnorm;
        cvx_end
        
        % Coefficients plot
        figure
        stem(vecnorm(Dict.r{rr}-Data.Sph.R0.'),abs(x)), grid on
        xlabel('r in m'), ylabel('Coefficients')
        applyAxisProperties(gca)
        
        [~,Idx] = max(abs(x));
        Range.Est{rr}(:,ii) = Dict.r{rr}(:,Idx);
        
        cc = cc + 1;
        waitbar(cc/(length(Dict.f)*Early.R),c,strcat("Loading... ",...
            string(round(100*cc/(length(Dict.f)*Early.R),2)),"\,\%"));
    end
    Range.Avg{rr} = mode(Range.Est{rr},2);
    
end
delete(c)

%% PLOT
if plotFlag    
    figure
    %scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3))
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3)), hold on
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    for rr = 1:Early.R, scatter3(Range.Avg{rr}(1),Range.Avg{rr}(2),Range.Avg{rr}(3),170,'filled'), end
    for rr = 1:Early.R, scatter3(Dict.SphEarly.r{rr}(1,:),Dict.SphEarly.r{rr}(2,:),Dict.SphEarly.r{rr}(3,:)), end
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Spherical Array','Source','Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

disp('Early reflections: RANGE - TV... OK')

end

