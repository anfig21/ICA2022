function Range = dirRange_TV(Data,Direct,Dict,plotFlag)
%Range = dirRange_TV(Data,Direct,Dict,plotFlag) Applies Total Variation to
%the Direct Sound Field to obtain the position of the source.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - Range        : Range estimation via via TV. Structure
%
% Author: Antonio Figueroa Durán
% Date: May 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('dirRange_TV Error: Not enough input parameters.'), end

%% MAIN CODE
Range.Est = nan(3,length(Dict.f));
NoiseMargin = 10;           % dB

% Stencil
% st = 3; mask = [1; -1; 0];        % First order
st = 3; mask = [1; -2; 1];        % Second order
% st = 3; mask = [0.5; -1; 0.5];    % Weighted Second order

% st = 5; mask = [0.5; 1; -3; 1; 0.5];  % Fourth order

Mask = padarray(mask, Dict.Sph.N-st,'post');
M = circshift(Mask(:),Dict.Sph.N-(st-1)/2);
D = zeros(Dict.Sph.N);
for ii=1:Dict.Sph.N
    D(ii,:)=circshift(M.',[0 ii-1]);
end

c = waitbar(0,'Loading...0\%','Name','dirRange_TV: CVX across frequencies...');
for ii = 1:length(Dict.f)
    Nnorm = 10^(NoiseMargin/20)*Direct.InnSph.Nnorm(Data.f==Dict.f(ii));
    Hii = squeeze(Dict.Sph.H(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dict.f(ii),:).';
%     pii = unwrap(angle(pii));
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
    variable x(Dict.Sph.N) complex;
    minimize norm(D*x,1);
    subject to
    norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    [~,Idx] = max(abs(x));
    Range.Est(:,ii) = Dict.Sph.r(:,Idx);
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

Range.Error = vecnorm(Data.Source.pos.'-Range.Est);
disp(Range.Error)

Range.Avg = mode(Range.Est,2);

%% PLOT
if plotFlag
    figure
    stem(abs(x)), grid on
    xlabel('Wave Index'), ylabel('Coefficients Amplitude')
    applyAxisProperties(gca)
    
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(Dict.Sph.r(1,:),Dict.Sph.r(2,:),Dict.Sph.r(3,:))
    scatter3(Range.Est(1,:),Range.Est(2,:),Range.Est(3,:),100,'filled')
    scatter3(Range.Avg(1),Range.Avg(2),Range.Avg(3),170,'filled')
    axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Atoms','Estimation','Mode')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

disp('Direct sound: RANGE - TV... OK')

end

