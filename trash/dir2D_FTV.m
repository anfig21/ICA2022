function TV = dir2D_FTV(Data,Direct,Dict,plotFlag)
%TV = dir2D_FTV(Data,Direct,Dict,plotFlag) Applies Fused Total
%Generalised Variation with a 3x3 stencil to the Direct Sound Field. The
%dictionary is comprised of point sources along a rectangular grid parallel
%to the XY-plan through the center of the spherical array.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - Dict      : dictionary of plane waves. Structure
%       - plotFlag  : 'true' to plot setup & DOA estimation
%                     'false' to avoid plotting. Default value
%   Output:
%       - TV        : Estimation and dictionary via FTV. Structure
%
% Author: Antonio Figueroa Durán
% Date: March 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('dir2D_FTV Error: Not enough input parameters.'), end

%% MAIN CODE
TV.Est = nan(3,length(Dict.f));

Nx = Dict.Grid.Nx;
Ny = Dict.Grid.Ny;

st = 3;     % Stencil dimension
mask = [0.5,1,0.5; 1,-6,1; 0.5,1,0.5];
Mask = padarray(mask, [Ny Nx]-st,'post');
M = circshift(Mask(:),Ny*Nx-((st-1)*0.5*Ny+(st-1)/2));
D = zeros(Ny*Nx);
for ii = 1:Ny*Nx
    D(ii,:) = circshift(M.',[0 ii-1]);
end
mu = 1;
Dcs = [D ; mu*eye(Nx*Ny)];

c = waitbar(0,'Loading...0\%','Name','CVX across frequencies...');
for ii = 1:length(Dict.f)
    Hii = squeeze(Dict.Grid.H(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dict.f(ii),:).';
    Nnorm = 1.1*Direct.InnSph.Nnorm(Data.f==Dict.f(ii));
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
    variable x(Dict.Grid.N) complex;
    minimize norm(Dcs*x,1);
    subject to
        norm((Hii*x-pii)) <= Nnorm;
    cvx_end
    
    [~,Idx] = max(abs(x));
    TV.Est(:,ii) = Dict.Grid.r(:,Idx);
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

TV.Error = vecnorm(Direct.DOA.'-TV.Est);

TV.Avg = mode(TV.Est,2);

%% PLOT
if plotFlag
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(Dict.Grid.r(1,:),Dict.Grid.r(2,:),Dict.Grid.r(3,:))
    scatter3(TV.Est(1,:),TV.Est(2,:),TV.Est(3,:),100,'filled')
    scatter3(TV.Avg(1),TV.Avg(2),TV.Avg(3),170,'filled')
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Atoms','Estimation','Mode')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

