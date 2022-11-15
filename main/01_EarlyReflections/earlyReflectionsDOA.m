function Early = earlyReflectionsDOA(Data,Early,N,f,DOAMethod,plotFlag)
%Direct = earlyReflectionsDOA(Data,Early,N,f,DOAMethod,plotFlag)
%Estimates the Direction-of-Arrival (DOA) of the early reflections using
%several optimisation methods.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure
%       - N         : number of plane waves. Integer
%       - f         : frequency span. 1 x Nf
%       - DOAMethod : specifies the optimisation method:
%                       'LS' - Least Squares
%                       'RLS' - Regularised Least Squares
%                       'CS' - Compressive Sensing
%                       'SOMP' - Simultaneous Orthogonal Matching Pursuit
%       - plotFlag  : 'true' to plot DOA estimation
%                     'false' (Default value)
%   Output:
%       - Direct    : direct sound field with DOA Estimation. Structure
%
% Author: Antonio Figueroa Durán
% Date: June 2022

%% ERROR HANDLING
if nargin < 5, error('earlyReflectionsDOA Error: Not enough input parameters.'), end
if N <= 0, error('earlyReflectionsDOA Error: The number of plane waves must be a positive integer.'), end
if nargin < 6, plotFlag = false; end

%% MAIN CODE
% Dictionary of plane waves
Dict.f = f;
Dict.N = N;                 % Number of plane waves

[Dict.H,Dict.uk] = dictionary(Data.c,Dict.f,Data.InnSph.pos',Dict.N);
disp('Plane Wave Dictionary... OK')

% Optimisation problem
switch DOAMethod
    case 'RLS'
        % DOA Estimation via Regularised Least-Squares
        Early.DOA = earlyDOA_RLS(Data,Early,Dict);
    case 'CS'
        % DOA Estimation via Compressive Sensing
        Early.DOA = earlyDOA_CS(Data,Early,Dict);
    case 'EN'
        % DOA Estimation via Elastic Net
        Early.DOA = earlyDOA_EN(Data,Early,Dict);
    case 'WL'
        % DOA Estimation via Weighted LASSO
        Early.DOA = earlyDOA_WL(Data,Early,Dict);
    otherwise
        error('directSoundDOA Error: DOA estimation method not valid.')
end

% Mode
[~,Idx] = maxk(Early.DOA.x,Early.R);
Early.DOA.Est = -Dict.uk(:,Idx);

%% PLOT
if plotFlag
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(Dict.uk(1,:)+Data.Sph.R0(1),Dict.uk(2,:)+Data.Sph.R0(2),...
        Dict.uk(3,:)+Data.Sph.R0(3),'MarkerEdgeColor', uint8([200 200 200]))
    scatter3(Early.DOA.Est(1,:)+Data.Sph.R0(1),Early.DOA.Est(2,:)+Data.Sph.R0(2),...
        Early.DOA.Est(3,:)+Data.Sph.R0(3),100,'filled','MarkerEdgeColor','k')
    
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Dictionary','DOA Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

