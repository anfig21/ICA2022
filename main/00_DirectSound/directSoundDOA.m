function Direct = directSoundDOA(Data,Direct,N,f,DOAMethod,plotFlag)
%Direct = directSoundDOA(Data,Direct,N,f,DOAMethod,plotFlag)
%Estimates the Direction-of-Arrival (DOA) of the direct sound using several
%optimisation methods.
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
if nargin < 5, error('directSoundDOA Error: Not enough input parameters.'), end
if N <= 0, error('directSoundDOA Error: The number of plane waves must be a positive integer.'), end
if nargin < 6, plotFlag = false; end

%% MAIN CODE
% True DOA
Direct.TrueDOA = Data.Source.pos-Data.Sph.R0;
Direct.TrueDOA = Direct.TrueDOA/vecnorm(Direct.TrueDOA);

% Dictionary of plane waves
Dict.f = f;
Dict.N = N;                 % Number of plane waves

[Dict.H,Dict.uk] = dictionary(Data.c,Dict.f,Data.InnSph.pos',Dict.N);
disp('Plane Wave Dictionary... OK')

% Optimisation problem
switch DOAMethod
    case 'LS'
        % DOA Estimation via Least-Squares
        Direct.DOA = dirDOA_LS(Data,Direct,Dict);
    case 'RLS'
        % DOA Estimation via Regularised Least-Squares
        Direct.DOA = dirDOA_RLS(Data,Direct,Dict);
    case 'CS'
        % DOA Estimation via Compressive Sensing
        Direct.DOA = dirDOA_CS(Data,Direct,Dict);
    case 'SOMP'
        % DOA Estimation via SOMP
        Direct.DOA = dirDOA_SOMP(Data,Direct,Dict);
    otherwise
        error('directSoundDOA Error: DOA estimation method not valid.')
end

% Mode & Error
Direct.DOA.Error = rad2deg(vecnorm(Direct.TrueDOA.'-Direct.DOA.Est));
Direct.DOA.Mode = mode(Direct.DOA.Est,2);

%% PLOT
if plotFlag
    if length(Dict.f) > 1
        % Mean Squared Error
        figure, plot(Dict.f,Direct.DOA.Error,'o-'), grid on
        xlabel('Frequency in Hz'), ylabel('DOA Error in Degrees'), ylim([0 180])
        applyAxisProperties(gca)
    end
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),...
        Direct.DOA.Mode(1),Direct.DOA.Mode(2),Direct.DOA.Mode(3),1.8,'Linewidth',4)
    
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','DOA Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

