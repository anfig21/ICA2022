function Direct = directSoundRange(Data,Direct,res,rMinMax,f,RangeMethod,plotFlag)
%Direct = directSoundRange(Data,Direct,res,rMinMax,f,RangeMethod,plotFlag)
%Estimates the Range of the source for the direct sound using several
%optimisation methods.
%   Input:
%       - Data          : raw data. Structure
%       - Direct        : direct sound field. Structure
%       - res           : spatial resolution in m. Double
%       - rMinMax       : minimum and maximum distances to the source. 2 x 1
%       - f             : frequency span. 1 x Nf
%       - RangeMethod   : specifies the regularisation method:
%                           'RLS' - Regularised Least Squares
%                           'CS' - Compressive Sensing
%                           'EN' - Elastic Net
%                           'TV' - Total Variation
%       - plotFlag      : 'true' to plot Range estimation
%                         'false' (Default value)
%   Output:
%       - Direct        : direct sound field with Range Estimation. Structure
%
% Author: Antonio Figueroa Durán
% Date: June 2022

%% ERROR HANDLING
if nargin < 5, error('directSoundRange Error: Not enough input parameters.'), end
if res <= 0, error('directSoundRange Error: The spatial resolution must be a positive double.'), end
if numel(rMinMax) ~=2 || rMinMax(1) > rMinMax(2), error('directSoundRange Error: rMinMax must be a vector of double in ascending order.'), end
if nargin < 6, plotFlag = false; end

%% MAIN CODE
% Dictionary of point sources
Dict.f = f;

[Dict.H,Dict.r,Dict.N] = dictionaryRange(Data.c,Dict.f,Data.InnSph.pos.',...
    Data.Sph.R0.',Direct.DOA.Mode.',rMinMax,res);
disp('Spherical Wave Dictionary... OK')

% Optimisation Problem
switch RangeMethod
    case 'RLS'
        % Range Estimation via Regularised Least-Squares
        Direct.Range = dirRange_RLS(Data,Direct,Dict);
    case 'CS'
        % Range Estimation via CS
        Direct.Range = dirRange_CS(Data,Direct,Dict);
    case 'EN'
        % Range Estimation via EN
        Direct.Range = dirRange_EN(Data,Direct,Dict);
    case 'TV'
        % Range Estimation via TV
        Direct.Range = dirRange_TV(Data,Direct,Dict);
    otherwise
        error('directSoundRange Error: Range estimation method not valid.')
end

% Mode & Error
Direct.Range.Error = vecnorm(Data.Source.pos.'-Direct.Range.Est);
Direct.Range.Mode = mode(Direct.Range.Est,2);

%% PLOT
if plotFlag
    if length(Dict.f) > 1
        % Mean Squared Error
        figure, plot(Dict.f,Direct.Range.Error,'o-'), grid on
        xlabel('Frequency in Hz'), ylabel('Mean Squared Error in meters'), ylim([0 +inf])
        applyAxisProperties(gca)
    end
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(Dict.r(1,:),Dict.r(2,:),Dict.r(3,:))
    scatter3(Direct.Range.Mode(1),Direct.Range.Mode(2),Direct.Range.Mode(3),150,'filled')
    
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Candidates','Range Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

