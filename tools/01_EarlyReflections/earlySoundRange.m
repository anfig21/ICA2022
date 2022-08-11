function Early = earlySoundRange(Data,Early,res,rMinMax,f,RangeMethod,plotFlag)
%Early = earlySoundRange(Data,Early,res,rMinMax,f,RangeMethod,plotFlag)
%Estimates the Range of the virtual sources for the early reflections using
%several optimisation methods.
%   Input:
%       - Data          : raw data. Structure
%       - Early         : early reflections. Structure
%       - res           : spatial resolution in m. Double
%       - rMinMax       : minimum and maximum distances to the source. 2 x 1
%       - f             : frequency span. 1 x Nf
%       - RangeMethod   : specifies the regularisation method:
%                           'RLS' - Regularised Least Squares
%                           'CS' - Compressive Sensing
%                           'TV' - Total Variation
%       - plotFlag      : 'true' to plot Range estimation
%                         'false' (Default value)
%   Output:
%       - Early         : Range Estimation for early reflections. Structure
%
% Author: Antonio Figueroa Durán
% Date: July 2022

%% ERROR HANDLING
if nargin < 5, error('earlySoundRange Error: Not enough input parameters.'), end
if res <= 0, error('earlySoundRange Error: The spatial resolution must be a positive double.'), end
if numel(rMinMax) ~=2 || rMinMax(1) > rMinMax(2), error('earlySoundRange Error: rMinMax must be a vector of double in ascending order.'), end
if nargin < 6, plotFlag = false; end

%% MAIN CODE
% Dictionary of point sources
Dict.f = f;

[Dict.H,Dict.r,Dict.N] = dictionaryRange(Data.c,Dict.f,Data.InnSph.pos.',...
    Data.Sph.R0.',Early.DOA.Est.',rMinMax,res);
disp('Spherical Wave Dictionary... OK')

% Optimisation Problem
switch RangeMethod
    case 'RLS'
        % Range Estimation via Regularised Least-Squares
        Early.Range = dirRange_RLS(Data,Early,Dict);
    case 'CS'
        % Range Estimation via CS
        Early.Range = dirRange_CS(Data,Early,Dict);
    case 'TV'
        % Range Estimation via TV
        Early.Range = dirRange_TV(Data,Early,Dict);
    otherwise
        error('earlySoundRange Error: Range estimation method not valid.')
end

% Mode
Early.Range.Mode = mode(Early.Range.Est,2);

%% PLOT
if plotFlag    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    scatter3(Dict.r(1,:),Dict.r(2,:),Dict.r(3,:))
    scatter3(Early.Range.Mode(1),Early.Range.Mode(2),Early.Range.Mode(3),150,'filled')
    
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source','Candidates','Range Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

