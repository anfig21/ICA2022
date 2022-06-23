function Plot = setupPlot(Data,flagS,flagH,flagT)
%setupPlot(Data,flagS,flagH,flagT) Plots the measurement setup, the
%frequency response and the RIR at the reference microphones.
%   Input:
%       - Data      : data structure.
%       - flagS     : plots setup
%                       'true'
%                       'false' (Default value)
%       - flagH     : plots frequency response
%                       'true'
%                       'false' (Default value)
%       - flagT     : plots RIR at reference line
%                       'true'
%                       'false' (Default value)
%
% Author: Antonio Figueroa Durán
% Date: June 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 1, error('setupPlot Error: Not enough input parameters.'), end
if nargin < 4, flagT = false; end
if nargin < 3, flagH = false; end
if nargin < 2, flagS = false; end

%% MAIN CODE
Plot.T = [5 25]*1e-3;       % Source near field
% Plot.T = [15 35]*1e-3;    % Source far field
Plot.N = Data.Fs*Plot.T;

% Time vector
Plot.t = Plot.T(1):1/Data.Fs:Plot.T(2)-(1/Data.Fs);

% Data downsizing
Plot.Ref.h = Data.Ref.h(Plot.N(1):Plot.N(2)-1,:);
Plot.InnSph.h = Data.InnSph.h(Plot.N(1):Plot.N(2)-1,:);

if flagS
    figure
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled'), hold on
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3))
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Source','Reference Line','Spherical Array')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

if flagH, plotFreqResponse(Data), end

% REFERENCE LINE RIR PLOT
if flagT
    figure
    s = surf(Data.Ref.pos(:,1),Plot.t*1e3,Plot.Ref.h);
    set(s,'edgecolor','none')
    xlabel('x in m'), ylabel('Time in ms')
    colormap hot
    view(2)
    c = colorbar;
    caxis([-0.04 0.04])
    applyColorbarProperties(c,'Room Impulse Response in Pa/V')
    applyAxisProperties(gca)
end

end

