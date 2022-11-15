function PWE = planeWaveExpansionTime_MiniBatching(Data,N,t,R,plotFlag)
%PWE = planeWaveExpansionTime_MiniBatching(Data,PWE,N,f,plotFlag) 
%   Input:
%       - Data      : raw data. Structure
%       - N         : number of plane waves. Integer
%       - t         : time vector. 1 x Nt
%       - R         : reconstruction points. 3 x Nr
%       - plotFlag  : 'true' to plot PW reconstruction
%                     'false' (Default value)
%   Output:
%       - PWE       : plane wave expansion in time domain. Structure
%
% Author: Antonio Figueroa Durán
% Date: July 2022

%% ERROR HANDLING
if nargin < 4, error('planeWaveExpansionTime_MiniBatching Error: Not enough input parameters.'), end
if N <= 0, error('planeWaveExpansionTime_MiniBatching Error: The number of plane waves must be a positive integer.'), end
if nargin < 5, plotFlag = false; end

%% MAIN CODE
% Coefficient estimation
PWE.Nt = length(t);
M = size(Data.InnSph.pos,1);
Mbatch = 7;

% Minibatching - random sampling
Idx = randperm(M);
Nbatches = ceil(M/Mbatch);
x = nan(Nt*N,Nbatches);

for ii = 1:Nbatches
    if ii == Nbatches
        micIdx = Idx(1+(ii-1)*Mbatch:end);
    else
        micIdx = Idx(1+(ii-1)*Mbatch:Mbatch+(ii-1)*Mbatch);
    end
    p = Data.InnSph.h(ismember(Data.t,t),micIdx);
    [H,~] = dictionaryPWTime(Data.c,t,Data.Fs,Data.InnSph.pos(micIdx,:)',Data.Sph.R0.',N);
    [x(:,ii),~] = reguLeastSquares(H,p(:));
    
    disp(strcat("Batch number... ",string(ii),"/",string(Nbatches)))
end
clear H

PWE.x = mean(x,2);

% RIR Reconstruction
PWE.h = PWTimeReconstruction(Data,PWE.x,t,N,R);

%% PLOT
% Validation vs reference line
if plotFlag
    if Nr == 1
        Mic = find(ismember(Data.Ref.pos,R','rows'));
        
        % Mic position
        setupPlot(Data,true,false,false);
        scatter3(R(1),R(2),R(3),150,'filled')

        % Reconstruction: impulse response
        figure, hold on
%         plot(t*1e3,Data.Ref.h(ismember(Data.t,t),Mic)/max(abs(Data.Ref.h(ismember(Data.t,t),Mic))))
%         plot(t*1e3,PWE.h/max(abs(Data.Ref.h(ismember(Data.t,t),Mic))))
        plot(t*1e3,Data.Ref.h(ismember(Data.t,t),Mic))
        plot(t*1e3,PWE.h)
        xlabel('Time in ms'), ylabel('Normalised RIR'), grid on
        legend('True','Reconstruction')
        applyAxisProperties(gca)
        applyLegendProperties(gcf)
    else
        % Reconstruction: reference line RIR
        figure
        s = surf(Data.Ref.pos(:,1),t*1e3,PWE.h);
        set(s,'edgecolor','none')
        xlabel('x in m'), ylabel('Time in ms')
        colormap hot
        view(2)
        c = colorbar;
%         caxis([-0.04 0.04])
        applyColorbarProperties(c,'Room Impulse Response in Pa/V')
        applyAxisProperties(gca)
    end
end
end

