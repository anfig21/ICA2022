function PWE = planeWaveExpansionTime(Data,N,t,Q,R,plotFlag)
%PWE = planeWaveExpansionTime(Data,N,t,Q,R,plotFlag) Decompose the sound
%field into plane waves in time domain, where windows are taken applying
%Overlap-Save.
%   Input:
%       - Data      : raw data. Structure
%       - N         : number of plane waves. Integer
%       - t         : time vector. 1 x Nt
%       - Q         : downsampling factor. Scalar
%       - R         : reconstruction points. 3 x Nr
%       - plotFlag  : 'true' to plot PW reconstruction
%                     'false' (Default value)
%   Output:
%       - PWE       : plane wave expansion in time domain. Structure
%
% Author: Antonio Figueroa Durán
% Date: July 2022

%% ERROR HANDLING
if nargin < 5, error('planeWaveExpansionTime Error: Not enough input parameters.'), end
if N <= 0, error('planeWaveExpansionTime Error: The number of plane waves must be a positive integer.'), end
if nargin < 6, plotFlag = false; end

%% MAIN CODE

% Downsampling
PWE.Fs = ceil(Data.Fs/Q);
PWE.t = t(1:Q:end);

% Low-pass Filter
Fc = 4e3;
lpFilt = designfilt('lowpassfir','PassbandFrequency',Fc, ...
            'StopbandFrequency',Fc*1.3,'PassbandRipple',0.5, ...
            'StopbandAttenuation',65,'DesignMethod','kaiserwin',...
            'SampleRate',PWE.Fs);
p_pre = filtfilt(lpFilt,Data.InnSph.h(ismember(Data.t,PWE.t),:));
        
% Array data:
%  d = 50 cm
%  DeltaT = d/c = 1.4461 ms
%  NDeltaT = DeltaT*Fs = 69.4 == 70 samples

% Pre zero-padding
M = size(Data.InnSph.pos,1);
overlapp = 15;      % samples
Ny = size(p_pre,1);

p = vertcat(zeros(overlapp-1,M),p_pre);
PWE.Nt = size(p,1);

% Overlapp-Save Parameters
Nblock = 41;
L = Nblock-overlapp+1;      % Convolution order

nFrames = ceil(PWE.Nt/L);

% Zero-pad at the end of the signal (last frame)
p = vertcat(p,zeros(nFrames*L-Ny,M));

% Time vector
t_padded = horzcat(zeros(1,overlapp-1),Data.t);
t_padded = t_padded(find(t_padded==PWE.t(1),1,'last')-(overlapp-1)+(1:size(p,1))-1);

% Coefficient estimation
PWE.x = nan(PWE.Nt,N);

for ii = 0:nFrames-1
    pm = p((1:L+overlapp-1) + ii*L,:);
    tm = t_padded((1:L+overlapp-1) + ii*L);
    
    [H,~] = dictionaryPWTime(Data.c,tm,PWE.Fs,Data.InnSph.pos',Data.Sph.R0.',N);

    [xm,~] = reguLeastSquares(H,pm(:),'lcurve',true);

    xm = reshape(xm,Nblock,N);
    PWE.x(ii*L+1:(ii+1)*L,:) = xm(overlapp:overlapp+L-1,:);
    
    disp(strcat("Frame number... ",string(ii+1),"/",string(nFrames)))
end
clear H

PWE.x = PWE.x(1:Ny,:);
PWE.x = PWE.x(:);

% RIR Reconstruction
PWE.h = PWTimeReconstruction(Data,PWE.x,PWE.t,N,R);

%% PLOT
% Validation vs reference line
if plotFlag
    if M == 1
        Mic = find(ismember(Data.Ref.pos,R','rows'));
        
        % Mic position
        setupPlot(Data,true,false,false);
        scatter3(R(1),R(2),R(3),150,'filled')

        % Reconstruction: impulse response
        figure, hold on
        plot(PWE.t*1e3,Data.Ref.h(ismember(Data.t,PWE.t),Mic))
        plot(PWE.t*1e3,PWE.h)
        xlabel('Time in ms'), ylabel('Normalised RIR'), grid on
        legend('True','Reconstruction')
        applyAxisProperties(gca)
        applyLegendProperties(gcf)
    else
        % Reconstruction: reference line RIR
        figure, hold on
        s = pcolor(Data.Ref.pos(:,1),PWE.t*1e3,PWE.h);
        set(s,'edgecolor','none')
        xlabel('x in m'), ylabel('Time in ms')
        colormap hot
        c = colorbar;
        caxis([-0.04 0.04])
        xline(min(Data.InnSph.pos(:,1))), xline(max(Data.InnSph.pos(:,1)))
        applyColorbarProperties(c,'Room Impulse Response in Pa/V')
        applyAxisProperties(gca)
    end
end
end

