%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       SOUND FIELD RECONSTRUCTION
%
% -------------------------------------------------------------------------
% E. Fernandez-Grande et al., "Reconstruction of room impulse responses
% over extended domains for navigable sound field reproduction", 2021
% -------------------------------------------------------------------------
%
% Antonio Figueroa Durán
% anfig@elektro.dtu.dk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all

addpath(genpath('../tools'))
% addpath(genpath('/Volumes/anfig/Data/room019/'))      % MacBook Pro
addpath(genpath('M:\Data'))                             % Windows 10

loadPlotParams

%% INITIAL PARAMETERS
Data.T = 5;
Data.Fs = 48e3;
Data.D = [6.266,9.357,2.977];       % Room dimensions
Data.Nsamples = Data.Fs*Data.T;
Data.f = Data.Fs/Data.Nsamples*(0:Data.Nsamples-1)/2;

Data.loudspeaker = 'rir_019_spk1.h5';
% Data.loudspeaker = 'rir_019_spk2.h5';

%% DATA ACQUISITION
Data = dataAcquisition(Data);

%% DATA HANDLING
Data.Ref.h = horzcat(Data.Line1.h,Data.Line2.h,Data.Line3.h);
Data.Ref.n = horzcat(Data.Line1.n,Data.Line2.n,Data.Line3.n);
Data.Ref.pos = vertcat(Data.Line1.pos,Data.Line2.pos,Data.Line3.pos);

% Inner Sphere (156 Microphones: samples 155-end)
Data.InnSph.h = Data.Sph.h(:,155:end);
Data.InnSph.pos = Data.Sph.pos(155:end,:);
Data.InnSph.p = Data.Sph.p(:,155:end);
Data.InnSph.M = size(Data.InnSph.pos,1);

% Frequency domain
Data.InnSph.H = fft(Data.InnSph.h,2*Data.Nsamples)/Data.Nsamples;
Data.InnSph.H = [Data.InnSph.H(1,:); 2*Data.InnSph.H(2:end/2,:)];
Data.Sph.N = fft(Data.Sph.n,2*Data.Nsamples)/Data.Nsamples;
Data.Sph.N = [Data.Sph.N(1,:); 2*Data.Sph.N(2:end/2,:)];
Data.InnSph.P = fft(Data.InnSph.p,2*Data.Nsamples)/Data.Nsamples;
Data.InnSph.P = [Data.InnSph.P(1,:); 2*Data.InnSph.P(2:end/2,:)];

% Noise norm
Data.Sph.Nnorm = mean(abs(Data.Sph.N),2);

%% SETUP PLOT
Plot.T = [5 25]*1e-3;               % Plot lenght
Plot.N = Data.Fs*Plot.T;

% Time vector
Plot.t = Plot.T(1):1/Data.Fs:Plot.T(2)-(1/Data.Fs);

% Data downsizing
Plot.Ref.h = Data.Ref.h(Plot.N(1):Plot.N(2)-1,:);
Plot.InnSph.h = Data.InnSph.h(Plot.N(1):Plot.N(2)-1,:);

% plotFreqResponse(Data)

% REFERENCE LINE RIR PLOT
% figure
% s = surf(Data.Ref.pos(:,1),Plot.t*1e3,Plot.Ref.h);
% set(s,'edgecolor','none')
% xlabel('x in m'), ylabel('Time in ms')
% colormap hot
% view(2)
% c = colorbar;
% caxis([-0.04 0.04])
% applyColorbarProperties(c,'Room Impulse Response in Pa/V')
% applyAxisProperties(gca)

%% ------------ DIRECT SOUND FIELD ------------
% Time window: 5-10 ms
Direct.T = 10*1e-3;
Direct.Nsamples = Data.Fs*Direct.T;

% Windowing - Hanning (hann) window
Direct.w = vertcat(repmat(hann(Direct.Nsamples),1,Data.InnSph.M),...
    zeros(Data.Nsamples-Direct.Nsamples,Data.InnSph.M));
Direct.InnSph.h = Direct.w.*Data.InnSph.h;
Direct.InnSph.n = Direct.w(:,size(Data.Sph.n,2)).*Data.Sph.n(1:2:end,:);

% Frequency Domain
Direct.InnSph.H = fft(Direct.InnSph.h,2*Data.Nsamples)/Data.Nsamples;
Direct.InnSph.H = [Direct.InnSph.H(1,:); 2*Direct.InnSph.H(2:end/2,:)];
Direct.InnSph.N = fft(Direct.InnSph.n,2*Data.Nsamples)/Data.Nsamples;
Direct.InnSph.N = [Direct.InnSph.N(1,:); 2*Direct.InnSph.N(2:end/2,:)];

% Noise norm
Direct.InnSph.Nnorm = mean(abs(Direct.InnSph.N),2);

% True DOA
Direct.DOA = Data.Sph.R0-Data.Source.pos;
Direct.DOA = Direct.DOA/vecnorm(Direct.DOA);

% Plot frequency response of mics [40, 59, 69]
% figure, plot(Data.f*1e-3,20*log10(abs(Direct.InnSph.H(:,[40 59 69])))), grid on
% xlabel('Frequency in kHz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
% legend('Mic 40','Mic 59','Mic 69')
% applyAxisProperties(gca)
% applyLegendProperties(gcf)

%% Dictionary of plane waves
% Dictionary
Dict.f = Data.f(200 <= Data.f & Data.f <= 300);
Dict.f = Dict.f(1:50:end);
% Dict.f = 2e2;                 % DOA estimation at 200 Hz
Dict.Plane.N = 1e3;             % Number of plane waves
Dict.Plane.K = 1;               % Number of sources (SOMP)

[Dict.Plane.H,Dict.Plane.uk] = dictionary(Dict.f,Data.InnSph.pos',Dict.Plane.N);

%% DOA Estimation via SOMP
% DOA = dirDOA_SOMP(Data,Direct,Dict,'true');

%% DOA Estimation via Least-Squares
% LS = dirDOA_LS(Data,Direct,Dict,'true');

%% DOA Estimation via Regularised Least-Squares
% RLS = dirDOA_RLS(Data,Direct,Dict,'true');
% Nnorm = 1.685*10-6 via L-curve method
Direct.InnSph.NnormLcurve = 1.685e-6;

%% DOA Estimation via Compressive Sensing
CS = dirDOA_CS(Data,Direct,Dict,'true');

%% Dictionary of spherical waves
Dict.Sph.Res = 5e-2;
Dict.Sph.rMinMax = [0.3 7];

[Dict.Sph.H,Dict.Sph.r,Dict.Sph.N] = dictionaryRange(Dict.f,Data.InnSph.pos.',...
    Data.Sph.R0.',CS.DOA.Avg.',Dict.Sph.rMinMax,Dict.Sph.Res);

%% Range Estimation via Regularised Least-Squares
% RLS = dirRange_RLS(Data,Direct,RLS,Dict,'true');

%% Range Estimation via CS
CS = dirRange_CS(Data,Direct,CS,Dict,'true');




