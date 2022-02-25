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
%addpath(genpath('/Volumes/anfig/Data/room019/'))   % MacBook Pro
addpath(genpath('M:\Data'))                         % Windows 10

loadPlotParams

%% INITIAL PARAMETERS
Data.T = 5;
Data.Fs = 48e3;
Data.D = [6.266,9.357,2.977];       % Room dimensions
Data.Nsamples = Data.Fs*Data.T;
Data.F = [20 20e3];
Data.f = linspace(Data.F(1),Data.F(2),Data.Nsamples);

Plot.T = [5 25]*1e-3;               % Plot lenght
Plot.N = Data.Fs*Plot.T;

% Data.loudspeaker = 'rir_019_spk1.h5';
Data.loudspeaker = 'rir_019_spk2.h5';

%% DATA ACQUISITION
% Source
Data.Source.pos = h5read(Data.loudspeaker,'/source/position');

% Line 1
Data.Line1.pos = h5read(Data.loudspeaker,'/dataset_lin2/position');
Data.Line1.h = h5read(Data.loudspeaker,'/dataset_lin2/impulse_response');

% Line 2
Data.Line2.pos = h5read(Data.loudspeaker,'/dataset_lin1/position');
Data.Line2.h = h5read(Data.loudspeaker,'/dataset_lin1/impulse_response');

% Line 3
Data.Line3.pos = h5read(Data.loudspeaker,'/dataset_lin3/position');
Data.Line3.h = h5read(Data.loudspeaker,'/dataset_lin3/impulse_response');

% Sphere
Data.Sph.pos = h5read(Data.loudspeaker,'/dataset_sph1/position');
Data.Sph.h = h5read(Data.loudspeaker,'/dataset_sph1/impulse_response');
Data.Sph.R0 = mean(Data.Sph.pos,1);

%% DATA MERGING
Data.Ref.h = horzcat(Data.Line1.h,Data.Line2.h,Data.Line3.h);
Data.Ref.pos = vertcat(Data.Line1.pos,Data.Line2.pos,Data.Line3.pos);

% Inner Sphere (156 Microphones: samples 155-end)
Data.InnSph.h = Data.Sph.h(:,155:end);
Data.InnSph.pos = Data.Sph.pos(155:end,:);
Data.InnSph.M = size(Data.InnSph.pos,1);

%% SETUP PLOT
% Time vector
Plot.t = Plot.T(1):1/Data.Fs:Plot.T(2)-(1/Data.Fs);

% Data downsizing
Plot.Ref.h = Data.Ref.h(Plot.N(1):Plot.N(2)-1,:);
Plot.InnSph.h = Data.InnSph.h(Plot.N(1):Plot.N(2)-1,:);

% Line plot
% figure
% scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
% scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
% scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3))
% axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
% xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
% legend('Reference Line','Spherical Array','Source')
% applyAxisProperties(gca)
% applyLegendProperties(gca)

%% REFERENCE LINE RIR PLOT
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
Direct.T = [5 10]*1e-3;
Direct.N = Data.Fs*Direct.T;
Direct.Nsamples = Direct.N(2)-Direct.N(1);
Direct.f = linspace(Data.F(1),Data.F(2),Direct.Nsamples);

% Data Downsizing
Direct.InnSph.h = Data.InnSph.h(Direct.N(1):Direct.N(2)-1,:);

% Frequency Domain
Direct.H = fft(Direct.InnSph.h,2*Direct.Nsamples)/Direct.Nsamples;
Direct.H = [Direct.H(1,:); 2*Direct.H(2:end/2,:)];

%% DOA Estimation via SOMP
DOA = directSOMP(Data,Direct);

%% DOA Estimation via TDOA with only 2 mics (Mics 59 & 69)
TDOA.Method = 'PHAT';
TDOA.p = (-Direct.Nsamples:Direct.Nsamples)/Data.Fs;

Direct.H = fft(Direct.InnSph.h,length(TDOA.p))/length(TDOA.p);

TDOA.H = Direct.H(:,[59 69]);

TDOA.R = mean(TDOA.H(:,1).*conj(TDOA.H(:,2)),2);

% Frequency-Domain Weighting Function
if strcmp(TDOA.Method,'CC')
    TDOA.theta = 1;
elseif strcmp(TDOA.Method,'PHAT')
    TDOA.theta = 1./abs(TDOA.R);
end

% Generalised Cross Spectrum
TDOA.GCS = TDOA.theta.*TDOA.R;

% Time Domain
TDOA.GCC = fftshift(ifft(TDOA.GCS,[],'symmetric'));

plot(TDOA.p*1e3,TDOA.GCC), grid on
disp(TDOA.p(TDOA.GCC==max(TDOA.GCC))*1e3)

%% DOA Estimation via TDOA
TDOA.Method = 'PHAT';
TDOA.p = (-Direct.Nsamples/2:Direct.Nsamples/2-1)/Data.Fs;

% TDOA.R = nan(Data.InnSph.M,Data.InnSph.M,Direct.Nsamples);
% for ii = 1:Data.InnSph.M
%     for jj = 1:Data.InnSph.M
%         TDOA.R(ii,jj,:) = mean(Direct.H(:,ii).*conj(Direct.H(:,jj)),2);
%     end
% end
TDOA.R = nan(2,2,Direct.Nsamples);
for ii = 1:2
    for jj = 1:2
        TDOA.R(ii,jj,:) = mean(TDOA.H(:,ii).*conj(TDOA.H(:,jj)),2);
    end
end

% Frequency-Domain Weighting Function
if strcmp(TDOA.Method,'CC')
    TDOA.theta = 1;
elseif strcmp(TDOA.Method,'PHAT')
    TDOA.theta = 1./abs(TDOA.R);
end

% Generalised Cross Spectrum
TDOA.GCS = TDOA.theta.*TDOA.R;

% Time Domain
TDOA.GCS_Full = cat(3,TDOA.GCS(:,:,1),TDOA.GCS(:,:,2:end)/2);
TDOA.GCS_Full = cat(3,TDOA.GCS_Full,flip(conj(TDOA.GCS_Full)));

TDOA.GCC = fftshift(ifft(TDOA.GCS_Full,Direct.Nsamples,3,'symmetric'),3);
% TDOA.GCC = fftshift(ifft(TDOA.GCS_Full,[],3,'symmetric'),3);

TDOA.Ra = nan(Direct.Nsamples,1);
for ii = 1:Direct.Nsamples
    TDOA.Ra(ii) = det(squeeze(TDOA.GCC(:,:,ii)));
end

figure, plot(TDOA.p*1e3,TDOA.Ra), grid on


