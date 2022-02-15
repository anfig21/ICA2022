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

Plot.T = [5 25]*1e-3;            % Plot lenght
Plot.Nsamples = Data.Fs*Plot.T;
Data.loudspeaker = 'rir_019_spk1.h5';

%% DATA ACQUISITION
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

%% DATA MERGING
Data.Ref.h = horzcat(Data.Line1.h,Data.Line2.h,Data.Line3.h);
Data.Ref.pos = vertcat(Data.Line1.pos,Data.Line2.pos,Data.Line3.pos);

% Inner Sphere (156 Microphones: samples 155-end)
Data.InnSph.h = Data.Sph.h(:,155:end);
Data.InnSph.pos = Data.Sph.pos(155:end,:);

%% SETUP PLOT
% Time vector
Plot.t = Plot.T(1):1/Data.Fs:Plot.T(2)-(1/Data.Fs);

% Data downsizing
Plot.Ref.h = Data.Ref.h(Plot.Nsamples(1):Plot.Nsamples(2)-1,:);
Plot.InnSph.h = Data.InnSph.h(Plot.Nsamples(1):Plot.Nsamples(2)-1,:);

% Line plot
figure
scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
legend('Reference Line','Spherical Array')
applyAxisProperties(gca)
applyLegendProperties(gca)

%% REFERENCE LINE RIR PLOT
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

%% DIRECT SOUND FIELD
% Time window: 5-10 ms
Direct.T = [5 10]*1e-3;
Direct.Nsamples = Data.Fs*Direct.T;

% Data Downsizing
Direct.InnSph.h = Data.InnSph.h(Direct.Nsamples(1):Direct.Nsamples(2)-1,:);

% DOA Estimation
DOA.N = 500;        % Number of plane waves
[DOA.H,DOA.uk] = dictionary(Data.f,Data.InnSph.pos',DOA.N);


