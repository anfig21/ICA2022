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
% anfig@dtu.dk
%
% January 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all

addpath(genpath('../tools'))
addpath(genpath('/Volumes/anfig/Data/room019/'))


%% INITIAL PARAMETERS
T = 5;
Fs = 48e3;
D = [6.266,9.357,2.977];    % Room dimensions

time = [5e-3,25e-3];   % Plot lenght
Nsamples = Fs*T;
Nplot = Fs*time;

loudspeaker = 'rir_019_spk1.h5';

%% DATA ACQUISITION
% Line 1
pos1 = h5read(loudspeaker,'/dataset_lin2/position');
h1 = h5read(loudspeaker,'/dataset_lin2/impulse_response');

% Line 2
pos2 = h5read(loudspeaker,'/dataset_lin1/position');
h2 = h5read(loudspeaker,'/dataset_lin1/impulse_response');

% Line 3
pos3 = h5read(loudspeaker,'/dataset_lin3/position');
h3 = h5read(loudspeaker,'/dataset_lin3/impulse_response');

%% PLOT
% Time vector
t = time(1):1/Fs:time(2)-(1/Fs);

% Data downsizing
h1_plot = h1(Nplot(1):Nplot(2)-1,:);
h2_plot = h2(Nplot(1):Nplot(2)-1,:);
h3_plot = h3(Nplot(1):Nplot(2)-1,:);

% Data merging
h_plot = horzcat(h1_plot,h2_plot,h3_plot);
pos = vertcat(pos1,pos2,pos3);

% Line plot
%scatter3(pos(:,1),pos(:,2),pos(:,3))
%axis([0 D(1) 0 D(2) 0 D(3)])

%% Impulse response plot
figure
s = surf(pos(:,1),t,h_plot);
set(s,'edgecolor','none')
colormap hot
view(2)
colorbar
caxis([-0.04 0.04])
