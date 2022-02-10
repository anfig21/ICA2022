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

% Sphere
posSph = h5read(loudspeaker,'/dataset_sph1/position');
hSph = h5read(loudspeaker,'/dataset_sph1/impulse_response');

%% DATA MERGING
h_ref = horzcat(h1,h2,h3);
posLin = vertcat(pos1,pos2,pos3);

% Inner Sphere (samples 155-end)
posInnSph = posSph(155:end,:);
hInnSph = hSph(:,155:end);

%% SETUP PLOT
% Time vector
t = time(1):1/Fs:time(2)-(1/Fs);

% Data downsizing
href2 = h_ref(Nplot(1):Nplot(2)-1,:);
hInnSph2 = hInnSph(Nplot(1):Nplot(2)-1,:);

% Line plot
figure
scatter3(posLin(:,1),posLin(:,2),posLin(:,3)), hold on
scatter3(posInnSph(:,1),posInnSph(:,2),posInnSph(:,3))
axis([0 D(1) 0 D(2) 0 D(3)])

%% REFERENCE LINE RIR PLOT
figure
s = surf(posLin(:,1),t,href2);
set(s,'edgecolor','none')
colormap hot
view(2)
colorbar
caxis([-0.04 0.04])

%% DIRECT SOUND FIELD
% Time window: 0-50 ms




