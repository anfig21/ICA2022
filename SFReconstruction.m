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

addpath(genpath('data'))
%% INITIAL PARAMETERS
folderName = ["robot_pos01","robot_pos02","robot_pos03"];
N = 152;
T = 10;         % Measurement length
time = 25e-3;   % Plot lenght
Fs = 48e3;      % Sampling frequency
Nsamples = Fs*T;
Nplot = Fs*time;

%% DATA ACQUISITION
% Reference grid
fileName = '/data_line_pos';

positions = nan(N,3);
data = nan(N,Nsamples);
for ff = 1:N
    load(strcat(folderName(1),fileName,string(ff)),'pos','dataMic')
    positions(ff,:) = pos;
    data(ff,:) = dataMic;
end
clear ff pos dataMic

%% PLOT
t = 0:1/Fs:time-1/Fs;
dataPlot = data(:,1:Nplot);

figure
s = surf(t,positions(:,1),dataPlot);
s.EdgeColor = 'none';
view(2)

