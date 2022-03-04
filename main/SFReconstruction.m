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

% Data.loudspeaker = 'rir_019_spk1.h5';
Data.loudspeaker = 'rir_019_spk2.h5';

%% DATA ACQUISITION
% Source
Data.Source.pos = h5read(Data.loudspeaker,'/source/position');

% Line 1
Data.Line1.pos = h5read(Data.loudspeaker,'/dataset_lin2/position');
Data.Line1.h = h5read(Data.loudspeaker,'/dataset_lin2/impulse_response');
Data.Line1.n = h5read(Data.loudspeaker,'/dataset_lin2/noise');

% Line 2
Data.Line2.pos = h5read(Data.loudspeaker,'/dataset_lin1/position');
Data.Line2.h = h5read(Data.loudspeaker,'/dataset_lin1/impulse_response');
Data.Line2.n = h5read(Data.loudspeaker,'/dataset_lin1/noise');

% Line 3
Data.Line3.pos = h5read(Data.loudspeaker,'/dataset_lin3/position');
Data.Line3.h = h5read(Data.loudspeaker,'/dataset_lin3/impulse_response');
Data.Line3.n = h5read(Data.loudspeaker,'/dataset_lin3/noise');

% Sphere
Data.Sph.pos = h5read(Data.loudspeaker,'/dataset_sph1/position');
Data.Sph.h = h5read(Data.loudspeaker,'/dataset_sph1/impulse_response');
Data.Sph.n = h5read(Data.loudspeaker,'/dataset_sph1/noise');
Data.Sph.R0 = mean(Data.Sph.pos,1);

%% DATA HANDLING
Data.Ref.h = horzcat(Data.Line1.h,Data.Line2.h,Data.Line3.h);
Data.Ref.n = horzcat(Data.Line1.n,Data.Line2.n,Data.Line3.n);
Data.Ref.pos = vertcat(Data.Line1.pos,Data.Line2.pos,Data.Line3.pos);

% Inner Sphere (156 Microphones: samples 155-end)
Data.InnSph.h = Data.Sph.h(:,155:end);
Data.InnSph.pos = Data.Sph.pos(155:end,:);
Data.InnSph.M = size(Data.InnSph.pos,1);

% Frequency domain
Data.InnSph.H = fft(Data.InnSph.h,2*Data.Nsamples)/Data.Nsamples;
Data.InnSph.H = [Data.InnSph.H(1,:); 2*Data.InnSph.H(2:end/2,:)];
Data.Sph.N = fft(Data.Sph.n,2*Data.Nsamples)/Data.Nsamples;
Data.Sph.N = [Data.Sph.N(1,:); 2*Data.Sph.N(2:end/2,:)];

% Noise norm
Data.Sph.Nnorm = vecnorm(Data.Sph.N);

%% SETUP PLOT
Plot.T = [5 25]*1e-3;               % Plot lenght
Plot.N = Data.Fs*Plot.T;

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
% applyLegendProperties(gcf)

% Frequency response
% figure, plot(Data.f*1e-3,20*log10(abs(Data.InnSph.H(:,[40 59 69])))), grid on
% xlabel('Frequency in kHz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
% legend('Mic 40','Mic 59','Mic 69')
% applyAxisProperties(gca)
% applyLegendProperties(gcf)

% figure, plot(Data.f*1e-3,20*log10(abs(Data.Sph.N)/20e-6)), grid on
% xlabel('Frequency in kHz'), ylabel('Magnitude $|N(j\omega)|$ in dB ref. 20\,$\mu Pa$')
% applyAxisProperties(gca)
% ylim([-5 40])

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
Direct.InnSph.Nnorm = vecnorm(Direct.InnSph.N);

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
Dict.f = Data.f(180 <= Data.f & Data.f <= 240);
Dict.f = Dict.f(1:20:end);
% Dict.f = 2e2;          % DOA estimation at 200 Hz
Dict.N = 1e3;           % Number of plane waves
Dict.K = 1;             % Number of sources (SOMP)
Dict.Est = nan(3,length(Dict.f));

[Dict.CA,Dict.uk] = dictionary(Dict.f,Data.InnSph.pos',Dict.N);

%% DOA Estimation via SOMP
%DOA = directSOMP(Data,Direct,Dict,'true');

%% DOA Estimation via Least-Squares
%LS = directLS(Data,Direct,Dict,'true');

%% DOA Estimation via Regularised Least-Squares
%RLS = directRLS(Data,Direct,Dict,'true');

%% DOA Estimation via Compressive Sensing
CS.Nnorm = mean(Direct.InnSph.Nnorm);

CS.Est = nan(3,length(Dict.f));
tic
for ii = 1:length(Dict.f)
    Hii = squeeze(Dict.CA(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dict.f(ii),:).';
    Nnorm = 1.1*CS.Nnorm;
    
    % CVX Formulation
    cvx_begin
    cvx_quiet(true)
    cvx_precision high
    variable x(Dict.N) complex;
    minimize norm(x,1);
    subject to
    norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    [~,Idx] = max(abs(x));
    CS.Est(:,ii) = Dict.uk(:,Idx);
end
toc
CS.Error = vecnorm(Direct.DOA.'-CS.Est);

CS.Avg = -mean(CS.Est,2);
CS.Avg = CS.Avg/vecnorm(CS.Avg);

if true
    % Mean Squared Error
    figure, plot(Dict.f,CS.Error), grid on
    xlabel('Frequency in Hz'), ylabel('DOA - Mean Squared Error')
    applyAxisProperties(gca)
    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
    scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3))
    quiver3(Data.Sph.R0(1),Data.Sph.R0(2),Data.Sph.R0(3),CS.Avg(1),CS.Avg(2),CS.Avg(3),2,'Linewidth',4)
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Spherical Array','Source')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

%% DOA Estimation via TDOA with only 2 mics (Mics 59 & 69)
TDOA.Method = 'PHAT';
TDOA.p = (-Direct.Nsamples:Direct.Nsamples)/Data.Fs;

Direct.InnSph.H = fft(Direct.InnSph.h,length(TDOA.p))/length(TDOA.p);

TDOA.H = Direct.InnSph.H(:,[59 69]);

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
%         TDOA.R(ii,jj,:) = mean(Direct.InnSph.H(:,ii).*conj(Direct.InnSph.H(:,jj)),2);
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


