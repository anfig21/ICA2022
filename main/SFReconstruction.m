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
Data.Temp = 23.2;
Data.Humi = 34.7;
Data.p0 = 100658;

[Data.rho,Data.c,~,Data.gamma,~,~,~,~,~,~] = amb2prop(Data.p0,Data.Temp,...
    Data.Humi,1);

Data.loudspeaker = 'rir_019_spk1.h5';
% Data.loudspeaker = 'rir_019_spk2.h5';

%% DATA ACQUISITION
Data = dataAcquisition(Data);

%% DATA HANDLING
Data.Ref.pos = vertcat(Data.Line1.pos,Data.Line2.pos,Data.Line3.pos);
Data.Ref.h = horzcat(Data.Line1.h,Data.Line2.h,Data.Line3.h);
Data.Ref.n = horzcat(Data.Line1.n,Data.Line2.n,Data.Line3.n);

% Inner Sphere (156 Microphones: samples 155-end)
Idx = 155:size(Data.Sph.pos,1);
Data.InnSph.M = length(Idx);

% Random sampling (N = 150)
% Data.InnSph.M = 150;
% Idx = randperm(size(Data.Sph.pos,1),Data.InnSph.M);

Data.InnSph.pos = Data.Sph.pos(Idx,:);
Data.InnSph.h = Data.Sph.h(:,Idx);
Data.InnSph.p = Data.Sph.p(:,Idx);

% Frequency domain
Data.InnSph.H = fft(Data.InnSph.h,2*Data.Nsamples)/(2*Data.Nsamples);
Data.InnSph.H = [Data.InnSph.H(1,:); 2*Data.InnSph.H(2:end/2,:)];
Data.InnSph.P = fft(Data.InnSph.p,2*Data.Nsamples)/(2*Data.Nsamples);
Data.InnSph.P = [Data.InnSph.P(1,:); 2*Data.InnSph.P(2:end/2,:)];
Data.Sph.N = fft(Data.Sph.n,2*Data.Nsamples)/(2*Data.Nsamples);
Data.Sph.N = [Data.Sph.N(1,:); 2*Data.Sph.N(2:end/2,:)];

%% SETUP PLOT
Plot.T = [5 25]*1e-3;       % Source near field
% Plot.T = [15 35]*1e-3;    % Source far field
Plot.N = Data.Fs*Plot.T;

% Time vector
Plot.t = Plot.T(1):1/Data.Fs:Plot.T(2)-(1/Data.Fs);

% Data downsizing
Plot.Ref.h = Data.Ref.h(Plot.N(1):Plot.N(2)-1,:);
Plot.InnSph.h = Data.InnSph.h(Plot.N(1):Plot.N(2)-1,:);

% plotFreqResponse(Data)

% REFERENCE LINE RIR PLOT
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
clear s c

%% Dictionary of plane waves
% Dictionary
% Dict.f = Data.f(6e2 <= Data.f & Data.f <= 1e3);
% Dict.f = Dict.f(1:50:end);
% Dict.f = Dict.f(1:20:end);
Dict.f = 8e2;                   % DOA estimation at 800 Hz
Dict.Plane.N = 1e3;             % Number of plane waves
Dict.Plane.K = 1;               % Number of sources (SOMP)

[Dict.Plane.H,Dict.Plane.uk] = dictionary(Data.c,Dict.f,Data.InnSph.pos',Dict.Plane.N);

%% ------------ DIRECT SOUND FIELD ------------ %%
% Time window: 5-10 ms
Direct.T = 8*1e-3;     % Source near field
% Direct.T = 22*1e-3;   % Source far field

Direct = windowRIR(Data,0,Direct.T);

% True DOA
Direct.TrueDOA = Data.Source.pos-Data.Sph.R0;
Direct.TrueDOA = Direct.TrueDOA/vecnorm(Direct.TrueDOA);

%% DOA Estimation via SOMP
% Direct.DOA.SOMP = dirDOA_SOMP(Data,Direct,Dict,'true');

%% DOA Estimation via Least-Squares
% Direct.DOA.LS = dirDOA_LS(Data,Direct,Dict,'true');

%% DOA Estimation via Regularised Least-Squares
% Direct.DOA.RLS = dirDOA_RLS(Data,Direct,Dict,'true');
% Nnorm = 1.685*10-6 via L-curve method
% Direct.InnSph.NnormLcurve = 1.0605e-7;       % Varies with frequency

%% DOA Estimation via Compressive Sensing
Direct.DOA.CS = dirDOA_CS(Data,Direct,Dict,'true');

%% Dictionary of spherical waves
Dict.Sph.Res = 0.1;
Dict.Sph.rMinMax = [0.3 3];

[Dict.Sph.H,Dict.Sph.r,Dict.Sph.N] = dictionaryRange(Data.c,Dict.f,Data.InnSph.pos.',...
    Data.Sph.R0.',Direct.DOA.CS.Avg.',Dict.Sph.rMinMax,Dict.Sph.Res);

% [Dict.Sph.H,Dict.Sph.r,Dict.Sph.N] = dictionaryRange(Data.c,Dict.f,Data.InnSph.pos.',...
%     Data.Sph.R0.',Direct.TrueDOA,Dict.Sph.rMinMax,Dict.Sph.Res);

disp('Direct sound: Spherical Wave Dictionary... OK')

%% Range Estimation via Regularised Least-Squares
% Direct.Range.RLS = dirRange_RLS(Data,Direct,Dict,'true');
% Direct.InnSph.NnormLcurve = 6.8e-4;       % Varies with frequency

%% Range Estimation via CS
% Direct.Range.CS = dirRange_CS(Data,Direct,Dict,'true');

%% Range Estimation via EN
% Direct.Range.EN = dirRange_EN(Data,Direct,Dict,'true');

%% Range Estimation via TV
Direct.Range.TV = dirRange_TV(Data,Direct,Dict,'true');
% figure, plot(Dict.f,Direct.Range.TV.Error), grid on
% OPTIMAL IDX: 19

%% ------------ EARLY REFLECTIONS ------------ %%
TotT = [8 10; 11 15; 16 18]*1e-3;
Early.totR = 3;

for ii = 1:Early.totR
    % Time window: 10-20 ms
    %     Early.T = [8 18]*1e-3;      % Source near field
    Early.T = TotT(ii,:);
    
    Early = windowRIR(Data,Early.T(1),Early.T(2));
    Early.R = 1;                % Number of early reflections
    
    %% DOA Estimation via Compressive Sensing
    Early.DOA.CS{ii} = earlyDOA_CS(Data,Early,Dict,'true');     % OUTPERFORMS LASSO
    
    %% DOA Estimation via Elastic Net
    % Early.DOA.EN = earlyDOA_EN(Data,Early,Dict,'true');   % DOES NOT WORK
    
    %% DOA Estimation via Weighted LASSO
    % Early.DOA.WL = earlyDOA_WL(Data,Early,Dict,'true');
    
    %% Dictionary of spherical waves
    Dict.SphEarly.Res = 0.1;
    Dict.SphEarly.rMinMax = [1 5];
    
    % Dictionary per DOA
    
    [Dict.SphEarly.H{1},Dict.SphEarly.r{1},Dict.SphEarly.N{1}] = dictionaryRange(Data.c,Dict.f,Data.InnSph.pos.',...
        Data.Sph.R0.',Early.DOA.CS{ii}.Est.',Dict.SphEarly.rMinMax,Dict.SphEarly.Res);
    
    
    %% Range Estimation via Regularised Least-Squares
    % Early.Range.RLS = earlyRange_RLS(Data,Early,Dict,'true');
    
    %% Range Estimation via CS
    % Early.Range.CS = earlyRange_CS(Data,Early,Dict,'true');
    
    %% Range Estimation via TV
    Early.Range.TV = earlyRange_TV(Data,Early,Dict,'true');
    
end
clear ii
disp('Early reflections: Spherical Wave Dictionary... OK')

%% ------------ PLANE WAVE EXPANSION ------------ %%
PWE.T = 30*1e-3;     % Source near field

PWE = windowRIR(Data,0,PWE.T);

Dict.PWE.N = 1e3;    % Number of plane waves
PWE.f = Data.f(0 <= Data.f & Data.f <= 1e3);
PWE.f = PWE.f(1:50:end);

[Dict.PWE.H,~] = dictionary(Data.c,PWE.f,Data.InnSph.pos',Dict.PWE.N);

%% Coefficient estimation
PWE.x = nan(Dict.PWE.N,length(PWE.f));
for ii = 1:length(PWE.f)
    [PWE.x(:,ii),~] = reguLeastSquares(squeeze(Dict.PWE.H(:,:,ii)),Direct.InnSph.H(Data.f==PWE.f(ii),:).');
end
Dict.PWE = rmfield(Dict.PWE,'H');

%% Reconstruction
[Dict.PWE.HR,~] = dictionary(Data.c,PWE.f,Data.Ref.pos',Dict.PWE.N);

PWE.P = nan(size(Data.Ref.pos,1),length(PWE.f));
for ii = 1:length(PWE.f)
    PWE.P(:,ii) = squeeze(Dict.PWE.HR(:,:,ii))*PWE.x(:,ii);
end
Dict.PWE = rmfield(Dict.PWE,'HR');

%% Inverse Fourier Transform
% Double-sided spectrum
PWE.P2 = [real(PWE.P(:,1)) PWE.P(:,2:end-1)/2 PWE.P(:,end)];
PWE.P2 = [PWE.P2 flip(conj(PWE.P2(:,2:end-1)),2)]*(2*Data.Nsamples);
PWE.p = ifft(PWE.P2,[],2,'symmetric');






