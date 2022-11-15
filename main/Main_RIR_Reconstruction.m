%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  ROOM IMPULSE RESPONSE RECONSTRUCTION
%
% -------------------------------------------------------------------------
% Figueroa-Duran, Fernandez-Grande, "Reconstruction of room impulse
% responses over an extended spatial domain using block-sparse and kernel
% regression methods", International Congress of Acoustics, 2022
% -------------------------------------------------------------------------
%
% Antonio Figueroa-Duran
% anfig@dtu.dk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all

addpath(genpath(pwd))
addpath(genpath('../../tools'))
% addpath(genpath('/Volumes/anfig/Data/room019/'))      % MacBook Pro
addpath(genpath('M:\Data'))                             % Windows 10

loadPlotParams

%% INITIAL PARAMETERS
Data.T = 1;                         % Pre-processing T = 1 s
Data.Fs = 48e3;
Data.t = 0:1/Data.Fs:Data.T-1/Data.Fs;
Data.D = [6.266,9.357,2.977];       % Room dimensions
Data.Nsamples = Data.Fs*Data.T;
Data.f2 = Data.Fs/Data.Nsamples*(-Data.Nsamples/2:Data.Nsamples/2-1);
Data.f = Data.Fs/Data.Nsamples*(0:Data.Nsamples/2);
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
Data = dataHandling(Data);

%% SETUP PLOT
% Flags
% - Setup
% - Frequency response
% - RIR reference line
Plot = setupPlot(Data,false,false,true);

%% ------------ DIRECT SOUND FIELD ------------ %%
% Windowing
Direct.T = 8*1e-3;      % Source near field
% Direct.T = 22*1e-3;     % Source far field

Direct = windowRIR(Data,0,Direct.T);

%% DOA Estimation
N = 1e3;                % Number of plane waves
DOAMethod = 'RLS';      % DOA Estimation Method
Direct.f = Data.f(1.5e3 <= Data.f & Data.f <= 2e3);
Direct.f = Direct.f(1:10:end);

Direct = directSoundDOA(Data,Direct,N,Direct.f,DOAMethod);
clear N DOAMethod

%% Range Estimation
res = 0.1;
rMinMax = [0.3 5];
RangeMethod = 'TV';    % Range Estimation Method

Direct = directSoundRange(Data,Direct,res,rMinMax,Direct.f,RangeMethod,true);
clear res rMinMax RangeMethod

%% Reconstruction
Rec.T = [5 7.5]*1e-3;
Rec.Direct = windowRIR(Data,Rec.T(1),Rec.T(2));

res = .5e-1;     % Spatial resolution (m)
L = 0.05;          % Cube side (m)
Rec.Direct.f = Data.f(1 <= Data.f & Data.f <= 5e3);
RecMethod = 'CS';

% Rec.Direct = reconstructReflection(Data,Rec.Direct,Direct.Range.Mode,L,res,RecMethod);
Rec.Direct = reconstructReflection(Data,Rec.Direct,Data.Source.pos,L,res,RecMethod);
clear res L RecMethod

% Plot
figure, hold on
s = pcolor(Data.Ref.pos(:,1),Plot.t*1e3,Rec.Direct.h(Plot.N(1):Plot.N(2),:));
set(s,'edgecolor','none')
xlabel('x in m'), ylabel('Time in ms')
colormap hot
c = colorbar;
caxis([-0.04 0.04])
xline(min(Data.InnSph.pos(:,1))), xline(max(Data.InnSph.pos(:,1)))
applyColorbarProperties(c,'Room Impulse Response in Pa/V')
applyAxisProperties(gca)

% Coefficients
figure, plot(abs(Rec.Direct.x).')

%% ------------ EARLY REFLECTIONS ------------ %%
% Windowing
% TotT = [8 10; 11 15; 16 18]*1e-3;
Early.T = [8 10]*1e-3;

Early = windowRIR(Data,Early.T(1),Early.T(2));

%% DOA Estimation
% Window & look for reflections individualy
N = 1e3;                % Number of plane waves
DOAMethod = 'RLS';      % DOA Estimation Method
Early.f = Data.f(1e3 <= Data.f & Data.f <= 1.5e3);
Early.f = Early.f(1:10:end);
Early.R = 1;

Early = earlyReflectionsDOA(Data,Early,N,Early.f,DOAMethod,true);
clear N DOAMethod

%% Range Estimation
res = 0.1;
rMinMax = [0.3 5];
RangeMethod = 'TV';    % Range Estimation Method

Early = earlySoundRange(Data,Early,res,rMinMax,Early.f,RangeMethod,true);
clear res rMinMax RangeMethod

%% Reconstruction
Rec.h = Rec.Direct.h;

PD.T0 = 6e-3;
PD.T = [9.58333e-3 11.4375e-3 17.6667e-3]';

PD.Range = vecnorm(Data.Sph.R0(:)-Data.Ref.pos(end/2,:)');
PD.Range = PD.Range+(PD.T-PD.T0)*Data.c;                        % Distance estimated from time delay

PD.TotT = PD.T+[-1.2 1.2]*1e-3;

Early.T = PD.TotT(1,:);

Rec.T = Early.T;
Rec.Early = windowRIR(Data,Rec.T(1),Rec.T(2));

res = 1e-1;     % Spatial resolution (m)
L = 0.1;          % Cube side (m)
Rec.Early.f = Data.f(1 <= Data.f & Data.f <= 3e3);
RecMethod = 'CS';

Rec.Early = reconstructReflection(Data,Rec.Early,Early.Range.Mode,L,res,RecMethod);
clear res L RecMethod

Rec.h = Rec.h+Rec.Early.h;

% Plot
figure, hold on
s = pcolor(Data.Ref.pos(:,1),Plot.t*1e3,Rec.Early.h(Plot.N(1):Plot.N(2),:));
set(s,'edgecolor','none')
xlabel('x in m'), ylabel('Time in ms')
colormap hot
c = colorbar;
caxis([-0.04 0.04])
xline(min(Data.InnSph.pos(:,1))), xline(max(Data.InnSph.pos(:,1)))
applyColorbarProperties(c,'Room Impulse Response in Pa/V')
applyAxisProperties(gca)

%% ------------ PLANE WAVE EXPANSION ------------ %%
% PWE in time domain does not work as expected

% Windowing
PWE.Freq.T = 1;
PWE.Freq = windowRIR(Data,0,PWE.Freq.T);

%% PWE in frequency
PWE.Freq.N = 1e3;    % Number of plane waves
PWE.Freq.f = Data.f(0 <= Data.f & Data.f <= 2e3);

PWE.Freq = planeWaveExpansionFrequency(Data,PWE.Freq,Plot,PWE.Freq.N,PWE.Freq.f,Data.Ref.pos.','RLS',true);

%% PWE in time
PWE.Time.N = 1e3;
PWE.Time.t = Data.t(Data.t >= 4e-3 & Data.t <= 25e-3);
PWE.Time.Q = 2;     % Downsampling factor

tic
PWE.Time = planeWaveExpansionTime(Data,PWE.Time.N,PWE.Time.t,PWE.Time.Q,Data.Ref.pos.',true);
toc

figure, hold on
plot(PWE.Time.t*1e3,PWE.Time.h(:,220)), grid on
plot(Data.t*1e3,Data.Ref.h(:,220))
xlabel('Time in ms'), ylabel('Impulse response')
legend('Rec','True')
xlim([PWE.Time.t(1) PWE.Time.t(end)]*1e3)
applyAxisProperties(gca)
applyLegendProperties(gcf)

%% ------------ BESSEL KERNEL ------------ %%
KRR.T = [69 70]*1e-3;
% KRR.T = [0 1];

KRR = windowRIR(Data,KRR.T(1),KRR.T(end));
KRR.t = Data.t(Data.t >= KRR.T(1) & Data.t <= KRR.T(2));
KRR.f = Data.f(0 <= Data.f & Data.f <= 5e3);

NoiseMargin = -15;           % dB
M = Data.InnSph.M;
R = Data.Ref.pos.';
Nr = size(R,2);
Nt = length(KRR.t);
KRR.H = zeros(Nr,length(Data.f));

% Gram matrices
d_MM = nan(M);
d_MNr = nan(M,Nr);
for ii = 1:M
    d_MM(ii,:) = vecnorm(repmat(Data.InnSph.pos(ii,:),M,1)-Data.InnSph.pos,2,2).';
    d_MNr(ii,:) = vecnorm(repmat(Data.InnSph.pos(ii,:),Nr,1)-R.',2,2).';
end

% SF Reconstruction
Nf = length(KRR.f);
k = (2*pi*KRR.f)/Data.c;

% Standard deviation
sigma_2 = std(KRR.InnSph.H(ismember(Data.f,KRR.f),:),[],2);

for ii = 1:Nf
    p = Data.InnSph.H(Data.f==KRR.f(ii),:).';

    % Multiply by k
    K_MM = sigma_2(ii)*sinc(k(ii)*d_MM/pi);
    K_MNr = sigma_2(ii)*sinc(k(ii)*d_MNr/pi).';
    
    % Regularisation parameter
    [U,s,V] = csvd(K_MM);
%     [lambda,~,~,~] = l_curve(U,s,p,'Tikh',[],[],false);
%     [lambda,~,~] = gcv(U,s,p,'Tikh');
    [lambda,~,~] = quasiopt(U,s,p,'Tikh');

%     lambda = 10^(NoiseMargin/20)*KRR.InnSph.Nnorm(Data.f==KRR.f(ii));
%     lambda = 2*KRR.InnSph.Nnorm(Data.f==KRR.f(ii));
    
    KRR.H(:,Data.f==KRR.f(ii)) = K_MNr*((K_MM+lambda*eye(M))\p);
end

mic = 220;
figure, hold on
plot(Data.f,20*log10(abs(Data.Ref.H(:,mic))))
plot(KRR.f,20*log10(abs(KRR.H(mic,ismember(Data.f,KRR.f))))), grid on
xlabel('Frequency in Hz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
legend('True','Rec')
xlim([KRR.f(1) KRR.f(end)])
applyAxisProperties(gca)
applyLegendProperties(gcf)

% Double-sided spectrum
KRR.P2 = [real(KRR.H(:,1)) KRR.H(:,2:end)/2];
KRR.P2 = [KRR.P2 flip(conj(KRR.P2(:,2:end)),2)];
KRR.h = ifft(KRR.P2*Data.Nsamples,[],2,'symmetric').';

figure, hold on
s = pcolor(Data.Ref.pos(:,1),KRR.t*1e3,KRR.h(ismember(Data.t,KRR.t),:));
set(s,'edgecolor','none')
xlabel('x in m'), ylabel('Time in ms')
colormap hot
c = colorbar;
caxis([-0.04 0.04])
xline(min(Data.InnSph.pos(:,1))), xline(max(Data.InnSph.pos(:,1)))
applyColorbarProperties(c,'Room Impulse Response in Pa/V')
applyAxisProperties(gca)

%% Aperture plot
t0 = 69.5e-3;
figure, hold on
plot(Data.Ref.pos(:,1),Data.Ref.h(find(Data.t > t0,1),:))
plot(Data.Ref.pos(:,1),KRR.h(find(Data.t > t0,1),:)), grid on
xline(min(Data.InnSph.pos(:,1))), xline(max(Data.InnSph.pos(:,1)))
xlabel('x in m'), ylabel('Room Impulse Response in Pa/V')
legend('True','Rec')
applyAxisProperties(gca)
applyLegendProperties(gcf)

%% Time plot
figure, hold on
plot(KRR.t*1e3,Data.Ref.h(ismember(Data.t,KRR.t),mic))
plot(KRR.t*1e3,KRR.h(ismember(Data.t,KRR.t),mic)), grid on
xlabel('Time in ms'), ylabel('Room Impulse Response in Pa/V')
legend('True','Rec')
applyAxisProperties(gca)
applyLegendProperties(gcf)




