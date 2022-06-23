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
Data.T = 1;                         % Pre-processing T = 1 s
Data.Fs = 48e3;
Data.D = [6.266,9.357,2.977];       % Room dimensions
Data.Nsamples = Data.Fs*Data.T;
Data.f2 = Data.Fs/Data.Nsamples*(-Data.Nsamples/2:Data.Nsamples/2-1);
Data.f = Data.f2(end/2+1:end);
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
Plot = setupPlot(Data,false,false,false);

%% ------------ DIRECT SOUND FIELD ------------ %%
% Parameters:
% Resolution SWs
% Range SWs
% Range Method

%% Windowing
Direct.T = 8*1e-3;      % Source near field
% Direct.T = 22*1e-3;     % Source far field

Direct = windowRIR(Data,0,Direct.T);

%% DOA Estimation
N = 1e3;                % Number of plane waves
DOAMethod = 'RLS';      % DOA Estimation Method
Direct.f = Data.f(6e2 <= Data.f & Data.f <= 2e3);

Direct = directSoundDOA(Data,Direct,N,Direct.f,DOAMethod,true);
clear N DOAMethod

%% Range Estimation
res = 0.1;
rMinMax = [0.3 6];
RangeMethod = 'TV';    % Range Estimation Method

Direct = directSoundRange(Data,Direct,res,rMinMax,Direct.f,RangeMethod,true);
clear res rMinMax RangeMethod

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
    Early.DOA.CS{ii} = earlyDOA_CS(Data,Early,Dict);     % OUTPERFORMS LASSO
    
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
    Early.Range.TV = earlyRange_TV(Data,Early,Dict);
    
end
clear ii
disp('Early reflections: Spherical Wave Dictionary... OK')

%% ------------ PLANE WAVE EXPANSION ------------ %%
PWE.T = 1;
PWE = windowRIR(Data,0,PWE.T);

Dict.PWE.N = 2.5e3;    % Number of plane waves
PWE.f = Data.f(1 <= Data.f & Data.f <= 1.5e3);
PWE.f = PWE.f(1:10:end);

% PWE.f = 500;
%% Coefficient estimation
PWE.x = nan(Dict.PWE.N,length(PWE.f));
for ii = 1:length(PWE.f)    
    [H,~] = dictionary(Data.c,PWE.f(ii),Data.InnSph.pos',Dict.PWE.N);
%     [PWE.x(:,ii),~] = reguLeastSquares(H,PWE.InnSph.H(Data.f==PWE.f(ii),:).');

    Nnorm = 10^(10/20)*PWE.InnSph.Nnorm(Data.f==PWE.f(ii));
    p = PWE.InnSph.H(Data.f==PWE.f(ii),:).';
    % COMPRESSIVE SENSING
    cvx_begin quiet
    cvx_precision high
        variable x(Dict.PWE.N) complex;
        minimize norm(x,1);
        subject to
            norm((H*x-p),2) <= Nnorm;
    cvx_end

    PWE.x(:,ii) = x;

    disp(strcat("Freq... ",string(round(ii/length(PWE.f)*100,2)),"%"))
end
clear H
% Dict.PWE = rmfield(Dict.PWE,'H');

%% Reconstruction
% Mic = 200;
% Mic = 230;
Mic = 10;

[Dict.PWE.HR,~] = dictionary(Data.c,PWE.f,Data.InnSph.pos(Mic,:)',Dict.PWE.N);

% PWE.P = zeros(1,length(Data.f));
PWE.P = zeros(1,length(PWE.f));
for ii = 1:length(PWE.f)
    PWE.P(:,ii) = squeeze(Dict.PWE.HR(:,:,ii))*PWE.x(:,ii);
end
Dict.PWE = rmfield(Dict.PWE,'HR');

%% Inverse Fourier Transform
% Double-sided spectrum
PWE.P2 = [real(PWE.P(:,1)) PWE.P(:,1:end-1)/2 PWE.P(:,end)];
PWE.P2 = [PWE.P2 flip(conj(PWE.P2(:,2:end-1)),2)];
PWE.h = ifft(PWE.P2*Data.Nsamples,[],2,'symmetric').';

%% Frequency Domain Plot
HTrue = 20*log10(abs(Data.InnSph.H(:,Mic)));
% HTrue = HTrue-max(HTrue);

HRec = 20*log10(abs(PWE.P2(:,1:end/2)));
% HRec = HRec - max(HRec);

figure, hold on
plot(PWE.f,HRec)
plot(Data.f,HTrue), grid on
xlabel('Frequency in Hz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
legend('Reconstruction','True')
% axis([0 PWE.f(end) -40 0])
applyAxisProperties(gca)
applyLegendProperties(gcf)

%% Microphone
setupPlot(Data,true,false,false);
scatter3(Data.Ref.pos(Mic,1),Data.Ref.pos(Mic,2),Data.Ref.pos(Mic,3),150,'filled')

%%
figure
plot(Plot.t*1e3,PWE.h(Plot.N(1):Plot.N(2)-1,:)/max(abs(PWE.h(Plot.N(1):Plot.N(2)-1,:)))), hold on
plot(Plot.t*1e3,Data.Ref.h(Plot.N(1):Plot.N(2)-1,Mic)/max(abs(Data.Ref.h(Plot.N(1):Plot.N(2)-1,Mic))))
legend('Reconstruction','True')

%%
figure
% plot(PWE.h/max(PWE.h(:,Mic))), hold on
% plot(Data.Ref.h(:,Mic)/max(Data.Ref.h(:,Mic)))
plot(PWE.h), hold on
plot(Data.Ref.h(:,Mic))
legend('Reconstruction','True')

%% Plot
figure
s = surf(Data.Ref.pos(:,1),Plot.t*1e3,PWE.h(Plot.N(1):Plot.N(2)-1,:));
set(s,'edgecolor','none')
xlabel('x in m'), ylabel('Time in ms')
colormap hot
view(2)
c = colorbar;
caxis([-0.04 0.04])
applyColorbarProperties(c,'Room Impulse Response in Pa/V')
applyAxisProperties(gca)




