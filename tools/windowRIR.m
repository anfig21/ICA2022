function Structure = windowRIR(Data,Tini,Tfin,plotFlag)
%Structure = windowRIR(Data,Tini,Tfin,plotFlag) Applies Hanning window to
%RIR between time given by Tini and Tfin. Obtains the frequency response
%and the noise spectrum. Estimates the noise power after windowing.
%   Input:
%       - Data      : raw data. Structure
%       - Tini      : initial time. Scalar
%       - Tfin      : end time. Scalar
%       - plotFlag  : 'true'  - RIR in time domain
%                     'false' - Default value
%   Output:
%       - Structure : structure with windowed data
%
% Author: Antonio Figueroa Durán
% Date: June 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('windowRIR Error: Not enough input parameters.'), end

%% MAIN CODE
tHann = 1e-3;
NHann = Data.Fs*tHann/2;
h = hann(2*NHann);
hHalf1 = h(1:end/2);
hHalf2 = h(end/2+1:end);

Structure.T = [Tini Tfin];
Structure.N = floor(Data.Fs*Structure.T);
Structure.Nsamples = Structure.N(2)-Structure.N(1);

% Windowing - 1/2 Hanning (hann) window on each side
Structure.w = vertcat(zeros(Structure.N(1),Data.InnSph.M),...
    repmat(hHalf1,1,Data.InnSph.M),...
    ones(Structure.Nsamples-2*NHann,Data.InnSph.M),...
    repmat(hHalf2,1,Data.InnSph.M),...
    zeros(Data.Nsamples-Structure.Nsamples-Structure.N(1),Data.InnSph.M));
Structure.InnSph.h = Structure.w.*Data.InnSph.h;
Structure.InnSph.n = Structure.w(:,size(Data.Sph.n,2)).*Data.Sph.n(1:2:end,:);

% Frequency Domain
Structure.InnSph.H = fft(Structure.InnSph.h,2*Data.Nsamples)/Data.Nsamples;
Structure.InnSph.H = [Structure.InnSph.H(1,:); 2*Structure.InnSph.H(2:end/2,:)];
Structure.InnSph.N = fft(Structure.InnSph.n,2*Data.Nsamples)/Data.Nsamples;
Structure.InnSph.N = [Structure.InnSph.N(1,:); 2*Structure.InnSph.N(2:end/2,:)];

% Noise norm
Structure.InnSph.Nnorm = mean(abs(Structure.InnSph.N),2);

% Plot frequency response of mics [40, 59, 69]
% figure, plot(Data.f*1e-3,20*log10(abs(Early.InnSph.H(:,[40 59 69])))), grid on
% xlabel('Frequency in kHz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
% legend('Mic 40','Mic 59','Mic 69')
% applyAxisProperties(gca)
% applyLegendProperties(gcf)

% Plot RIR Time domain
if plotFlag
    figure
    subplot(211), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Data.InnSph.h), grid on
    xlim([5 35]), title('RIR')
    applyAxisProperties(gca)
    subplot(212), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Structure.InnSph.h), grid on
    xlim([5 35]), title('Early Reflections'), xlabel('Time in ms')
    applyAxisProperties(gca)
end

end

