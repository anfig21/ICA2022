function Early = earlyWindow(Data,Early)
%Early = earlyWindow(Data,Early) Applies Hanning window to sound field
%composed of the early reflections. Calculates the spectrum and the noise
%norm.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early reflections. Structure with time limits
%   Output:
%       - Early     : early reflections. Structure
%
% Author: Antonio Figueroa Durán
% Date: April 2022

tHann = 1e-3;
NHann = Data.Fs*tHann/2;
h = hann(2*NHann);
hHalf1 = h(1:end/2);
hHalf2 = h(end/2+1:end);

Early.N = Data.Fs*Early.T;
Early.Nsamples = Early.N(2)-Early.N(1);

% Windowing - 1/2 Hanning (hann) window on each side
Early.w = vertcat(zeros(Early.N(1),Data.InnSph.M),...
    repmat(hHalf1,1,Data.InnSph.M),...
    ones(Early.Nsamples-2*NHann,Data.InnSph.M),...
    repmat(hHalf2,1,Data.InnSph.M),...
    zeros(Data.Nsamples-Early.Nsamples-Early.N(1),Data.InnSph.M));
Early.InnSph.h = Early.w.*Data.InnSph.h;
Early.InnSph.n = Early.w(:,size(Data.Sph.n,2)).*Data.Sph.n(1:2:end,:);

% Frequency Domain
Early.InnSph.H = fft(Early.InnSph.h,2*Data.Nsamples)/Data.Nsamples;
Early.InnSph.H = [Early.InnSph.H(1,:); 2*Early.InnSph.H(2:end/2,:)];
Early.InnSph.N = fft(Early.InnSph.n,2*Data.Nsamples)/Data.Nsamples;
Early.InnSph.N = [Early.InnSph.N(1,:); 2*Early.InnSph.N(2:end/2,:)];

% Noise norm
Early.InnSph.Nnorm = mean(abs(Early.InnSph.N),2);

% Plot frequency response of mics [40, 59, 69]
% figure, plot(Data.f*1e-3,20*log10(abs(Early.InnSph.H(:,[40 59 69])))), grid on
% xlabel('Frequency in kHz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
% legend('Mic 40','Mic 59','Mic 69')
% applyAxisProperties(gca)
% applyLegendProperties(gcf)

% Plot RIR Time domain
% figure
% subplot(211), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Data.InnSph.h)
% xlim([5 35]), title('RIR')
% applyAxisProperties(gca)
% subplot(212), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Early.InnSph.h)
% xlim([5 35]), title('Early Reflections')
% applyAxisProperties(gca)

end

