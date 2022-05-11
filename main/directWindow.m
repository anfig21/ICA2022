function Direct = directWindow(Data,Direct)
%Direct = directWindow(Data, Direct) Applies Hanning window to direct field
%from the RIR. Calculates the spectrum and the noise norm.
%   Input:
%       - Data      : raw data. Structure
%       - Direct    : direct sound field. Structure with time limits
%   Output:
%       - Direct    : direct field data. Structure
%
% Author: Antonio Figueroa Durán
% Date: April 2022

tHann = 1e-3;
NHann = Data.Fs*tHann/2;
h = hann(2*NHann);
hHalf = h(end/2+1:end);

Direct.Nsamples = Data.Fs*Direct.T;

% Windowing - 1/2 Hanning (hann) window on each side
Direct.w = vertcat(ones(Direct.Nsamples-NHann,Data.InnSph.M),...
    repmat(hHalf,1,Data.InnSph.M),...
    zeros(Data.Nsamples-Direct.Nsamples,Data.InnSph.M));

Direct.InnSph.h = Direct.w.*Data.InnSph.h;
Direct.InnSph.n = Direct.w(:,size(Data.Sph.n,2)).*Data.Sph.n(1:2:end,:);

% Frequency Domain
Direct.InnSph.H = fft(Direct.InnSph.h,2*Data.Nsamples)/Data.Nsamples;
Direct.InnSph.H = [Direct.InnSph.H(1,:); 2*Direct.InnSph.H(2:end/2,:)];
Direct.InnSph.N = fft(Direct.InnSph.n,2*Data.Nsamples)/Data.Nsamples;
Direct.InnSph.N = [Direct.InnSph.N(1,:); 2*Direct.InnSph.N(2:end/2,:)];

% Noise norm
Direct.InnSph.Nnorm = mean(abs(Direct.InnSph.N),2);

% True DOA
Direct.TrueDOA = Data.Source.pos-Data.Sph.R0;
Direct.TrueDOA = Direct.TrueDOA/vecnorm(Direct.TrueDOA);

% Plot frequency response of mics [40, 59, 69]
% figure, plot(Data.f*1e-3,20*log10(abs(Direct.InnSph.H(:,[40 59 69])))), grid on
% xlabel('Frequency in kHz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
% legend('Mic 40','Mic 59','Mic 69')
% applyAxisProperties(gca)
% applyLegendProperties(gcf)

% Plot RIR Time domain
figure
subplot(211), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Data.InnSph.h), grid on
xlim([5 35]), title('RIR'),
applyAxisProperties(gca)
subplot(212), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Direct.InnSph.h), grid on
xlim([5 35]), title('Direct Sound'), xlabel('Time in ms')
applyAxisProperties(gca)

disp('Direct sound: Windowing... OK')

end

