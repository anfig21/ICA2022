function [] = plotFreqResponse(Data)
%plotFreqResponse(Data) Plots the measurement setup, the magnitude of the
%transfer function and the frequency response of signal and noise.
%   Input:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: March 2022

% Line plot
figure
scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on
scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3))
axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
legend('Reference Line','Spherical Array','Source')
applyAxisProperties(gca)
applyLegendProperties(gcf)

% Frequency response
figure, plot(Data.f*1e-3,20*log10(abs(Data.InnSph.H(:,[40 59 69])))), grid on
xlabel('Frequency in kHz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
legend('Mic 40','Mic 59','Mic 69')
applyAxisProperties(gca)
applyLegendProperties(gcf)

figure
subplot(221)
plot(Data.f*1e-3,20*log10(abs(Data.InnSph.P(:,[40 59 69]))/20e-6)), grid on
xlabel('Frequency in kHz'), ylabel('$|P(j\omega)|$ in dB SPL')
legend('Mic 40','Mic 59','Mic 69')
applyAxisProperties(gca)
applyLegendProperties(gcf)
ylim([-5 70])

subplot(222)
plot(Data.f,20*log10(abs(Data.InnSph.P(:,[40 59 69]))/20e-6)), grid on
xlabel('Frequency in Hz'), ylabel('$|P(j\omega)|$ in dB SPL')
legend('Mic 40','Mic 59','Mic 69')
applyAxisProperties(gca)
applyLegendProperties(gcf)
axis([0 4e3 -5 70])

subplot(223)
plot(Data.f*1e-3,20*log10(abs(Data.Sph.N)/20e-6)), grid on
xlabel('Frequency in kHz'), ylabel('$|N(j\omega)|$ in dB SPL')
applyAxisProperties(gca)
ylim([-5 70])

subplot(224)
plot(Data.f,20*log10(abs(Data.Sph.N)/20e-6)), grid on
xlabel('Frequency in Hz'), ylabel('$|N(j\omega)|$ in dB SPL')
applyAxisProperties(gca)
axis([0 4e3 -5 70])

end

