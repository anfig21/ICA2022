mic = 220;
figure
subplot(211), hold on
plot(Plot.t*1e3,Data.Ref.h(Plot.N(1):Plot.N(2),mic))
plot(Plot.t*1e3,Rec.Direct.h(Plot.N(1):Plot.N(2),mic)), grid on
xlabel('Time in ms'), ylabel('RIR in Pa/V')
legend('True','Rec'), title('Interpolation')
applyAxisProperties(gca)
applyLegendProperties(gcf)

mic = 120;
subplot(212), hold on
plot(Plot.t*1e3,Data.Ref.h(Plot.N(1):Plot.N(2),mic))
plot(Plot.t*1e3,Rec.Direct.h(Plot.N(1):Plot.N(2),mic)), grid on
xlabel('Time in ms'), ylabel('RIR in Pa/V')
legend('True','Rec'), title('Extrapolation')
applyAxisProperties(gca)
applyLegendProperties(gcf)