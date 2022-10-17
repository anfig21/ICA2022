h = Data.Ref.h(:,end/2);

figure
plot(Data.t*1e3,h)
xlim([0 25])

fd = h(1:end-1)-h(2:end);

figure, grid on
plot(Data.t(1:end-1)*1e3,abs(fd))
xlim([0 25])

%%
R0 = Data.Source.pos;
r0 = Data.Ref.pos(end/2,:);

T0 = vecnorm(R0-r0)/Data.c;             % Direct. T = 5.9494 ms

R1 = [R0(1) R0(2) -R0(3)];
T1 = vecnorm(R1-r0)/Data.c;             % Floor. T = 9.4784 ms

R2 = [R0(1) -R0(2) R0(3)];
T2 = vecnorm(R2-r0)/Data.c;             % Back. T = 17.6203 ms

R3 = [R0(1) R0(2) 2*Data.D(3)-R0(3)];   % Ceiling. T = 11.4793 ms
T3 = vecnorm(R3-r0)/Data.c;

%%
% Delays
t0 = 6e-3;
t1 = 9.58333e-3;
t2 = 11.4375e-3;
t3 = 17.6667e-3;

% Range
s0 = vecnorm(Direct.Range.Mode-r0(:));
s1 = s0+(t1-t0)*Data.c;
s2 = s0+(t2-t0)*Data.c;
s3 = s0+(t3-t0)*Data.c;

%%
SlopeThreshold = 0;
AmpThreshold = 5e-3;
smoothwidth = 0;
peakgroup = 1;
smoothtype = 1;

P=findpeaksx(Data.t,h,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);

figure, grid on
plot(Data.t*1e3,h), hold on
for ii = 1:size(P,1), xline(P(ii,2)*1e3), end
xlim([0 25])
