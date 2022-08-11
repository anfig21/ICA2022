
%% Coefficients
N = 500;

mics = Data.Line2.pos(1,:);
t = Data.t(Data.t >= 0e-3 & Data.t <= 10e-3);
L = length(t);
M = size(mics,1);

Idx = 1;
x = zeros(L,N);
x(:,Idx) = 2;
x = x(:);

%% Synthesis
[HR,uk] = dictionaryPWTime(Data.c,t,Data.Fs,mics.',Data.Sph.R0.',N);
h = HR*x;

phi = rad2deg(atan(uk(2,Idx)/uk(1,Idx)));
disp(string(phi))

clear HR
%% Plot
h_m = reshape(h,[L M]);

figure('units','normalized','outerposition',[0 0 1 1])
s = surf(mics(:,1),t*1e3,h_m);
set(s,'edgecolor','none')
xlabel('x in m'), ylabel('Time in ms')
colormap hot
view(2)
c = colorbar;
% caxis([-0.04 0.04])
applyColorbarProperties(c,'Room Impulse Response in Pa/V')
applyAxisProperties(gca)

%% Plot sources
figure
scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled'), hold on
scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3))
scatter3(Data.InnSph.pos(:,1),Data.InnSph.pos(:,2),Data.InnSph.pos(:,3))
scatter3(uk(1,:)+mics(M/2,1),uk(2,:)+mics(M/2,2),uk(3,:)+mics(M/2,3),200,'filled')
drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
legend('Source','Reference Line','Spherical Array','Sources')
applyAxisProperties(gca)
applyLegendProperties(gcf)
