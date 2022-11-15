function KRR = kernelRidgeRegression(Data,t,f,R,plotFlag)
%KRR = kernelRidgeRegression(Data,t,f,R,plotFlag) Reconstruct the late
%reverberation of the RIR assuming a dense sound field and therefore
%applying the Bessel kernel for a ridge regression.
%   Input:
%       - Data      : raw data. Structure
%       - t         : time vector. 1 x Nt
%       - f         : frequency span. 1 x Nf
%       - R         : reconstruction points. 3 x Nr
%       - plotFlag  : 'true' to plot PW reconstruction
%                     'false' (Default value)
%   Output:
%       - KRR       : kernel ridge regression in frequency domain. Structure
%
% Author: Antonio Figueroa Durán
% Date: July 2022

%% ERROR HANDLING
if nargin < 4, error('kernelRidgeRegression Error: Not enough input parameters.'), end
if nargin < 5, plotFlag = false; end

%% MAIN CODE
% Array data:
%  d = 50 cm
%  DeltaT = d/c = 1.4461 ms
%  NDeltaT = DeltaT*Fs = 69.4 == 70 samples

% Pre zero-padding
M = size(Data.InnSph.pos,1);
overlapp = 25;      % samples
Ny = size(Data.InnSph.h(ismember(Data.t,t),:),1);

p = vertcat(zeros(overlapp-1,M),Data.InnSph.h(ismember(Data.t,t),:));
KRR.Nt = size(p,1);

% Overlapp-Save Parameters
[Nblock,~] = optimalN(PWE.Nt,overlapp);
% Nblock = 75;
L = Nblock-overlapp+1;      % Convolution order

nFrames = ceil(KRR.Nt/L);

% Zero-pad at the end of the signal (last frame)
p = vertcat(p,zeros(nFrames*L-Ny,M));

% Time vector
t_padded = horzcat(zeros(1,overlapp-1),Data.t);
t_padded = t_padded(find(t_padded==t(1),1,'last')-(overlapp-1)+(1:size(p,1))-1);

% Coefficient estimation
Nr = size(R,2);
KRR.x = nan(KRR.Nt,Nr);

% Gram matrices
d_MM = nan(M);
d_MNr = nan(M,Nr);
for ii = 1:M
    d_MM(ii,:) = vecnorm(repmat(Data.InnSph.pos(ii,:),M,1)-Data.InnSph.pos,2,2).';
    d_MNr(ii,:) = vecnorm(repmat(Data.InnSph.pos(ii,:),Nr,1)-R.',2,2).';
end

for ii = 0:nFrames-1
    pm = p((1:L+overlapp-1) + ii*L,:);
    tm = t_padded((1:L+overlapp-1) + ii*L);
    
    % Frequency domain analysis
    
    % Multiply by k
    K_MM = sinc(d_MM);
    K_MNr = sinc(d_MNr).';
    
    % Kernel ridge regression (try lambda = 0 first, ask Efren )
    
    
    
    % Output: xm
    
    KRR.x(ii*L+1:(ii+1)*L,:) = xm(overlapp:overlapp+L-1,:);
    
    disp(strcat("Frame number... ",string(ii+1),"/",string(nFrames)))
end
clear H

KRR.x = KRR.x(1:Ny,:);
KRR.x = KRR.x(:);

% RIR Reconstruction
KRR.h = PWTimeReconstruction(Data,KRR.x,t,N,R);

%% PLOT
% Validation vs reference line
if plotFlag
    if Nr == 1
        Mic = find(ismember(Data.Ref.pos,R','rows'));
        
        % Mic position
        setupPlot(Data,true,false,false);
        scatter3(R(1),R(2),R(3),150,'filled')

        % Reconstruction: impulse response
        figure, hold on
        plot(t*1e3,Data.Ref.h(ismember(Data.t,t),Mic))
        plot(t*1e3,KRR.h)
        xlabel('Time in ms'), ylabel('Normalised RIR'), grid on
        legend('True','Reconstruction')
        applyAxisProperties(gca)
        applyLegendProperties(gcf)
    else
        % Reconstruction: reference line RIR
        figure
        s = surf(Data.Ref.pos(:,1),t*1e3,KRR.h);
        set(s,'edgecolor','none')
        xlabel('x in m'), ylabel('Time in ms')
        colormap hot
        view(2)
        c = colorbar;
        caxis([-0.04 0.04])
        applyColorbarProperties(c,'Room Impulse Response in Pa/V')
        applyAxisProperties(gca)
    end
end
end

