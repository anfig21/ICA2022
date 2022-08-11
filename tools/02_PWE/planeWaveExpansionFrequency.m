function PWE = planeWaveExpansionFrequency(Data,PWE,Plot,N,f,R,PWEMethod,plotFlag)
%PWE = planeWaveExpansionFrequency(Data,PWE,Plot,N,f,PWEMethod,plotFlag)
%Description of the function
%   Input:
%       - Data      : raw data. Structure
%       - PWE       : windowed RIR. Structure
%       - Plot      : plot parameters. Structure
%       - N         : number of plane waves. Integer
%       - f         : frequency span. 1 x Nf
%       - R         : reconstruction points. 3 x Nr
%       - PWEMethod : specifies the optimisation method:
%                       'RLS' - Regularised Least Squares
%                       'CS' - Compressive Sensing
%       - plotFlag  : 'true' to plot PW reconstruction
%                     'false' (Default value)
%   Output:
%       - PWE       : plane wave expansion in frequency domain. Structure
%
% Author: Antonio Figueroa Durán
% Date: June 2022

%% ERROR HANDLING
if nargin < 7, error('planeWaveExpansionFrequency Error: Not enough input parameters.'), end
if N <= 0, error('planeWaveExpansionFrequency Error: The number of plane waves must be a positive integer.'), end
if nargin < 8, plotFlag = false; end

%% MAIN CODE
% Coefficient estimation
Nf = length(f);
PWE.x = nan(N,Nf);

for ii = 1:Nf
    [H,~] = dictionary(Data.c,f(ii),Data.InnSph.pos',N);
    
    switch PWEMethod
        case 'RLS'
            % PWE via Regularised Least-Squares
            [PWE.x(:,ii),~] = reguLeastSquares(H,PWE.InnSph.H(Data.f==f(ii),:).');
        case 'CS'
            % PWE via Compressive Sensing
            Nnorm = 10^(10/20)*PWE.InnSph.Nnorm(Data.f==f(ii));
            p = PWE.InnSph.H(Data.f==f(ii),:).';
            cvx_begin quiet
            cvx_precision high
                variable x(N) complex;
                minimize norm(x,1);
            subject to
                norm((H*x-p),2) <= Nnorm;
            cvx_end
            PWE.x(:,ii) = x;
        otherwise
            error('planeWaveExpansionFrequency Error: PWE method not valid.')
    end
    disp(strcat("Estimating coefficients... ",string(round(ii/Nf*100,2)),"%"))
end

% RIR Reconstruction
Nr = size(R,2);
PWE.P = zeros(Nr,length(Data.f));

for ii = 1:Nf
    [HR,~] = dictionary(Data.c,f(ii),R,N);
    PWE.P(:,Data.f==f(ii)) = squeeze(HR)*PWE.x(:,ii);
    disp(strcat("Reconstructing sound field... ",string(round(ii/Nf*100,2)),"%"))
end

% Double-sided spectrum
PWE.P2 = [real(PWE.P(:,1)) PWE.P(:,2:end)/2];
PWE.P2 = [PWE.P2 flip(conj(PWE.P2(:,2:end)),2)];
PWE.h = ifft(PWE.P2*Data.Nsamples,[],2,'symmetric').';

%% PLOT
% Validation vs reference line
if plotFlag
    if Nr == 1
        Mic = find(ismember(Data.Ref.pos,R','rows'));
        
        % Mic position
        setupPlot(Data,true,false,false);
        scatter3(R(1),R(2),R(3),150,'filled')

        % Reconstruction: frequency response
        HTrue = 20*log10(abs(Data.Ref.H(:,Mic)));
        HRec = 20*log10(abs(PWE.P));

        figure, hold on
        plot(Data.f,HTrue), grid on
        plot(Data.f,HRec)
        xlabel('Frequency in Hz'), ylabel('Magnitude $|H(j\omega)|$ in dB')
        legend('True','Reconstruction')
        axis([0 PWE.f(end) max([HTrue.' HRec])-40 inf])
        applyAxisProperties(gca)
        applyLegendProperties(gcf)

        % Reconstruction: impulse response
        figure, hold on
        plot(Plot.t*1e3,Data.Ref.h(Plot.N(1):Plot.N(2)-1,Mic)/max(abs(Data.Ref.h(Plot.N(1):Plot.N(2)-1,Mic))))
        plot(Plot.t*1e3,PWE.h(Plot.N(1):Plot.N(2)-1,:)/max(abs(Data.Ref.h(Plot.N(1):Plot.N(2)-1,Mic))))
        xlabel('Time in ms'), ylabel('Normalised RIR'), grid on
        legend('True','Reconstruction')
        applyAxisProperties(gca)
        applyLegendProperties(gcf)
    else
        % Reconstruction: reference line RIR
        figure, hold on
        s = pcolor(Data.Ref.pos(:,1),Plot.t*1e3,PWE.h(Plot.N(1):Plot.N(2)-1,:));
        set(s,'edgecolor','none')
        xlabel('x in m'), ylabel('Time in ms')
        colormap hot
        c = colorbar;
        caxis([-0.04 0.04])
        xline(min(Data.InnSph.pos(:,1))), xline(max(Data.InnSph.pos(:,1)))
        applyColorbarProperties(c,'Room Impulse Response in Pa/V')
        applyAxisProperties(gca)
    end
end
end

