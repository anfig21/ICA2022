function Rec = reconstructReflection(Data,Rec,R0,L,res,RecMethod,plotFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% ERROR HANDLING
% plotFlag default value
if nargin < 7, plotFlag = false;
elseif nargin < 6, error('reconstructReflection Error: Not enough input parameters.'), end

%% MAIN CODE
% Point sources - cubic grid
l = -L/2:res:L/2;
[X,Y,Z] = meshgrid(l,l,l);
rs = [X(:) Y(:) Z(:)]+R0(:)';

% figure, scatter3(rs(:,1),rs(:,2),rs(:,3)), axis equal

M = size(rs,1);                 % Number of point sources
N = size(Data.InnSph.pos,1);    % Number of microphones
Nf = length(Rec.f);             % Number of frequency bins

k = (2*pi*Rec.f)/Data.c;        % Propagation vector

% Euclidean distance matrix
d = nan(N,M);       % point sources x mics
for nn = 1:N
    d(nn,:) = vecnorm(repmat(Data.InnSph.pos(nn,:),M,1)-rs,2,2);
end

% Coefficient estimation
NoiseMargin = 0;
x = nan(M,Nf);
c = parcluster('local');
parpool(c,4);
% for ii = 1:Nf
parfor ii = 1:Nf
    H = (1./d).*exp(-1i*d*k(ii));       % Dictionary (point sources)
    %H = exp(-1i*d*k(ii));
    pii = Rec.InnSph.H(Data.f==Rec.f(ii),:).';
    
    switch RecMethod
        case 'RLS'
            % Regularised Least-Squares
            [x(:,ii),~] = reguLeastSquares(H,pii);
        case 'CS'
            Nnorm = 10^(NoiseMargin/20)*Rec.InnSph.Nnorm(Data.f==Rec.f(ii));
            
%             [U,s,V] = csvd(H);
%             [lambda,~,~,~] = l_curve(U,s,pii,'Tikh',[],[],false);
            
            x(:,ii) = compressiveSensing(H,M,pii,Nnorm);
            
            disp(strcat("Reconstructing reflection... ",string(ii-1),'/',string(Nf-1)," Hz"))            
        otherwise
            error('reconstructReflection Error: Regularisation method not valid.')
    end
end
Rec.x = x;
delete(gcp('nocreate'))

% Reconstruction
rr = Data.Ref.pos;
R = size(rr,1);

d = nan(M,R);       % point sources x mics
for nn = 1:M
    d(nn,:) = vecnorm(repmat(rs(nn,:),R,1)-rr,2,2);
end

% Dictionary
Rec.P = zeros(R,length(Data.f));
for ii = 1:Nf
    %H = (1./d).*exp(-1i*d*k(ii));       % Dictionary (point sources)
    H = exp(-1i*d*k(ii));
    Rec.P(:,ii) = H.'*x(:,ii);
end

% Double-sided spectrum
P2 = [real(Rec.P(:,1)) Rec.P(:,2:end)/2];
P2 = [P2 flip(conj(P2(:,2:end)),2)];
Rec.h = ifft(P2*Data.Nsamples,[],2,'symmetric').';

%% PLOT
% Validation vs reference line
if plotFlag
    % Reconstruction: reference line RIR
    figure, hold on
    s = pcolor(Data.Ref.pos(:,1),Plot.t*1e3,Rec.h(Plot.N(1):Plot.N(2)-1,:));
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

