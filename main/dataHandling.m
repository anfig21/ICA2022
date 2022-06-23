function Data = dataHandling(Data)
%Data = dataHandling(Data) Processes the corresponding data of the given
%loudspeaker.
%   Input:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: June 2022

% Pre-processing
Data.Ref.pos = vertcat(Data.Line1.pos,Data.Line2.pos,Data.Line3.pos);
Data.Ref.h = horzcat(Data.Line1.h,Data.Line2.h,Data.Line3.h);
Data.Ref.n = horzcat(Data.Line1.n,Data.Line2.n,Data.Line3.n);

% Recwin up to T = 1s
Data.Ref.h = Data.Ref.h(1:Data.Nsamples,:);
Data.Ref.n = Data.Ref.n(1:2*Data.Nsamples,:);
Data.Sph.h = Data.Sph.h(1:Data.Nsamples,:);
Data.Sph.n = Data.Sph.n(1:2*Data.Nsamples,:);
Data.Sph.p = Data.Sph.p(1:2*Data.Nsamples,:);

% Inner Sphere (156 Microphones: samples 155-end)
Idx = 155:size(Data.Sph.pos,1);
Data.InnSph.M = length(Idx);
Data.InnSph.pos = Data.Sph.pos(Idx,:);
Data.InnSph.h = Data.Sph.h(:,Idx);
Data.InnSph.p = Data.Sph.p(:,Idx);

% Frequency domain
[Data.InnSph.H,~] = fftUniBi(Data.InnSph.h);
[Data.InnSph.P,~] = fftUniBi(Data.InnSph.p);
[Data.Sph.N,~] = fftUniBi(Data.Sph.n);
[Data.Ref.H,~] = fftUniBi(Data.Ref.h);

end

