function Data = dataAcquisition(Data)
%Data = dataAcquisition(Data) Reads the corresponding data of the given
%loudspeaker.
%   Input:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: March 2022

% Source
Data.Source.pos = h5read(Data.loudspeaker,'/source/position');

% Line 1
Data.Line1.pos = h5read(Data.loudspeaker,'/dataset_lin2/position');
Data.Line1.h = h5read(Data.loudspeaker,'/dataset_lin2/impulse_response');
Data.Line1.n = h5read(Data.loudspeaker,'/dataset_lin2/noise');

% Line 2
Data.Line2.pos = h5read(Data.loudspeaker,'/dataset_lin1/position');
Data.Line2.h = h5read(Data.loudspeaker,'/dataset_lin1/impulse_response');
Data.Line2.n = h5read(Data.loudspeaker,'/dataset_lin1/noise');

% Line 3
Data.Line3.pos = h5read(Data.loudspeaker,'/dataset_lin3/position');
Data.Line3.h = h5read(Data.loudspeaker,'/dataset_lin3/impulse_response');
Data.Line3.n = h5read(Data.loudspeaker,'/dataset_lin3/noise');

% Sphere
Data.Sph.pos = h5read(Data.loudspeaker,'/dataset_sph1/position');
Data.Sph.h = h5read(Data.loudspeaker,'/dataset_sph1/impulse_response');
Data.Sph.p = h5read(Data.loudspeaker,'/dataset_sph1/pressure');
Data.Sph.n = h5read(Data.loudspeaker,'/dataset_sph1/noise');
Data.Sph.R0 = mean(Data.Sph.pos,1);

disp('Reading data... OK')

end

