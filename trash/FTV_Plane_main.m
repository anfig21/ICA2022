%% Range Estimation via Fused Total Generalised Variation - Laplacian/Stencil 3x3
% Dictionary
Dict.Grid.Res = [.25 .25];
[Dict.Grid.H,Dict.Grid.r,Dict.Grid.N,Dict.Grid.Nx,Dict.Grid.Ny] = ...
    dictionaryGrid(Dict.f,Data.InnSph.pos.',Data.Sph.R0.',Data.D,...
                    Dict.Grid.Res(1),Dict.Grid.Res(2));

% Total Variation (Over a plane)
TV = dir2D_FTV(Data,Direct,Dict,'true');