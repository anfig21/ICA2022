function Range = dirRange_TV(Data,Direct,Dict)
%Range = dirRange_TV(Data,Direct,Dict) Applies Total Variation to the
%Direct Sound Field to obtain the position of the source.
%   Input:
%       - Data          : raw data. Structure
%       - Direct        : direct sound field. Structure
%       - Dict          : dictionary of spherical waves. Structure
%   Output:
%       - Range        : Range estimation via TV. Structure
%
% Author: Antonio Figueroa Durán
% Date: May 2022

%% ERROR HANDLING
if nargin < 3, error('dirRange_TV Error: Not enough input parameters.'), end

%% MAIN CODE
Range.Est = nan(3,length(Dict.f));
Range.Idx = nan(1,length(Dict.f));
NoiseMargin = 10;           % dB

% Stencil
% st = 3; mask = [1; -1; 0];        % First order
% st = 3; mask = [1; -2; 1];        % Second order
st = 3; mask = [0.5; -1; 0.5];    % Weighted Second order
% st = 5; mask = [0.5; 1; -3; 1; 0.5];  % Fourth order

Mask = padarray(mask, Dict.N-st,'post');
M = circshift(Mask(:),Dict.N-(st-1)/2);
D = zeros(Dict.N);
for ii=1:Dict.N
    D(ii,:)=circshift(M.',[0 ii-1]);
end

c = waitbar(0,'Loading...0\%','Name','dirRange_TV: CVX across frequencies...');
for ii = 1:length(Dict.f)
    Nnorm = 10^(NoiseMargin/20)*Direct.InnSph.Nnorm(Data.f==Dict.f(ii));
    Hii = squeeze(Dict.H(:,:,ii));
    pii = Direct.InnSph.H(Data.f==Dict.f(ii),:).';
%     pii = unwrap(angle(pii));
    
    % CVX Formulation
    cvx_begin quiet
    cvx_precision high
    variable x(Dict.N) complex;
    minimize norm(D*x,1);
    subject to
    norm((Hii*x-pii),2) <= Nnorm;
    cvx_end
    
    [~,Range.Idx(ii)] = max(abs(x));
    Range.Est(:,ii) = Dict.r(:,Range.Idx(ii));
    
    waitbar(ii/length(Dict.f),c,strcat("Loading... ",...
        string(round(100*ii/length(Dict.f),2)),"\,\%"));
end
delete(c)

disp('Direct sound: RANGE - TV... OK')

end

