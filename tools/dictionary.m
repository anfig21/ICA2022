function [H,uk] = dictionary(f,r0,N)
%[H,uk] = dictionary(f,N) Obtains the H matrix for a sound field at a
% position r0 comprised of N plane waves with propagation vector sampled
% along the surface of a unit sphere.
%   Input:
%       - f         : frequency (Hz). Nf x 1
%       - r0        : measurement point. 3 x M
%       - N         : number of plane waves. Scalar
%   Output:
%       - H         : dictionary. M x N x Nf
%       - uk        : propagation unit vector. 3 x N
%
% Author: Antonio Figueroa Durán
% Date: February 2022

c = 343;
Nf = length(f);
[~,M] = size(r0);

% Spatial sampling
uk = RandSampleSphere(N,'uniform').';

% Propagation vector
k = c./(2*pi*f);
kv = nan(3,N,Nf);
for ii = 1:Nf, kv(:,:,ii) = k(ii)*uk; end

% Dictionary
H = nan(M,N,Nf);
c = waitbar(0,'Loading...0%','Name','Building dictionary...');
for ii = 1:Nf
    H(:,ii) = exp(-1i*r0'*kv(:,:,ii));
    waitbar(ii/Nf,c,sprintf('Loading... %.2f%%',100*ii/Nf));
end
delete(c)

% Normalisation
H = H./vecnorm(H,2,2);

end