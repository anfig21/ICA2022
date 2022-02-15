function x = somp(H, p, K)
%x = somp(CA, b, k) Solve the sparse coding problem via SOMP
%   min_x sum(||p - CAx||_2^2) s.t. ||x||_0 <= K
%   
%   Input:
%       - H         : dictionary. M x N x Nf
%       - p         : measurements. M x Nf
%       - K         : number of possible sources. Scalar
%   Output:
%       - x         : estimated beamformer output. N x L
%
% Author: Antonio Figueroa Durán
% Date: February 2022

H = H./vecnorm(H,2,2);  % Normalise columns
[~,N,Nf] = size(H);     % Dictionary dimentions

s = zeros(K,1);         % Support
r = p;                  % Initialise residual

for ik = 1:K            % Possible sources
    x = zeros(N,Nf);
    c = zeros(N,Nf);
    for jj = 1:Nf        % Frequency
        c(:,jj) = squeeze(H(:,:,jj))'*r(:,jj);    % (N,M)x(M,1) = (N,1)
    end
    
    [~,i] = max(sum(abs(c),2));
    s(ik) = i;
    
    for jj = 1:Nf
        CAs = squeeze(H(:,s(1:ik),jj));
        x(s(1:ik),jj) = pinv(CAs)*p(:,jj);
        r(:,jj) = p(:,jj) - CAs*x(s(1:ik),jj);
    end
end

end

