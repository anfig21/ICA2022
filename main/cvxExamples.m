%% CVX EXAMPLES
% Efren Fernandez Grande

%% SOLVE SYSTEM OF LINEAR EQUATIONS
%Conventional regularised inversion:
% figure; % A=(1j*1.2*343*k)*(Ghp); % [U,s,V]=csvd(A); % [lambda] =l_curve(U,s,ph);%[q_l2] = tikhonov(U,s,V,ph,lambda);
 
A=(1j*1.2*343*k)*(Ghp);
%Estimate relavite noise level
snr=SNR; rnl = 10^-(snr/20)*norm(ph);
n0=nx0*ny0;
 
%% CVX solution l2:
tic
cvx_begin
cvx_quiet(true)
cvx_precision high
    variable q_l2(n0) complex;
    minimize norm(q_l2,2);
    subject to
        norm((A*q_l2-ph),2) <= rnl;
cvx_end
toc
 
%% CVX solution CS:
tic
cvx_begin
cvx_quiet(true)
cvx_precision high
    variable q_l1(n0) complex;
    minimize norm(q_l1,1);
    subject to
        norm((A*q_l1-ph),2) <= rnl;
cvx_end
% toc
 
%% ELASTIC NET - CVX solution EN:
% tic
% cvx_begin
% cvx_quiet(true)
% cvx_precision high
%     variable q_en(n0) complex;
%     minimize norm(q_en,1) + norm(q_en,2);
%     subject to
%         norm((A*q_en-ph),2) <= rnl;
% cvx_end
% toc
 
%% TOTAL VARIATION L-1
D=-1*eye(ny0)+circshift(eye(ny0),-1);
Dx = kron(D,eye(nx0));
Dy = kron(eye(ny0),D);
% D=zeros(ny0*nx0);
% for ii=1:ny0*nx0
%     D(ii,:)=circshift(M.',[0 ii-1]);
% end
Dcs=[D];
%figure;surf(D);shading flat;colormap(flipud(gray))
clear mask Mask M
tic
cvx_begin quiet
    variable q_tv(n0) complex;
%    minimize norm(D*q_tv,1)
%    minimize norm(Dx*q_tv,1) + norm(Dy*q_tv,1)
%   minimize norm(Dx*q_tv,2) + norm(Dy*q_tv,2); %works well
    minimize norm(Dx*q_tv,1) + norm(Dy*q_tv,1) %
 
     subject to
        (norm(A*q_tv-ph)) <= rnl;
cvx_end
toc
 
 
%% TOTAL GENERALIZED VARIATION 1 - x-y 5 point stencil; mask= [0,1,0; 1,-4,1; 0,1,0];
%LAPLACIAN TGV;
st=3; % stencil dimension
mask= [0,1,0; 1,-4,1; 0,1,0];%/d0^2;
%mask= [0.5,1,0.5; 1,-6,1; 0.5,1,0.5];%/d0^2; % mask= [1,1,1; 1,-8,1; 1,1,1]/d0^2;
 
Mask=padarray(mask, [ny0 nx0]-st,'post');
M=circshift(Mask(:),ny0*nx0-((st-1)*0.5*ny0+(st-1)/2));clear st
D=zeros(ny0*nx0);
for ii=1:ny0*nx0
    D(ii,:)=circshift(M.',[0 ii-1]);
end
Dcs=[D];
%figure;surf(D);shading flat;colormap(flipud(gray))
clear mask Mask M
tic
cvx_begin quiet
    variable q_tgv1(n0) complex;
%    minimize norm(D*q_tv,1)
    minimize norm(Dcs*q_tgv1,1)
     subject to
        (norm(A*q_tgv1-ph)) <= rnl;
cvx_end
toc
 
 
%% TOTAL GENERALIZED VARIATION 2 - 5 point stencil; compact stencil mask= [0.5,1,0.5; 1,-6,1; 0.5,1,0.5];
% %LAPLACIAN TGV;
% st=3; % stencil dimension
% %mask= [0,1,0; 1,-4,1; 0,1,0];%/d0^2;
% mask= [0.5,1,0.5; 1,-6,1; 0.5,1,0.5];%/d0^2; % mask= [1,1,1; 1,-8,1; 1,1,1]/d0^2;
%
% Mask=padarray(mask, [ny0 nx0]-st,'post');
% M=circshift(Mask(:),ny0*nx0-((st-1)*0.5*ny0+(st-1)/2));clear st
% D=zeros(ny0*nx0);
% for ii=1:ny0*nx0
%     D(ii,:)=circshift(M.',[0 ii-1]);
% end
% Dcs=[D];
% %figure;surf(D);shading flat;colormap(flipud(gray))
% clear mask Mask M
% tic
% cvx_begin quiet
%     variable q_tgv2(n0) complex;
% %    minimize norm(D*q_tv,1)
%     minimize norm(Dcs*q_tgv2,1)
%      subject to
%         (norm(A*q_tgv2-ph)) <= rnl;
% cvx_end
% toc
 
%% TOTAL GENERALIZED VARIATION 3 - stencil dimension 5 x 5
% %LAPLACIAN TGV;
% % st=3; % stencil dimension
% % % mask= [0,1,0; 1,-4,1; 0,1,0];%/d0^2;
% % mask= [0.5,1,0.5; 1,-6,1; 0.5,1,0.5];%/d0^2; % mask= [1,1,1; 1,-8,1; 1,1,1]/d0^2;
% st=5; % stencil dimension - nine point stencil
% %mask= [0, 0, -1/12, 0, 0; 0, 0, 4/3, 0, 0; -1/12, 4/3, -5/2, 4/3, -1/12; 0, 0, 4/3, 0, 0; 0, 0, -1/12, 0, 0];
% mask= [0, 0, 1, 0, 0; 0, 0, -16, 0, 0; 1, -16, 60, -16, 1; 0, 0, -16, 0, 0; 0, 0, 1, 0, 0]/60;
%
% % % 4th DERIVATIVE;
% % st=5; % stencil dimension
% % mask= [0,0,1,0,0; 0,0,-4,0,0; 1,-4,12,-4,1; 0,0,-4,0,0; 0,0,1,0,0]/d0^2;
%
% Mask=padarray(mask, [ny0 nx0]-st,'post');
% M=circshift(Mask(:),ny0*nx0-((st-1)*0.5*ny0+(st-1)/2));clear st
% D=zeros(ny0*nx0);
% for ii=1:ny0*nx0
%     D(ii,:)=circshift(M.',[0 ii-1]);
% end
% Dcs=[D];
% %figure;surf(D);shading flat;colormap(flipud(gray))
% clear mask Mask M
% tic
% cvx_begin quiet
%     variable q_tgv3(n0) complex;
% %    minimize norm(D*q_tv,1)
%     minimize norm(Dcs*q_tgv3,1)
%      subject to
%         (norm(A*q_tgv3-ph)) <= rnl;
% cvx_end
% toc
 
%% TOTAL GENERALIZED VARIATION 4 - 4th derivative; 9-point stencil (dimension 5 x 5); mask= [0,0,1,0,0; 0,0,-4,0,0; 1,-4,12,-4,1; 0,0,-4,0,0; 0,0,1,0,0];
% %LAPLACIAN TGV;
% % st=3; % stencil dimension
% % % mask= [0,1,0; 1,-4,1; 0,1,0];%/d0^2;
% % mask= [0.5,1,0.5; 1,-6,1; 0.5,1,0.5];%/d0^2; % mask= [1,1,1; 1,-8,1; 1,1,1]/d0^2;
% % st=5; % stencil dimension
% % mask= [0, 0, -1/12, 0, 0; 0, 0, 4/3, 0, 0; -1/12, 4/3, -5/2, 4/3, -1/12; 0, 0, 4/3, 0, 0; 0, 0, -1/12, 0, 0];
% % 4th DERIVATIVE;
% st=5; % stencil dimension
% mask= [0,0,1,0,0; 0,0,-4,0,0; 1,-4,12,-4,1; 0,0,-4,0,0; 0,0,1,0,0];%/d0^2;
% Mask=padarray(mask, [ny0 nx0]-st,'post');
% M=circshift(Mask(:),ny0*nx0-((st-1)*0.5*ny0+(st-1)/2));clear st
% D=zeros(ny0*nx0);
% for ii=1:ny0*nx0
%     D(ii,:)=circshift(M.',[0 ii-1]);
% end
% Dcs=[D];
% %figure;surf(D);shading flat;colormap(flipud(gray))
% clear mask Mask M
% tic
% cvx_begin quiet
%     variable q_tgv4(n0) complex;
% %    minimize norm(D*q_tv,1)
%     minimize norm(Dcs*q_tgv4,1)
%      subject to
%         (norm(A*q_tgv4-ph)) <= rnl;
% cvx_end
% toc
 
%% FUSED TOTAL GENERALIZED VARIATION - FUSED Laplacian; stencil dimension 3 x 3
% %First derivative, centered
% st=3; % stencil dimension
% mask= [0,-0.5,0; -0.5,0,0.5; 0,0.5,0];
 
%LAPLACIAN TGV;
st=3; % stencil dimension
% mask= [0,1,0; 1,-4,1; 0,1,0];%/d0^2;
mask= [0.5,1,0.5; 1,-6,1; 0.5,1,0.5];%/d0^2; % mask= [1,1,1; 1,-8,1; 1,1,1]/d0^2;
% st=5; % stencil dimension
% mask= [0, 0, -1/12, 0, 0; 0, 0, 4/3, 0, 0; -1/12, 4/3, -5/2, 4/3, -1/12; 0, 0, 4/3, 0, 0; 0, 0, -1/12, 0, 0];
% % 4th DERIVATIVE;
% st=5; % stencil dimension
% mask= [0,0,1,0,0; 0,0,-4,0,0; 1,-4,12,-4,1; 0,0,-4,0,0; 0,0,1,0,0]/d0^2;
 
Mask=padarray(mask, [ny0 nx0]-st,'post');
M=circshift(Mask(:),ny0*nx0-((st-1)*0.5*ny0+(st-1)/2));clear st
D=zeros(ny0*nx0);
for ii=1:ny0*nx0
    D(ii,:)=circshift(M.',[0 ii-1]);
end
mu=1;
Dcs=[D ; mu*eye(nx0*ny0)];
%figure;surf(D);shading flat;colormap(flipud(gray))
clear mask Mask M
tic
cvx_begin quiet
    variable q_ftgv(n0) complex;
%    minimize norm(D*q_tv,1)
    minimize norm(Dcs*q_ftgv,1)
     subject to
        (norm(A*q_ftgv-ph)) <= rnl;
cvx_end
toc
 
%% FUSED TOTAL GENERALIZED VARIATION - FUSED Laplacian; 9 point stencil dimension 5 x 5
mask= [0, 0, 1, 0, 0; 0, 0, -16, 0, 0; 1, -16, 60, -16, 1; 0, 0, -16, 0, 0; 0, 0, 1, 0, 0];
st=5;
Mask=padarray(mask, [ny0 nx0]-st,'post');
M=circshift(Mask(:),ny0*nx0-((st-1)*0.5*ny0+(st-1)/2));clear st
D=zeros(ny0*nx0);
for ii=1:ny0*nx0
    D(ii,:)=circshift(M.',[0 ii-1]);
end
mu=1;
Dcs=[D ; mu*eye(nx0*ny0)];
%figure;surf(D);shading flat;colormap(flipud(gray))
clear mask Mask M
tic
cvx_begin quiet
    variable q_ftgv9(n0) complex;
%    minimize norm(D*q_tv,1)
    minimize norm(Dcs*q_ftgv9,1)
     subject to
        (norm(A*q_ftgv9-ph)) <= rnl;
cvx_end
toc