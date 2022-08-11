function [C, Cpsi, i_qr] = qr_sampling(p,n,psi)
% Selection of samples using QR decompostion as in Manohar2018.
% IN:
%  p: number of selected positions
%  n: number of candidate positions
%  psi: sensing matrix (n x l) 
% OUT:
%  C: selection matrix 
%  Cpsi: sensing matrix
%  i_qr: index of selected positions

[Q, R, pivot] = qr(psi*psi','vector');
i_qr = pivot(1:p);

C = sparse(1:p, i_qr, ones(1,p), p, n);
C = full(C); 

Cpsi = C*psi;
end