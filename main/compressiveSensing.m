function x = compressiveSensing(H,M,pii,Nnorm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cvx_begin quiet
cvx_precision high
variable p(M) complex;
minimize norm(p,1);
subject to
norm((H*p-pii),2) <= Nnorm;
cvx_end

x = p;
end

