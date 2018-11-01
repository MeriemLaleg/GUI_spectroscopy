function [sigma,bnd] = bdsqr(alpha,beta)
  
% [SIGMA,BND] = BDSQR(ALPHA,BETA)
%
  

k = length(alpha);
if min(size(alpha)') ~= 1  | min(size(beta)') ~= 1
  error('alpha and beta must be vectors')
elseif length(beta) ~= k
  error('alpha and beta must have the same lenght')
end    
B = spdiags([alpha(:),beta(:)],[0,-1],k+1,k);
[U,S,V] = svd(full(B),0);
sigma = diag(S);
bnd = U(end,1:k)';


