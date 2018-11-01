function EEC = toepinit(colT,rowT);
%
%  EEC = toepinit(colT,rowT);
%
%  EEC contains eigenvalues for matrix-vector 
%  multiplication with a Toeplitz matrix.
%
%  Input:
%  colT, contains first column of Toeplitz matrix
%  rowT, contains first row of Toeplitz matrix
%
%  Output:
%  EEC, contains eigenvalues of extended circulant matrix, 
%  for use in toepmult.m to perform a Toeplitz 
%  matrix-vector multiplication.

% J. Nagy, SMU.
% modified by Rick Fierro, CSUSM, March 9, 1999.

% ensure colT is a *column* vector
[mcolT,ncolT] = size(colT);
if mcolT == length(colT)
  colT = colT(:,1);
else
   % use non-conjugate transpose
  colT  = colT(1,:).';
  mcolT = ncolT;
end

% ensure rowT is a *column* vector
[mrowT,nrowT] = size(rowT);
if mrowT == length(rowT)
  rowT = rowT(:,1);
else
   % use non-conjugate transpose
  rowT  = rowT(1,:).';
  mrowT = nrowT;
end

% find the next power of 2
p = nextpow2(mcolT + mrowT - 1);
n = 2^p;

% initialize
c = zeros(n,1);
% create the first column of the n-by-n circulant matrix
c(1:mcolT,1) = colT;
c(n-mrowT+2:n,1) = flipud(rowT(2:mrowT,1));
% compute the eigenvalues of the circulant matrix
EEC = fft(c);

% END OF TOEPINIT.M