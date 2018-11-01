function y = tprod(A1,A2,x,value,m,n)
%
% y = tprod(A1,A2,x,value,m,n)
%
% This routine computes 
% T*x for value = 0 or T'*x for value = 1 
% where T is an m-by-n Toeplitz matrix. There are 
% no restrictions on the nonzero dimensions m and n. 
% For efficiency, due to the FFT, m+n-1 should be a power of 2. 
% A1 sets the stage  for T*x: A1 = toepinit(colT,rowT)
% A2 sets the stage for T'*x: A2 = toepinit(rowT,colT).
% Note: TPROD.M is used in TKSVD.M through TPRITZP.M
%
% Ricardo Fierro, CSUSM
% last revised: June 15, 2000.

 if value == 0
    % compute the m-by-1 vector y = T*x
    lena1 = length(A1);
    y = ifft(A1 .* fft([x; zeros(lena1-n,1)]));
    y = y(1:m,1);
 elseif value == 1
    % compute y = T'*x
    lena2 = length(A2);

    y = ifft(A2 .* fft([x; zeros(lena2-m,1)]));
    y = y(1:n,1);
 end

 % END OF TPROD.M
