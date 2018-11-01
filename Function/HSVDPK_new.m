%% Function to suppress water from time domain signal 

function [fit_sig,amp,Ind_sig,freq_temp,sit,damp] = HSVDPK_new(pc_fid,fres0,Fs,modelorder) 

N = length(pc_fid);
n = 1:N;
k=modelorder;  % No. of selecteed peaks (model order)
cur_pol = fres0/Fs;
y1=exp(-1i*2*pi*cur_pol*n');
y = pc_fid;%.*exp(-1i*2*pi*cur_pol*n');  % Signal used for decomposition 
% L = 513;
L = N/2 + 1; % length of decomposition window
M = N+1-L;
colhmat=y(N-L+1:N);
rowhmat=y(N-L+1:-1:1);
%flops(0);
[U,sv,V,info1,info2]=lansvdMNT(colhmat,rowhmat,k,'L');
%fl =flops;
% h_mat1 = hankel(y,y(L:end));  % Hankel matrix of the signal using the particualr window length
% h_r = h_mat1(1:L,:);  % conversion fo hankel matrix to a square matrix
% 
% [U,sv,V] = svd(h_r);   % SVD of the hankel matrix performed 
% [U_r,sv_r,V_r] = svd(h_r);

% U_top = U(2:end,1:k);   % top portion( all the rows except the 1st one and considering the 1st 5 columns) of the left singular vector matrix 
% U_bottom = U(1:end-1,1:k);  %bottom portion( all the rows except the last one and considering the 1st 5 columns) of the left singular vector matrix
U_top = U(2:L,1:k);   % top portion( all the rows except the 1st one and considering the 1st 5 columns) of the left singular vector matrix 
U_bottom = U(1:L-1,1:k);  %bottom portion( all the rows except the last one and considering the 1st 5 columns) of the left singular vector matrix
X_til = pinv(U_bottom)*U_top;  % Z vector obtained from multiplication of pseudoinverse of U bottom and U top
[E_Vect,tmp] = eig(X_til);    % eigen values and eigen vectors of the Z matrix
evs = diag(tmp);        % All the eigen values of Z matrix
ph=angle(evs);          
freq=ph*Fs./(2*pi);     % frequency components of the peaks 
% ind=find(freq>=-100 & freq<=100);
ind=freq;
no_ind=length(ind);
a=Fs/(2*pi);
for i=1:k
    freq_temp(i)=imag(log(evs(i))*a);
    damp(i)=-real(log(evs(i))*Fs);
end
for l=1:no_ind
    sit(:,l) = evs(l).^n;  % reconstructed signal from the peaks
end
% pc_r = real(pc_fid);
% pc_i = imag(pc_fid);
% sit1 = [real(sit);imag(sit)];
% amp = pinv(sit1)*[pc_r;pc_i];
amp = pinv(sit)*pc_fid;   % Amplitude estimates of the peaks 
% amp=ones(size(sit,2),1);
fit_sig = sit*amp;        % Final reconstructed signal
Ind_sig = sit*diag(amp);  % Each reconstructed peaks




