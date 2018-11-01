function [yf_wp,squaredEIGF0,sel_Eig2,Nh]=Supress_undesired_peaks(undesired_peak, h,gm,y,fs)

N=max(size(y));
desired_peak=[undesired_peak(end)+1:N];
TOP_peak=find(y==max(y));
Lcl = (1/(2*sqrt(pi)))*(gamma(gm+1)/gamma(gm+(3/2)));

[h, yf_wp0,Nh,psinnor,kappa,Ymin]= SCSA_1D(y,fs,h,gm);

squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;
squaredEIGF=(psinnor.^2)*kappa;

%% Method 1

% EIGF_th= max(max(squaredEIGF(undesired_peak,:)));
% squaredEIGF=squaredEIGF-EIGF_th;
% squaredEIGF_binary=squaredEIGF>0;
% selected_eigenfunctions = ~all(squaredEIGF_binary == 0);
% new_kappa=kappa.*diag(selected_eigenfunctions);
% yws =((h/Lcl)*sum((psinnor.^2)*new_kappa,2)).^(2/(1+2*gm));
% yws = yws';
% yws = yws +  Ymin;
% yws=yws-min(yws)+min(y);
% yf_wp11=yws;
% 
% squaredEIGF=(psinnor.^2)*kappa;
% sel_Eig=selected_eigenfunctions;

%% Method 2
new_kappa2=0*kappa;
sel_Eig2=diag(kappa)*0;
for k=1:Nh
  idx=find(squaredEIGF(1:TOP_peak,k)==max(squaredEIGF(1:TOP_peak,k)));
  
  if undesired_peak(end)<idx %| idx<undesired_peak(1) 
      new_kappa2(k,k)=kappa(k,k);
      sel_Eig2(k)=1;

  end
  
end

sel_Eig2;
yws2 =((h/Lcl)*sum((psinnor.^2)*new_kappa2,2)).^(2/(1+2*gm));
yws2 = yws2';
yws2 = yws2 +  Ymin;
yf_wp=yws2-min(yws2)+min(y);

