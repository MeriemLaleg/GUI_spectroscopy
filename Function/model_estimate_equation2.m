function [err] = model_estimate_equation2(x,t,y,model_select)

amp=x(:,1);
damp=x(:,2);
fre=x(:,3);
fun_type=ones(1,length(fre))*model_select;
K = length(fre);
C = zeros(length(t),K);
count = 1;
for j = 1:K
    switch fun_type(j)
        case 1
            d_l = damp(count);
            C(:,j) = amp(j).*exp((-d_l + 1i*2*pi*fre(j)).*t);
            count = count + 1;
        case 2
%             C(:,j) = exp((-damp(count)*t + 1i*2*pi*fre(j)).*t);
            d_g =(damp(count)).^2/(2*sqrt(log(2)));
            C(:,j) = amp(j)*exp((-d_g*t + 1i*2*pi*fre(j)).*t);
            count = count+1;
        case 3
            d_l = damp(count);
%             d_g = damp(count);
            d_g = (damp(count)).^2/(2*sqrt(log(2)));
            C(:,j) = amp(j)*exp((-d_l -d_g*t + 1i*2*pi*fre(j)).*t);  
            count = count+1;
    end
end
sig_sum=sum(C,2);
% Z=abs(fftshift(fft(y)));
% err1=Z-abs(fftshift(fft(sig_sum)));
err1=y-sig_sum;
err =[real(err1);imag(err1)];
% err =err1;% [real(err1);imag(err1)];
