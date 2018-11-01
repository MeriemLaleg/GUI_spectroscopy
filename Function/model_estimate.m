function [sig_sum,C,D,Z_sig] = model_estimate(x,t,model_select)

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
            D(:,j) = 1.*exp((-d_l + 1i*2*pi*fre(j)).*t);
            count = count + 1;
        case 2
%             C(:,j) = exp((-damp(count)*t + 1i*2*pi*fre(j)).*t);
            d_g = (damp(count)).^2/(2*sqrt(log(2)));
            C(:,j) = amp(j)*exp((-d_g*t + 1i*2*pi*fre(j)).*t);
            D(:,j) = 1*exp((-d_g*t + 1i*2*pi*fre(j)).*t);
            count = count+1;
        case 3
            d_l = damp(count);
%             d_g = damp(count);
            d_g = (damp(count)).^2/(2*sqrt(log(2)));
            C(:,j) = amp(j)*exp((-d_l -d_g*t + 1i*2*pi*fre(j)).*t);      
            D(:,j) = 1*exp((-d_l -d_g*t + 1i*2*pi*fre(j)).*t);  
            count = count+1;
    end
end
sig_sum=sum(C,2);
Z_sig=sum(D,2);
