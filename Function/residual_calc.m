function [res] = residual_calc(FD_sig_wbc, win_size, hw,freq,FWHM,pk_loc)
global details_para;
global info_struct;


ini_f = freq;
total_size = length(FD_sig_wbc);
n_seg = total_size/win_size;

% Finding Noise standard deviation
for i=1:n_seg
    temp = FD_sig_wbc((i-1)*win_size+1:i*win_size);
    mu(i) = mean(temp);
    sigma(i) = (1/(win_size-1))*sum((temp - mean(temp)).^2);
end
sigma_noise = abs(sqrt(min(sigma)));

% Marking residual points
res = FD_sig_wbc(1:hw-1);
for i=hw:total_size-hw
    temp = FD_sig_wbc(i+1-hw:i+hw);
    difference = max(temp) - min(temp);
    if abs(difference) <= 15*sigma_noise % vary this to change points in baseline and signal-decrease to increase points in signal, decrease points in baseline; increase to reduce points in signal
        res = [res; FD_sig_wbc(i)];
    else
        res = [res; NaN];
    end   
end
res = [res; FD_sig_wbc(end-hw+1:end)];

% Removing selected peak points from residual
fres_dif = details_para.fres(2) - details_para.fres(1);
for i=1:length(ini_f)
    hw = round(FWHM(i)*1/fres_dif);
%     pk_loc = floor((ini_f(i)+ details_para.Fs/2)*details_para.N/details_para.Fs);
    res(pk_loc-hw:pk_loc+hw,1) = NaN;
end

% interpolating over signal points
discont = find(isnan(res));
next = discont(1);
for i=1:length(discont)-1
    temp = discont(i+1)-discont(i);
    if (temp~=1)
        res(next:discont(i))= linspace(FD_sig_wbc(next),FD_sig_wbc(discont(i)),discont(i)-next+1)';
        next = discont(i+1);
    end
end
res(next:discont(end))= linspace(FD_sig_wbc(next),FD_sig_wbc(discont(end)),discont(end)-next+1)';

% x_scale = details_para.ppm_referenced;
% min_p = min(FD_sig_wbc);
% max_p = max(FD_sig_wbc);
% plot(x_scale, FD_sig_wbc)
% hold on
% plot(x_scale, res,'r-')
% axis(gca,[x_scale(1), x_scale(end), min_p - 0.1*max_p, max_p + 0.1*max_p])
% set(gca,'Xdir', 'reverse')

