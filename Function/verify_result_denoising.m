clear all;
clc;

addpath (strcat(pwd,'/Input_data/'));
addpath (strcat(pwd,'/Output_data/'));
load('dataset_2.mat');
complex_fid_met=complex_fid_unsuppressed;
load('dataset_2_output_HLSVD.mat');
hsvd_result=data_value.FID_TD;
NAA_Cr_ratio_HSVD=data_value.NAA_Cr;
Cho_Cr_ratio_HSVD=data_value.Cho_Cr;
snr_HSVD=data_value.snr;
load('dataset_2_output_SCSA.mat');
scsa_result=data_value.FID_TD;
Cho_Cr_ratio=data_value.Cho_Cr;
NAA_Cr_ratio_SCSA=data_value.NAA_Cr;
Cho_Cr_ratio_SCSA=data_value.Cho_Cr;
snr_SCSA=data_value.snr;
load('dataset_2_output.mat');
org_result=data_value.FID_TD;
Cho_Cr_ratio=data_value.Cho_Cr;
NAA_Cr_ratio_org=data_value.NAA_Cr;
Cho_Cr_ratio_org=data_value.Cho_Cr;
snr_org=data_value.snr;
zerofill=2;
close all;
for count=1:4
    X(:,count)=interp((real(fftshift(fft(hsvd_result(:,count)))))-(max(real(fftshift(fft(hsvd_result(:,count)))))-max(real(fftshift(fft(scsa_result(:,count)))))),zerofill);
    Y(:,count)=interp((real(fftshift(fft(scsa_result(:,count)))))-(max(real(fftshift(fft(scsa_result(:,count)))))-max(real(fftshift(fft(scsa_result(:,count)))))),zerofill);
    Z(:,count)=interp((real(fftshift(fft(org_result(:,count)))))-(max(real(fftshift(fft(org_result(:,count)))))-max(real(fftshift(fft(scsa_result(:,count)))))),zerofill);
%     figure();
%     plot(ppm,real(fftshift(fft(complex_fid_met(:,count)))));
%     hold on
%     plot(ppm,X);
%     plot(ppm,Y);
%     plot(ppm,Z);
%     axis([0, 10, -10, 50]);
%     xlabel('PPM');
%     ylabel('Amplitude');
%     legend('Original without water peak','HSVD','SCSA','Original with water peak');
%     set(gca,'XDir','reverse');
%     title('Spectra');
    fres= (-(Fs)/2 + ((Fs)/(size(X,1)*1))*(0:((size(X,1)*1-1))));
    details_para.fres=fres;
    ppm = (fres)*(1E6/Tf);
    ppm=ppm + 4.7 ;
    seg1=find(ppm>=0);
    seg1=seg1(1);
    seg2=find(ppm<=4.2);
    seg2=seg2(end);
end

figure;
subplot(1,4,1);
plot(ppm,Z(:,1),'linewidth',1.5);
hold on 
plot(ppm,X(:,1)+10000,'linewidth',1.5);
plot(ppm,Y(:,1)+30000,'linewidth',1.5);
axis([0, 6, -10000, 80000]);
xlabel('PPM');
ylabel('Amplitude');
set(gca,'XDir','reverse');
set(gca,'fontsize', 16,'fontweight','bold');
% legend('data', 'SVD', 'SCSA');
title('(a)');

subplot(1,4,2);
plot(ppm,Z(:,2),'linewidth',1.5);
hold on 
plot(ppm,X(:,2)+10000,'linewidth',1.5);
plot(ppm,Y(:,2)+30000,'linewidth',1.5);
axis([0, 6, -10000, 80000]);
xlabel('PPM');
% ylabel('Amplitude');
set(gca,'XDir','reverse');
set(gca,'fontsize', 16,'fontweight','bold');
% legend('data', 'SVD', 'SCSA');
title('(b)');

subplot(1,4,3);
plot(ppm,Z(:,3),'linewidth',1.5);
hold on 
plot(ppm,X(:,3)+10000,'linewidth',1.5);
plot(ppm,Y(:,3)+30000,'linewidth',1.5);
axis([0, 6, -10000, 80000]);
xlabel('PPM');
% ylabel('Amplitude');
set(gca,'XDir','reverse');
set(gca,'fontsize', 16,'fontweight','bold');
% legend('data', 'SVD', 'SCSA');
title('(c)');

subplot(1,4,4);
plot(ppm,Z(:,4),'linewidth',1.5);
hold on 
plot(ppm,X(:,4)+10000,'linewidth',1.5);
plot(ppm,Y(:,4)+30000,'linewidth',1.5);
axis([0, 6, -10000, 80000]);
xlabel('PPM');
% ylabel('Amplitude');
set(gca,'XDir','reverse');
set(gca,'fontsize', 16,'fontweight','bold');
legend('Original data', 'SVD', 'SCSA');
title('(d)');




