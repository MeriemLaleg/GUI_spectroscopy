clear all;
clc;

addpath (strcat(pwd,'/Input_data/'));
addpath (strcat(pwd,'/Output_data/'));
load('dataset_1.mat');
complex_fid_met=complex_fid_unsuppressed;
load('dataset_1_output_HLSVD.mat');
hsvd_result=data_value.FID_TD;
NAA_Cr_ratio_HSVD=data_value.NAA_Cr;
Cho_Cr_ratio_HSVD=data_value.Cho_Cr;
load('dataset_1_output_SCSA.mat');
scsa_result=data_value.FID_TD;
Cho_Cr_ratio=data_value.Cho_Cr;
NAA_Cr_ratio_SCSA=data_value.NAA_Cr;
Cho_Cr_ratio_SCSA=data_value.Cho_Cr;
zerofill=2;
close all;
for count=1:40
    X(:,count)=interp((real(fftshift(fft(hsvd_result(:,count)))))-(max(real(fftshift(fft(hsvd_result(:,count)))))-max(real(fftshift(fft(scsa_result(:,count)))))),zerofill);
    Y(:,count)=interp((real(fftshift(fft(scsa_result(:,count)))))-(max(real(fftshift(fft(scsa_result(:,count)))))-max(real(fftshift(fft(scsa_result(:,count)))))),zerofill);
    Z(:,count)=interp((real(fftshift(fft(complex_fid_unsuppressed(:,count)))))-(min(real(fftshift(fft(complex_fid_unsuppressed(:,count)))))-min(real(fftshift(fft(scsa_result(:,count)))))),zerofill);
    met_signal(:,count)=(real(fftshift(fft(complex_fid_met(:,count)))));
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
    seg1=find(ppm>=4.3);
    seg1=seg1(1);
    seg2=find(ppm<=5.1);
    seg2=seg2(end);
    error_HSVD_seg(:,count)=(X(seg1:seg2,count));%(real(fftshift(fft(complex_fid_met(seg1:seg2,count))))-X(seg1:seg2));%./max((real(fftshift(fft(complex_fid_met(:,count))))));
    error_SCSA_seg(:,count)=(Y(seg1:seg2,count));%(real(fftshift(fft(complex_fid_met(seg1:seg2,count))))-Y(seg1:seg2));%./max((real(fftshift(fft(complex_fid_met(:,count))))));
    error_HSVD_seg_val(count)=var(error_HSVD_seg(:,count))-var(X(end-50:end,count));%.*100./norm(real(fftshift(fft(complex_fid_met(:,count)))));
    error_SCSA_seg_val(count)=var(error_SCSA_seg(:,count))-var(Y(end-50:end,count));%.*100./norm(real(fftshift(fft(complex_fid_met(:,count)))));
end
figure();
boxplot([abs(error_SCSA_seg_val)',abs(error_HSVD_seg_val)'],'Notch','on','Labels',{'SCSA','HLSVD-PRO'},'colors','rb');
% title('Residual error');
set(gca,'fontsize', 20,'fontweight','bold');
set(findobj(gca,'type','line'),'linew',2);
h = findobj(gcf,'tag','Outliers');
set(h(2),'MarkerEdgeColor','r');
set(h(1),'MarkerEdgeColor','b');
SCSA_invivo=abs(error_SCSA_seg_val);
HLSVD_invivo=abs(error_HSVD_seg_val);


figure();
subplot(1,2,1);
boxplot([NAA_Cr_ratio_SCSA',NAA_Cr_ratio_HSVD'],'Notch','on','Labels',{'SCSA','HLSVD-PRO'},'colors','rb');
title('NAA/Cr ratio');
set(gca,'fontsize', 20,'fontweight','bold');
set(findobj(gca,'type','line'),'linew',2);
h = findobj(gcf,'tag','Outliers');
set(h(2),'MarkerEdgeColor','r');
set(h(1),'MarkerEdgeColor','b');
SCSA_NAA_Cr=NAA_Cr_ratio_SCSA;
HLSVD_NAA_Cr=NAA_Cr_ratio_HSVD;


subplot(1,2,2);
boxplot([Cho_Cr_ratio_SCSA',Cho_Cr_ratio_HSVD'],'Notch','on','Labels',{'SCSA','HLSVD-PRO'},'colors','rb');
title('Cho/Cr ratio');
set(gca,'fontsize', 20,'fontweight','bold');
set(findobj(gca,'type','line'),'linew',2);
h = findobj(gcf,'tag','Outliers');
set(h(2),'MarkerEdgeColor','r');
set(h(1),'MarkerEdgeColor','b');
SCSA_Cho_Cr=Cho_Cr_ratio_SCSA;
HLSVD_Cho_Cr=Cho_Cr_ratio_HSVD;




figure;
subplot(1,3,1);
plot(ppm,Z(:,3),'linewidth',1.5);
% axis([0, 6, min(Z(:,1)), max(Z(:,1))]);
axis([0, 6, -2, 5]);
xlabel('PPM');
ylabel('Amplitude');
set(gca,'XDir','reverse');
set(gca,'fontsize', 20,'fontweight','bold');
title('(a)');
set(gca,'fontsize', 20,'fontweight','bold');

subplot(1,3,2);
plot(ppm,Y(:,3),'linewidth',1.5);

axis([0, 6, -2, 5]);
xlabel('PPM');
% ylabel('Amplitude');
set(gca,'XDir','reverse');
set(gca,'fontsize', 20,'fontweight','bold');
title('(b)');
set(gca,'fontsize', 20,'fontweight','bold');

subplot(1,3,3);
plot(ppm,X(:,3),'linewidth',1.5);

axis([0, 6, -2, 5]);
xlabel('PPM');
% ylabel('Amplitude');
set(gca,'XDir','reverse');
set(gca,'fontsize', 20,'fontweight','bold');
title('(c)');
set(gca,'fontsize', 20,'fontweight','bold');




