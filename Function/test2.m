clear all
clc

global details_para


% load('Example_MRS_data.mat') 
load_MRS_data

details_para.fla=0;

details_para.seg = [1 length(complex_fid_unsuppressed_FD)];
details_para.ppm_referenced=ppm;

Skip_W_Estim=0;     % Skip the water estimation if the baseline of the  metabolites is small 
gm=0.5;fs=1;        % SCSA parameters
Nh_WP=8;            % the minimum number of eigenfunction
WR_Ratio = 1 ;      % Water Reduction Ratio 
MaxIter= 100;       % The maximum iteration for water residual reduction loop
yf_SCSA=SCSA_MRS_Water_Suppression_GUI(Skip_W_Estim,ppm,freq,freq_vector, complex_fid_unsuppressed_FD,gm,fs,Nh_WP,WR_Ratio, MaxIter);

%% Location where there is almost noise

global shif show_plot
shif=0; close all
Denoising_coeff=100;
width_peaks=5;          %  The ratio w.r.t metabolite amplitude to set the threshold from which the peak can be selected
Th_peaks_ratio=2;        %  The with of the peak from its center(max values)


[yf_SCSA_Denoised, h_op, Nh]=SCSA_MRS_Denoising(500, ppm, yf_SCSA, gm , fs , Th_peaks_ratio, width_peaks);

cnt=1;shif=0; show_plot=1;close all
lgndDnzd{1}='Noisy MRS signal';
for h=h_op*[1:0.5:2]
    [yscsa,Nh]= Apply_SCSA_for_denoising(yf_SCSA, ppm,  h,gm,fs)
    shif=shif+0.1;pause(0.3)
   legend({'Noisy input spectrum ', 'Denoised Spectrum ','Residue'},'Location','northwest');
    cnt=cnt+1;lgndDnzd{cnt}=strcat('SCSA denosing using h= ',num2str(h),', N_h=',num2str(Nh));
end
legend(lgndDnzd);















% %% Location where there is almost noise
% 
% [No_pks1,pk_loc1,FWHM1,signal_power1(1,:)] = find_peaks_new(ifft(ifftshift(yf_SCSA)),freq,'');
% 
% m1= pk_loc1(1)-((length(yf_SCSA)-pk_loc1(1))./2); m2=pk_loc1(1)+ (length(yf_SCSA)-pk_loc1(1));
% Noise=[round(m1),round(m2)];
% 
% %% Location of metabolites 
% 
% [No_pks2,pk_loc2,FWHM2,signal_power2(1,:)] = find_peaks_new(ifft(ifftshift(yf_SCSA)),freq,'Please enter the number of metabolites in your MRS signal');
% 
% for i=1:No_pks2
%     Metabolite(i,:)=[round(pk_loc2(i)-(FWHM2(1,i)*length(yf_SCSA)./Fs)) round(pk_loc2(i)+(FWHM2(1,i)*length(yf_SCSA)./Fs))]; 
% end
% gm=0.5; fs=1;
% [yf_SCSA_Denoised, h_op, Nh]=SCSA_MRS_Denoising(ppm, yf_SCSA, Metabolite, Noise, gm , fs );
% clearvars -except ppm freq freq_vector  complex_fid_unsuppressed_FD yf_SCSA_Denoised  h_op  Nh

