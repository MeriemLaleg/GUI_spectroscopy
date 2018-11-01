
function  [yf_ws,yf_wp0,yf_wp]=Water_peak_estimation(f,yf,fs,Nh_WP,gm,ind_Metabolite)
        
%% Second methods eigenfundtion seletion
fprintf('\n--> Water peak estimation using eigenfunction selection')            

[hwp, yf_wp0,Nh_O]= SCSA_1D_Nh(yf,fs,Nh_WP,gm);

%%  Remove the small peak  caused by the metabolites

[yf_wp,squaredEIGF0,Sel_Eigf]=Supress_undesired_peaks(ind_Metabolite, hwp,gm,yf,fs);

yf_ws=yf-yf_wp;
close all ;figure; plot(f,yf);hold on ;plot(f,yf_wp0);hold on ;   plot(f,yf_wp);hold on;   plot(f,yf_ws);hold off;  legend('yf','Water peak phase 1','water peak refined',' suppressed water peak')

d=1;           

