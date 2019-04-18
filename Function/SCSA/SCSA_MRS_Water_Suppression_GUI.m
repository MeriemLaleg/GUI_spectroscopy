%% Water Suppression using SCSA Method V2: 
% The Semi-Classical Signal Analysis (SCSA)method is based for MRS signal 
% water suppression  consists of two steps: 
% the first step: is to estimate the water peak formulated as 
% an optimization problem. It reconstructs the water peak of a minimum number 
% of eigenfunction using an  optimal choice of the value of $h$. 
% The second step: removes some eigenfunctions that contribute to 
% the metabolites to recover only the water peak.
%
% Important: the input data are stored in "./Input_data/*.mat' files  containning:
% [img : reference WS-MRS signal] , [original_file : name of source file] 
% [noisy_img :  MRS signal], [noisy_file : name of source file] 
%% ######################  PARAMETERS ######################################
% Input 
% ppm         : The MRS spectrum fequancies in ppm
% freq        : The MRS spectrum fequancies in Hz
% freq_vector : The MRS metabolites  fequancy range in ppm
% complex_fid_unsuppressed_FD: Unsuppressed  complex MRS FID signal
% gm, fs  : SCSA parameter: gm=0.5, fs=1
% Nh_WP : Number of eigenfunction to be used to estimate the water peak

% Global variables:
% Loop_Water_reduction: It reduces the water residual. The values range is   from 1 to 4
% Nh_WP_reduction: It make the water residual more flat . The values range is   from 1 to 3

% Output
% yf_SCSA: Water suppressed  complex MRS FID signal

%% ###########################################################################
%  Author:
%  Abderrazak Chahid (abderrazak.chahid@gmail.com)
%  Adviser:
%  Taous-Meriem Laleg (taousmeriem.laleg@kaust.edu.sa)
% Done: August,  2018
% King Abdullah University of Sciences and Technology (KAUST)

function yf_SCSA=SCSA_MRS_Water_Suppression_GUI(ppm,freq,freq_vector, complex_fid_unsuppressed_FD,gm,fs,Nh_WP)
global ind_Metabolite ind_Water  Plot_fig show_plot Loop_Water_reduction; Plot_fig=0;%show_plot=1;

%% Generate signals
fprintf('\n_____________________________________________________________')
fprintf('\n Water suppression using SCSA  (By Abderrazak - KAUST- 2018)')
fprintf('\n_____________________________________________________________\n')

tic
yf=real(complex_fid_unsuppressed_FD)';
N=max(size(yf));f=ppm;


freq_desired=freq_vector;
ind_Metabolite= find( freq_desired(1)<= freq & freq<=freq_desired(2));
ind_Water= find( -50<= freq & freq<=50) ;%  

%% Setup the  Wait bar
wait_bar = waitbar(0.12,'Water peak estimation','Name','MRS water suppression using SCSA');


%% Step1: Water Peak Estimation from  the input MRS  singal 

[yf_SCSA1,yf_wp0,yf_wp]=Water_peak_estimation(f,yf,fs,Nh_WP,gm,ind_Metabolite);

waitbar(0.50,wait_bar, 'Water peak refinement using eigenfunctions selection ')                   % Wait bar progression 
% close all ;  figure; plot(f,yf_SCSA1);hold off; legend('yf','Water peak phase 1','water peak refined',' suppressed water peak')

d=1;
    
%% Step 2:  Water Residual suppression by  reducing  its  amplitude by SCSA 
fprintf('\n--> Water residual Reduction.')
waitbar(0.56,wait_bar, 'Water Residual Reduction')                   % Wait bar progression

for jj=1:Loop_Water_reduction
    yf_SCSA1=Remove_base_no_metabolites(yf_SCSA1,ind_Water, fs,gm);
end
%close all ;figure; plot(f,yf);hold on;   plot(f,yf_SCSA1);%hold on;   plot(f,yf-yf_SCSA1);

d=1;  
% close all ;figure; plot(f,yf_SCSA1);hold on;   plot(f,yf_SCSA1);%hold on;   plot(f,yf-yf_SCSA1);
yf_SCSA=yf_SCSA1;

d=1;   
fprintf('\n--> Water suppression using SCSA is Done!!.\n\n')
waitbar(1,wait_bar)                   % Wait bar progression

         
%% plot the results 
if show_plot==1
    yf=yf-min(yf); ymin=min(yf) ;
    yf_SCSA=yf_SCSA -min(yf_SCSA) + max(min(yf(1:20))); 
    ymax=1.1*max([max(yf(ind_Metabolite)), max(yf_SCSA(ind_Metabolite))]);
    figure;
    plot(ppm,yf,'b','LineWidth',2);
    hold on
    plot(ppm, yf_SCSA ,'r','LineWidth',2)

    legend( 'Input MRS signal with water peak','Output MRS signal with suppressed water peak  using SCSA');
    xlabel('ppm')
    ylabel('Intensity')
%     set(gca,'YTickLabel',[])
    ylim([ymin ymax])
    title([ ' SCSA MRS-water suppression '])%with : h = ' num2str(hwp) , '   Nh = ' num2str(Nh_O) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
    set(gcf,'color','w') 
    set(gca,'Xdir','reverse');
    set(gca,'fontsize',16)
    % text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
    box 
end


close(wait_bar)

function yp=padd_signal(y,padd)
 yp=[y(1)*ones(1,padd) y y(end)*ones(1,padd)];
   
function y=unpadd_signal(yp,padd)
 y=yp(padd+1:end-padd);
   
function yf_SCSA1=Remove_base_no_metabolites(yf_SCSA1,area,fs,gm)
 global Nh_WP_reduction  
 
padd=floor(max(size(area))/5); 

if padd > 10; padd=10; end

 N=max(size(yf_SCSA1));
 yf_No_Mtab=yf_SCSA1(area);  
 yff=padd_signal(yf_No_Mtab,padd);
 
% Nh_WP_reduction= 2;% floor(N/100); % = 10; the reconstruct the global behavior     
[h3, yff_base,Nh_O]= SCSA_1D_Nh(yff,fs,Nh_WP_reduction,gm);
yf_22=unpadd_signal(yff_base,padd);
yf_22=yf_22-yf_22(1)+yff(1);
% yf_SCSA11(W0-5:N)=yf_22;

Res=abs(yf_No_Mtab-yf_22);
% close all ;figure; plot(yf_No_Mtab);hold on;   plot(yf_22);

% Idx0=find(Res==min(Res))
yf_SCSA11=yf_SCSA1;
yf_SCSA11(area)=Res+yf_SCSA1(area(1));
%  close all ;figure; plot(yf_SCSA1);hold on;   plot(yf_SCSA11);hold on;   plot(yf_SCSA1-yf_SCSA11);
   
d=1;

yf_SCSA1=yf_SCSA11;

function  [yf_SCSA1,yf_wp0,yf_wp]=Water_peak_estimation(f,yf,fs,Nh_WP,gm,ind_Metabolite)
        
%% Second methods eigenfundtion seletion
fprintf('\n--> Water peak estimation using eigenfunction selection')            

[hwp, yf_scsa,Nh_O]= SCSA_1D_Nh(yf,fs,Nh_WP,gm);


%%  Remove the small peak  caused by the metabolites

[yf_wp0,yf_wp,yf_wp2,squaredEIGF0,Sel_Eigf]=Supress_undesired_peaks(ind_Metabolite, hwp,gm,yf,fs);

yf_wp0=yf_wp0-min(yf_wp0)+min(yf);
yf_wp=yf_wp-min(yf_wp)+min(yf);

yf_SCSA=yf-yf_wp;yf_SCSA11=yf_SCSA-min(yf_SCSA(ind_Metabolite));
yf_SCSA1=yf_SCSA11;

% close all ;
% figure; plot(f,yf);hold on ;plot(f,yf_wp0);hold on ;plot(f,yf_wp);hold on;   plot(f,yf_SCSA1);hold off; legend('yf','Water peak phase 1','water peak refined',' suppressed water peak')
%  d=1;           

    
function [y_h, yf_SCSA,yf_SCSA2,squaredEIGF0,sel_Eig]=Supress_undesired_peaks(undesired_peak, h,gm,y,fs)

N=max(size(y));

TOP_peak=find(y==max(y));
Lcl = (1/(2*sqrt(pi)))*(gamma(gm+1)/gamma(gm+(3/2)));

[h, y_h,Nh,psinnor,kappa,Ymin]= SCSA_1D(y,fs,h,gm);

squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;
squaredEIGF=(psinnor.^2)*kappa;
EIGF_th= max(max(squaredEIGF(undesired_peak,:)));
squaredEIGF=squaredEIGF-EIGF_th;
squaredEIGF_binary=squaredEIGF>0;
selected_eigenfunctions = ~all(squaredEIGF_binary == 0);
new_kappa=kappa.*diag(selected_eigenfunctions);
yws =((h/Lcl)*sum((psinnor.^2)*new_kappa,2)).^(2/(1+2*gm));
yws = yws';
yws = yws +  Ymin;
yws=yws-min(yws)+min(y);
yf_SCSA=yws;
squaredEIGF=(psinnor.^2)*kappa;

sel_Eig=selected_eigenfunctions;

%% Method 2
new_kappa2=0*new_kappa;
selected_2=selected_eigenfunctions*0;
for k=1:Nh
  idx=find(squaredEIGF(1:TOP_peak,k)==max(squaredEIGF(1:TOP_peak,k)));
  
  if undesired_peak(end)<idx %| idx<undesired_peak(1) 
      new_kappa2(k,k)=kappa(k,k);
      selected_2(k)=1;

  end
  
end

selected_2;
yws2 =((h/Lcl)*sum((psinnor.^2)*new_kappa2,2)).^(2/(1+2*gm));
yws2 = yws2';
yws2 = yws2 +  Ymin;
yws2=yws2-min(yws2)+min(y);
yf_SCSA2=yws2;
d=1;
% figure;plot(yf_SCSA);hold on ;plot(yf_SCSA2);
% yf_SCSA=yf_SCSA-yf_SCSA(1)+y(1);


