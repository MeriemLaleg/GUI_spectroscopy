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
% Skip_W_Estim: set to 1 for invitro data
% ppm         : The MRS spectrum fequancies in ppm
% freq        : The MRS spectrum fequancies in Hz
% freq_vector : The MRS metabolites  fequancy range in ppm
% complex_fid_unsuppressed_FD: Unsuppressed  complex MRS FID signal
% gm, fs  : SCSA parameter: gm=0.5, fs=1
% Nh_WP   : SCSA parameter: gm=0.5, fs=1
% WR_Ratio: The  desired ratio betwwn the water residue and the metabolites  
% MaxIter:  The maximuim allowed iteration to achieve the desired WR_Ratio   

% Output
% yf_SCSA: Water suppressed  complex MRS FID signal

%% ###########################################################################
%  Author:
%  Abderrazak Chahid (abderrazak.chahid@gmail.com)
%  Adviser:
%  Taous-Meriem Laleg (taousmeriem.laleg@kaust.edu.sa)
% Done: May,  2018
% King Abdullah University of Sciences and Technology (KAUST)

function yf_SCSA=SCSA_MRS_Water_Suppression_GUI(ppm,freq,freq_vector, complex_fid_unsuppressed_FD,gm,fs,Nh_WP)
global ind_Metabolite ind_Water  Plot_fig show_plot ; Plot_fig=0;%show_plot=1;

hwp=-1;Nh_O=-1;
%% Generate signals
fprintf('\n_____________________________________________________________')
fprintf('\n Water suppression using SCSA  (By Abderrazak - KAUST- 2018)')
fprintf('\n_____________________________________________________________\n\n')

tic
yf=real(complex_fid_unsuppressed_FD)';
N=max(size(yf));f=ppm;


 freq_desired=freq_vector;
ind_Metabolite= find( freq_desired(1)<= freq & freq<=freq_desired(2));
ind_Water= find( -50<= freq & freq<=50) ;%  
peaks_pos=sort([ind_Metabolite(1) ind_Metabolite(end) ind_Water(1) ind_Water(end)]);


% peaks_pos_ref=peaks_pos;Counting=0;
% name_peaks={'Lac-NAA','Water'};
% peaks_time=ppm(peaks_pos);N_peaks=max(size(peaks_pos))/2 -1;
ppm_Max=ppm(peaks_pos(end-2)) ;ppm_zoom=max(find(f<=ppm_Max));ppm_zoom2=floor(peaks_pos(end-2)*1.3);

%% #Evaluation of SNR of the noisy image/signal
fprintf('\n--> Loading Data' )


%% Setup the  Wait bar
wait_bar =waitbar(0.12,'Water peak estimation','Name','MRS water suppression using SCSA');


%% Step1: Water Peak Estimation from  the input MRS  singal 
[yf_SCSA1,yf_wp0,yf_wp]=Water_peak_estimation(f,yf,fs,Nh_WP,gm,ind_Metabolite);

%close all ;figure; plot(yf);hold on ;plot(yf_SCSA1);hold on ;  plot(yf_wp0);plot(yf_wp);hold on
d=1;

%% Step 2:  Water Residual suppression by  reducing  its  amplitude by SCSA 
fprintf('\n--> Water residual Reduction.')

%% find the start position of the  Water peak
MW0=max(find(f<4)); MW1=min(find(f>5)); W1=MW1;

%% Find W0 and W1
dif_yf=diff(yf(MW0:end));W0=min(find(dif_yf>0))+MW0;
if W0<MW0,   W0=MW0; end;% if W1>MW1,   W1=MW1; end;
Left_W=1:W0;Right_W=W1+1:N; waterPeak=W0+1:W1;
yf_SCSAWL=yf_SCSA1(Left_W); yf_SCSAW=yf_SCSA1(waterPeak);yf_SCSAWR=yf_SCSA1(Right_W);

%close all ;figure; plot(yf_SCSAWL);hold on ;plot(yf_SCSAW);hold on ;  plot(yf_SCSAWR)
d=1;

%% Remove the water residue  
waitbar(0.56,wait_bar, 'Water Residual Reduction')                   %%   wait bar progression
yf_SCSA1=Remove_base_no_metabolites(yf_SCSA1,waterPeak, fs,gm);

         
%% Remove the non metabolite nase line        
area_no_Metb=W0-5:N;
waitbar(0.75,wait_bar, 'Water Residual Reduction')                   %%   wait bar progression

yf_SCSA=Remove_base_no_metabolites(yf_SCSA1,area_no_Metb, fs,gm);
      
d=1;   
fprintf('\n--> Done.')
%waitbar(1,wait_bar)                   %%   wait bar progression
         
%% plot the results 
if show_plot==1
    yf=yf-min(yf); ymin=min(yf) ;
    yf_SCSA=yf_SCSA -min(yf_SCSA) + max(min(yf(1:20))); 
    ymax=1.1*max([max(yf(Left_W)), max(yf_SCSA(Left_W))]);
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
    
padd=floor(max(size(area))/5); 
if padd > 10; padd=10; end


 N=max(size(yf_SCSA1));
 yf_No_Mtab=yf_SCSA1(area);  
 yff=padd_signal(yf_No_Mtab,padd);
 
Nh_No_Mtb= floor(N/50); % = 10; the reconstruct the global behavior     
[h3, yff_base,Nh_O]= SCSA_1D_Nh(yff,fs,Nh_No_Mtb,gm);
yf_22=unpadd_signal(yff_base,padd);
yf_22=yf_22-yf_22(1)+yff(1);
% yf_SCSA11(W0-5:N)=yf_22;

Res=abs(yf_No_Mtab-yf_22);
% %close all ;figure; plot(yf_No_Mtab);hold on;   plot(yf_22);

% Idx0=find(Res==min(Res))
yf_SCSA11=yf_SCSA1;
yf_SCSA11(area)=Res+yf_SCSA1(area(1));
%  %close all ;figure; plot(yf_SCSA1);hold on;   plot(yf_SCSA11);hold on;   plot(yf_SCSA1-yf_SCSA11);
   
d=1;

yf_SCSA1=yf_SCSA11;

function  [yf_SCSA1,yf_wp0,yf_wp]=Water_peak_estimation(f,yf,fs,Nh_WP,gm,ind_Metabolite)
        
%% Second methods eigenfundtion seletion
fprintf('\n--> Water peak estimation using eigenfunction selection')            

[hwp, yf_scsa,Nh_O]= SCSA_1D_Nh(yf,fs,Nh_WP,gm);


%%  Remove the small peak  caused by the metabolites

[yf_wp0,yf_wp,yf_wp2,squaredEIGF0,Sel_Eigf]=Supress_undesired_peaks(ind_Metabolite, hwp,gm,f,yf,fs);

yf_wp0=yf_wp0-min(yf_wp0)+min(yf);
yf_wp=yf_wp-min(yf_wp)+min(yf);

yf_SCSA=yf-yf_wp;yf_SCSA11=yf_SCSA-min(yf_SCSA(ind_Metabolite));
yf_SCSA1=yf_SCSA11;

% %close all ;
% figure; plot(f,yf);hold on ;plot(f,yf_wp0);hold on ;plot(f,yf_wp);hold on;   plot(f,yf_SCSA1);hold off; legend('yf','Water peak phase 1','water peak refined',' suppressed water peak')
%  d=1;           

function [Z,h]=Water_Res_reduction( WR_Ratio, MaxIter, AmpW0, y,fs,Nh0,tol,gm)
global  Plot_fig

y_inR=y;
fprintf('\n   --> Water reduction using %d Eigenfunction (+/- %d). Please wait.. \n',Nh0,tol)
Nh=1;
k=0;
if Nh0 < 9, Nh0=9;end
AmpW=10*AmpW0;AmpW_old=AmpW;


while AmpW>AmpW0/WR_Ratio && k<MaxIter %(Nh~=0) ||
    y_in=y_inR;k=k+1;
%     [h0, y_base,Nh]= SCSA_1D_Nh_Tolereance(y_in,fs,Nh0,tol,gm); %
    [h0, y_base,Nh]= SCSA_1D_Nh(y_in,fs,Nh0,gm);

    h(k)=h0;

    y_inR=y_in-y_base;
    y_inR=y_inR-min(y_inR)+min(y_in);
    if Plot_fig==1
    figure; plot(y);hold on;
            plot(y_in);hold on;
            plot(y_inR);hold on;
    end
    
    AmpW=peak2peak(y_inR);
    if abs(AmpW_old-AmpW)<0.5 
        
        if Nh0<max(size(y))/2
            Nh0=Nh0+2;  
        end
    end %&& Nh0 >7
    AmpW_old=AmpW;
    
    d=1;
   fprintf('.')
end

AmpWin=peak2peak(y_in);
% AmpW0
% AmpW

if abs(AmpW0-AmpWin)>abs(AmpW0-AmpW)
    Z=y_inR;
else
    Z=y_in;
end

% AmpZ=peak2peak(Z);
% if AmpZ >AmpW0
%     Z=(AmpW0/AmpZ)*Z;
% end


function [Z,h]=Extact_baseline_SCSA_Nh(loop, y,fs,Nh0,tol,gm)
global  Plot_fig

y_inR=y;
fprintf('\n   --> Optimization: Reconstruction of signal using %d Eigenfunction (+/- %d) in %d iterations. Please wait.. ',Nh0,tol,loop)

for k=1:loop
    y_in=y_inR;
%     [h0, y_base]= SCSA_1D_Nh_Tolereance(y_in,fs,Nh0,tol,gm); 
    [h0, y_base]= SCSA_1D_Nh(y_in,fs,Nh0,gm);

    
    h(k)=h0;
    y_inR=y_in-y_base;
    y_inR=y_inR-min(y_inR)+min(y_in);
%     if Plot_fig==1
%     figure; plot(y_in);hold on;
%             plot(y_inR);hold on;
%     end
    d=1;
end
Z=y_inR;


function [h, yscsa,Nh,psinnor,kappa,Ymin]= SCSA_1D_Nh_Tolereance(y,fs,Nh0,tol,gm)

global  Plot_fig


[h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D_Nh(y,fs,Nh0,gm)


h=max(y);N=max(size(y));
[h, yscsa,Nh,psinnor,kappa,Ymin]= SCSA_1D(y,fs,h,gm);

Lcl = (1/(2*sqrt(pi)))*(gamma(gm+1)/gamma(gm+(3/2)));
squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;

%% Searching for h corresponding to Nh
% fprintf('\n   --> Optimization: Reconstruction of signal using %d Eigenfunction. Please wait.. ',Nh0)
eta=2*Nh0/10;
% if Plot_fig==1
%     %close all
%     lgnd{1}='Input';
%     figure(1);subplot(211);plot(y);hold on
%     legend(lgnd)
% end

kk=1;
Delta_Nh=0.1*h;
while abs(Nh-Nh0)>=tol && kk<100
    Nh_old=Nh;kk=kk+1;
    h=abs(h+Delta_Nh);
    h_histo(kk)=h;
    [h, yscsa,Nh,psinnor,kappa,Ymin]= SCSA_1D(y,fs,h,gm);
    h_old=h;
    Delta_Nh=Nh-Nh0;
   
    if   Plot_fig==1
%         lgnd{kk}=strcat('Nh=',num2str(Nh));
%         figure(1);subplot(211);plot(yscsa);hold on
%         legend(lgnd)
%         subplot(212);plot(h_histo); 

        d=1;
    end
end

 d=1;








function [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D_Nh(y,fs,Nh0,gm)
global  Plot_fig
% h=max(y);
N=max(size(y));
h=max(y)-min(y);
[h, yscsa,Nh,psinnor,kappa,Ymin]= SCSA_1D(y,fs,h,gm);

eta=5;
if Plot_fig==1
%     lgnd{1}='Input';
%     figure(1);plot(y);hold on
%     legend(lgnd)
end
kk=1;
Delta_Nh=0.1*h;


while Nh/Nh0~=1
    Nh_old=Nh;
    kk=kk+1;
    if Nh==0
        h=h/2;
    else
    h=h*Nh/Nh0;
    end
    [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm);
    h_old=h;
    Delta_Nh=Nh-Nh0;
%    if Plot_fig==1
%         h
%         Nh
%         Delta_Nh
%         lgnd{kk}=strcat('Nh=',num2str(Nh));
%         figure(1);plot(yscsa);hold on
%         legend(lgnd)
%    end
    
    d=1;
end



% while abs(Nh-Nh0)~=0 && kk<100
%     Nh_old=Nh;kk=kk+1;
%     h=h+eta*Delta_Nh;
%     [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm);
%     h_old=h;
%     Delta_Nh=Nh-Nh0;
%    if Plot_fig==1
%         h
%         Nh
%         Delta_Nh
%         lgnd{kk}=strcat('Nh=',num2str(Nh));
%         figure(1);plot(yscsa);hold on
%         legend(lgnd)
%    end
    
%     d=1;
% end

    d=1;

 




function [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm)

Lcl = (1/(2*sqrt(pi)))*(gamma(gm+1)/gamma(gm+(3/2)));
N=max(size(y));
%% remove the negative part
Ymin=min(y);

 y_scsa = y -Ymin;
%% Build Delta metrix for the SC_hSA
feh = 2*pi/N;
D=delta(N,fs,feh);

%% start the SC_hSA
Y = diag(y_scsa);
SC_h = -h*h*D-Y; % The Schrodinger operaor

% = = = = = = Begin : The eigenvalues and eigenfunctions
[psi,lamda] = eig(SC_h); % All eigenvalues and associated eigenfunction of the schrodinger operator
% Choosing  eigenvalues
All_lamda = diag(lamda);
ind = find(All_lamda<0);


%  negative eigenvalues
Neg_lamda = All_lamda(ind);
kappa = diag((abs(Neg_lamda)).^gm); 
Nh = size(kappa,1); %%#ok<NASGU> % number of negative eigenvalues



if Nh~=0
    
% Associated eigenfunction and normalization
psin = psi(:,ind(:,1)); % The associated eigenfunction of the negarive eigenvalues
I = simp(psin.^2,fs); % Normalization of the eigenfunction 
psinnor = psin./sqrt(I);  % The L^2 normalized eigenfunction 

squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;

%yscsa =4*h*sum((psinnor.^2)*kappa,2); % The 1D SC_hSA formula
yscsa1 =((h/Lcl)*sum((psinnor.^2)*kappa,2)).^(2/(1+2*gm));
else
    
  psinnor = 0*y;  % The L^2 normalized eigenfunction 
  kappa=0;
  squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;

  yscsa1=0*y;
  yscsa1=yscsa1-10*abs(max(y));
%   disp('There are no negative eigenvalues. Please change the SCSA parameters: h, gm ')
  fprintf('\0')

end


if size(y_scsa) ~= size(yscsa1)
yscsa1 = yscsa1';
end
 
 %% add the removed negative part
 yscsa = yscsa1 + Ymin;






    %**********************************************************************
    %*********              Numerical integration                 *********
    %**********************************************************************

    % Author: Taous Meriem Laleg

    function y = simp(f,dx);
    %  This function returns the numerical integration of a function f
    %  using the Simpson method

    n=length(f);
    I(1)=1/3*f(1)*dx;
    I(2)=1/3*(f(1)+f(2))*dx;

    for i=3:n
        if(mod(i,2)==0)
            I(i)=I(i-1)+(1/3*f(i)+1/3*f(i-1))*dx;
        else
            I(i)=I(i-1)+(1/3*f(i)+f(i-1))*dx;
        end
    end
    y=I(n);
    

    %**********************************************************************
    %*********             Delata Metrix discretization           *********
    %**********************************************************************
    
    
%Author: Zineb Kaisserli

function [Dx]=delta(n,fex,feh)
    ex = kron([(n-1):-1:1],ones(n,1));
    if mod(n,2)==0
        dx = -pi^2/(3*feh^2)-(1/6)*ones(n,1);
        test_bx = -(-1).^ex*(0.5)./(sin(ex*feh*0.5).^2);
        test_tx =  -(-1).^(-ex)*(0.5)./(sin((-ex)*feh*0.5).^2);
    else
        dx = -pi^2/(3*feh^2)-(1/12)*ones(n,1);
        test_bx = -0.5*((-1).^ex).*cot(ex*feh*0.5)./(sin(ex*feh*0.5));
        test_tx = -0.5*((-1).^(-ex)).*cot((-ex)*feh*0.5)./(sin((-ex)*feh*0.5));
    end
    Ex = full(spdiags([test_bx dx test_tx],[-(n-1):0 (n-1):-1:1],n,n));
    
    Dx=(feh/fex)^2*Ex;

    
    
function [y_h, yf_SCSA,yf_SCSA2,squaredEIGF0,sel_Eig]=Supress_undesired_peaks(undesired_peak, h,gm,f,y,fs)

N=max(size(y));
desired_peak=[undesired_peak(end)+1:N];
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


function A=peak2peak_amp(y)
A = max(y) - min(y);

function invitro=Detect_invitro(f,y)
    
peak_y0=find(abs(y)==max(abs(y)));
y_slope0=diff(y);

%% check if the water peak is down not up
if y_slope0(peak_y0-1)<0
    y=-y;
end

y=normalized_vector(y);
y_slope=diff(y);


peak_y=find(y==max(y));
peak_diff=find(y_slope==max(y_slope));

slop_water=peak_y-peak_diff;

if slop_water==2
    invitro=1;
else
    invitro=0;
end

% figure(2);
% plot(y,'LineWidth',2);hold on
% plot(y_slope,'LineWidth',2);hold off
% d=1;


function y_n=normalized_vector(y)
y=y-min(y);
y_n=y./max(y);