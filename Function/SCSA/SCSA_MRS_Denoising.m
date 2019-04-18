%% MRS signal Denosiing  using SCSA Method : 
% This function denoised the real part of the MRS signal using and optimal
% value of the h the ensures the preservation of the metabolite area,
% determined by <Metabolite>, while remove as maximuim as possible the noise from the
% region determined by <Noise>

%% ######################  PARAMETERS ######################################
% Input 
% ppm         : The MRS spectrum fequancies in ppm
% yf          : Noisy real art of the  complex MRS FID signal
% Th_peaks_ratio  : The ratio w.r.t metabolite amplitude to set the threshold from which the peak can be selected
% width_peaks  : The with of the peak from its center(max values)
% gm, fs  : SCSA parameter: gm=0.5, fs=1

% Output
% yscsa: Densoied real art of the  complex MRS FID signal
% h_op : The  real art of the  complex MRS FID signal
% Nh   : Densoied real art of the  complex MRS FID signal

%% ###########################################################################
%  Author:
%  Abderrazak Chahid (abderrazak.chahid@gmail.com)
%  Adviser:
%  Taous-Meriem Laleg (taousmeriem.laleg@kaust.edu.sa)
% Done: June,  2018
% King Abdullah University of Sciences and Technology (KAUST)

function [yscsa, h_op, Nh]=SCSA_MRS_Denoising(SNR_regularization, ppm, yf, gm , fs , Th_peaks_ratio, width_peaks)
% close all
global shif show_plot  

%% Get the noise and metabolites locations
[Noise,Metabolite]=find_metabolite_Noisy_Area(yf, Th_peaks_ratio, width_peaks);

Noise(2)=Noise(1)+10;
%% Generate signals
fprintf('\n_____________________________________________________________')
fprintf('\n MRS signal denoising using SCSA  (By Abderrazak - KAUST- 2018)')
fprintf('\n_____________________________________________________________\n\n')


%% MRS signal 
y = real(yf) ;%/max(real(yf))*76;

y_positive=y-min(y);
h_max=0.5*(max(y_positive) + mean(y_positive));
if max(y_positive)>1
    h_min=0.1;
else
    h_min=max(y_positive)/10;
end

fprintf('\n-->Searching for the optimal h. Please wait...')

% you can choose either a vecotr for h or one value 
% h_list=h_min : .1 :h_max ;  % search in this interval
% h_list=10:0.1:17;
eps=0.2;
% shif=0;close all;
Max_empty_iter=10;     % the the maximuim iteration allows with increasing cost function [H goes far from the optimal h]


%% Setup the  Wait bar
global cnt_wait Tot_iter   wait_bar
Tot_iter=0;

M =3;                       % number of refinement loop to find the optimal h for denoising
for k=1:M
    Tot_iter= Tot_iter+floor(200/2^k);
end

cnt_wait=0;
wait_bar = waitbar(cnt_wait,'Searching for the optimal value of h .... ','Name','MRS signal denoising using SCSA');


%% Search for the optimal h 
for k=1:M
    h_list= linspace(h_min,h_max,floor(200/2^k));
    [h_op,Min_Cost]=search_4_optimal_h_denoising(Metabolite, Noise,SNR_regularization, y,fs,h_list,gm);
    h_min=h_op*(1-eps);
    h_max=h_op*(1+eps);
    %% plot the results     
    h_op_list(k)=h_op;
    Cost_function(k)=Min_Cost;
    
end
h_op=min(h_op_list(find(Cost_function==min(Cost_function))));
[h_op, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h_op,gm);

fprintf('\n--> MRS denoising is completed h=%f!!',h_op)
close(wait_bar)


 %% Get the optimal     
% figure(2)
% plot(Cost_function); hold on
% plot(Cost_func_SNR); hold on
% plot(Cost_func_MSE); hold on
% 
% xlabel('Iterations')
% ylabel('Cost function ')
% title([ 'The Iterative MRS signal denoising']);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])


%% Plot the denoising results 


if show_plot==1
      figure(1);
    if shif==0
        plot(ppm,y,'b','LineWidth',1);
    end
    hold on
    plot(ppm, yscsa+shif ,'LineWidth',1)

    % hold on
    % plot(ppm, y-yscsa-5,'g','LineWidth',1)
    legend({'Noisy input spectrum ', 'Denoised Spectrum '},'Location','northwest');
    xlabel('ppm')
    ylabel('Intensity')
    set(gca,'YTickLabel',[])
%     xlim([0 5])
    title([ ' SCSA MRS denoising with : h = ' num2str(h_op) , '   Nh = ' num2str(Nh) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
    set(gcf,'color','w') 
    set(gca,'Xdir','reverse');
    set(gca,'fontsize',16)
    % text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
    box 
end
    


function  [h_op,Min_Cost]=search_4_optimal_h_denoising(Metabolite, Noise,SNR_regularization, y,fs,h_list,gm)
global Tot_iter  wait_bar   cnt_wait
% start the loop for several values of h 
 for i=1:length(h_list)
 
    h = h_list(i);%0.045;%15.2;

    [Cost_SNR,Cost_peak,h, yscsa,Nh]=Cost_function_4_denoising(Metabolite, Noise,SNR_regularization, y,fs,h,gm);
    Cost_func_SNR(i)=Cost_SNR;
    Cost_func_MSE(i)=Cost_peak;
    Cost_function(i)=Cost_SNR+Cost_peak;

    fprintf('.')
    cnt_wait=cnt_wait+1;
    % Update waitbar and message
    waitbar(cnt_wait/Tot_iter,wait_bar)
%     figure(154);
%     plot(ppm,y,'LineWidth',2.5);hold on
%     plot(ppm, yscsa ,'LineWidth',1.2);hold off
% 
%     % hold on
%     % plot(f, y-yscsa-5,'g','LineWidth',1)
%     legend({'Noisy input spectrum ', 'Denoised Spectrum ','Residue'},'Location','northwest');
%     xlabel('ppm')
%     ylabel('Intensity')
%     set(gca,'YTickLabel',[])
%     xlim([0 5])
%     title([ ' SCSA MRS densoing with : h = ' num2str(h) , '   Nh = ' num2str(Nh) ]);% ', fe = ' num2str(fe) ', amp noise = ' num2str(amp_noise)])
%     set(gcf,'color','w') 
%     set(gca,'Xdir','reverse');
%     % text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
%     box 
    
    
 end
Min_Cost=min(Cost_function);
h_op=min(h_list(find(Cost_function==min(Cost_function))));

function [noise_area,peaks_Locations]=find_metabolite_Noisy_Area(y0, Th_peaks_ratio, width_peaks)

%% Find tw peaks location 
N=max(size(y0));yf=y0(1:floor(0.48*N));
TH=median(yf)+0.12*mean(yf(1:10));
yf_peaks=0*yf+TH;
yf_peaks(find(yf>TH))=yf(find(yf>TH));
yf_smooth   = sgolayfilt(yf_peaks, 2, 7);
Th=(max(yf_smooth)-min(yf_smooth))/Th_peaks_ratio +min(yf_smooth);



[pks,locs,w,p] = findpeaks(yf_peaks);

peaks_location=locs(find(yf_peaks(locs) > 0.7*mean(yf_peaks(locs))));
Peaks_areas=0;

Peak_TH=peaks_location(find( yf_smooth(peaks_location)>Th));

for i=1:max(size(Peak_TH))
   
    peaks_Locations(i,:)=[Peak_TH(i)-width_peaks,Peak_TH(i)+width_peaks];    
end


for i=1:max(size(Peak_TH))
   
    Peaks_areas=[Peaks_areas,peaks_Locations(i,1):peaks_Locations(i,2)];
    
end

Metabolite_area=1:peaks_location(end);

Peaks_areas=unique(Peaks_areas(2:end));
noise_areas=setdiff(Metabolite_area,Peaks_areas);
noise_area=[noise_areas(1),noise_areas(end)];


function [Cost_SNR,Cost_peak,h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]=Cost_function_4_denoising(Metabolite, Noise,SNR_regularization, y,fs,h,gm)

        % yscsa = scsa(h);
[h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm);

% coco1(i) = norm(y(a1:b1) - yscsa(a1:b1));
% coco2(i) =  norm(y(a2:b2) - yscsa(a2:b2));
% coco3(i) =  norm(y(a3:b3) - yscsa(a3:b3));
% % coco(i) =norm(y - yscsa);

y_res=y-yscsa;

SNR_y=max(yscsa)/std(y(Noise(1):Noise(2)));   %second term retaed to the  SNR
SNR_scsa=max(yscsa)/std(yscsa(Noise(1):Noise(2)));   %second term retaed to the  residual SNR

Cost_SNR=(SNR_regularization/abs(SNR_scsa));%*abs(SNR_y);

Cost_peak=0;
for k=1:size(Metabolite,1) 
    Cost_peak=Cost_peak + norm(y_res(Metabolite(k,1:2)));
end


    

% figure;plot(yf_peaks);hold on;plot(yf_smooth);hold on;plot(diff(yf_smooth));hold on;plot(yf_smooth*0+Th);

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


%yscsa =4*h*sum((psinnor.^2)*kappa,2); % The 1D SC_hSA formula
yscsa1 =((h/Lcl)*sum((psinnor.^2)*kappa,2)).^(2/(1+2*gm));
else
    
  psinnor = 0*psi;  % The L^2 normalized eigenfunction 

  yscsa1=0*y;
  yscsa1=yscsa1-10*abs(max(y));
  disp('There are no negative eigenvalues. Please change the SCSA parameters: h, gm ')
end


if size(y_scsa) ~= size(yscsa1)
yscsa1 = yscsa1';
end
 
 %% add the removed negative part
 yscsa = yscsa1 + Ymin;


squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;




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

    


