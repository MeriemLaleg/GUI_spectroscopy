function [yf_wsUP,h_op,Nh_op]=SCSA_MRS_Water_Suppression_Scaling(Scale_domain,Nh_list,f,freq,freq_vector, complex_fid_unsuppressed_FD,complex_fid_suppressed_FD,gm,fs)
tic
global  show_plot yf0 

yf=real(complex_fid_unsuppressed_FD)';
yf0=real(complex_fid_suppressed_FD)';

%% Scale doen the Speectrum 
ydown=scale_signal(yf,Scale_domain);

% figure; plot(f, yf);hold on ;plot(f, yup+1000);hold on ;
% figure;plot(f, ydown);hold on 
% 

%% Search for the optimal number of eigenfunctions to remove the base line of metabolites 
yf2=ydown;
[ind_Metabolite,ind_Water]=Get_Metabolite_water_Areas(freq,freq_vector);

PSNR_op=-inf;

[h_list,Nh_list]= SCSA_1D_Nh_increase(yf2,fs,Nh_list,gm);

for h=h_list

    [yf_wp,squaredEIGF0,Sel_Eigf,Nh]=Supress_undesired_peaks(ind_Metabolite, h,gm,yf2,fs);

    %% Remore the water peak 
    yf_wp=yf_wp-min(yf_wp);
    yf2=yf2-min(yf2);
    yf_ws=yf2-yf_wp;

    PSNR=psnr(yf2,yf_wp);
%     yf_wpUP=scale_signal(yf_wp,yf);
%     yf_wpUP=yf_wpUP-min(yf_wpUP);
% 
%     PSNR=psnr(yf_wpUP,yf0);

    
    if PSNR_op<PSNR
        figure; plot(f, yf2);hold on  ; plot(f, yf_wp);hold on;hold on; plot( f, yf_ws);hold on;% plot( f, squaredEIGF0)% ;plot(SetInFrame(yf_scsa,yf));hold off; 
        legend('yf','Water peak phase','Suppressed Water peak  [Framed]')

        figure; plot(f, yf0);hold on  ; plot( f, yf_ws);hold on;% plot( f, squaredEIGF0)% ;plot(SetInFrame(yf_scsa,yf));hold off; 
        legend('yf0','Water peak phase 1','Suppressed Water peak  [Framed]')

        title(strcat('Estimatino using Nh_{wp}=',num2str(Nh),', PSNR=',num2str(PSNR)))

        PSNR_op=PSNR;
        h_op=h;
        Nh_op=Nh;
        yf_wp_op=yf_wp;
        
    end

    d=1;
end

yf_wpUP=scale_signal(yf_wp_op,yf);
yf_wpUP=yf_wpUP-min(yf_wpUP);
yf=yf-min(yf);
yf_wsUP=yf -yf_wpUP;

4.3- 5.1 
waterPeak=find(4.3<f & f<5.1);
yf_ws=Remove_base_no_metabolites(yf_wsUP,waterPeak, fs,gm);



%     figure; plot(f, yf);hold on  ; plot(f, yf_wpUP);hold on;hold on; plot( f, yf_ws);hold on;% plot( f, squaredEIGF0)% ;plot(SetInFrame(yf_scsa,yf));hold off; 
%     legend('yf','Water peak phase 1','Water peak phase [Framed]')
%     title(strcat('Estimatino using Nh_{wp}=',num2str(Nh),', PSNR=',num2str(PSNR)))
   
    
time=toc;
if show_plot==1

    plot(f,yf,'b','LineWidth',1.5); hold on
    plot(f, yf_ws ,'r','LineWidth',1.5)

    legend( 'Input MRS signal with water peak','Output MRS signal with suppressed water peak  using SCSA');
    xlabel('ppm')
    ylabel('Intensity')
%     set(gca,'YTickLabel',[])

    ymin=min(min([yf;yf_ws]'));  ymax=1.5*max(yf(ind_Metabolite));   
    ylim([ymin ymax])
    xlim([0 7])
    title(strcat(' SCSA MRS-water suppressionNh_{wp}=',num2str(Nh_op)))
    set(gcf,'color','w') 
    set(gca,'Xdir','reverse');
    set(gca,'fontsize',16)
    % text(1.9,90,'NAA');text(1.40,63,'Lac(1)');text(1.2,55,'Lac(2)');text(4.7,45,'Water residue');
    box 
    
    
end

 function yf_wp1=Remove_base_no_metabolites(yf_wp1,area,fs,gm)
    
padd=floor(max(size(area))/5); 
if padd > 10; padd=10; end


 N=max(size(yf_wp1));
 yf_No_Mtab=yf_wp1(area);  
 yff=padd_signal(yf_No_Mtab,padd);
 
Nh_No_Mtb= floor(N/100); % = 10; the reconstruct the global behavior     
[h3, yff_base,Nh_O]= SCSA_1D_Nh(yff,fs,Nh_No_Mtb,gm);
yf_22=unpadd_signal(yff_base,padd);
yf_22=yf_22-yf_22(1)+yff(1);
% yf_wp11(W0-5:N)=yf_22;

Res=abs(yf_No_Mtab-yf_22);
% %close all ;figure; plot(yf_No_Mtab);hold on;   plot(yf_22);

% Idx0=find(Res==min(Res))
yf_wp11=yf_wp1;
yf_wp11(area)=Res+yf_wp1(area(1));
%  close all ;figure; plot(yf_wp1);hold on;   plot(yf_wp11);hold on;   plot(yf_wp1-yf_wp11);
   
d=1;

yf_wp1=yf_wp11;
       

function [yd]=scale_signal(y,y0)
y=y-min(y);
y=y/max(y);
% min(y)
% max(y)
yd=y*(max(y0)-min(y0));
yd=yd+min(y0);
% min(yd)
% max(yd)
d=1;
