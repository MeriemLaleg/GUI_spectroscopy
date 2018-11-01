
 global   cname ind_Metabolite ind_Water
%% Generate signals

%% #########################    Load data   ################################
% ext = './Input_data/Sourav/18-03-31-New_MRSdata_1_4_17/*.mat';  
% ext = './Input_data/Sourav/18-04-08_DataSet_Invivo_Simulated_Sourav/*.mat';  
% ext = './Input_data/Sourav/18-04-19-DataSet3_Invivo_Simulated_Sourav/*.mat';  
ext = './Input_data/*.mat';  

[filename rep]= uigetfile({ext}, 'File selector')  ;
chemin = fullfile(rep, ext);  list = dir(chemin);  
cname=strcat(rep, filename);  
load(cname); 

%% Prepare data

% if real_abs==0
    yf=real(complex_fid_unsuppressed_FD)';
% else
%     yf=abs(complex_fid_unsuppressed_FD)';
%   
% end
if exist('complex_fid_suppressed')==1
    
%     if real_abs==0
        yf0=real(complex_fid_suppressed_FD)';
%     else
%         yf0=abs(complex_fid_suppressed_FD)';
% 
%     end



else
yf0=yf;
complex_fid_suppressed_FD=complex_fid_unsuppressed_FD;
    
end


yt=complex_fid_unsuppressed';

if exist('complex_fid_suppressed')==1
    yt0=complex_fid_suppressed';
else
    yt0=yt;
    complex_fid_suppressed=complex_fid_unsuppressed;

end
%% Change notations
N=max(size(yf));
f=ppm;

if exist('peaks_pos')~=0
ind_Water= peaks_pos(end-1):peaks_pos(end) ;% 
ind_Metabolite=peaks_pos(1):peaks_pos(end-2);

else
    freq_desired=freq_vector;
    ind_Metabolite= find( freq_desired(1)<= freq & freq<=freq_desired(2));
    ind_Water= find( -50<= freq & freq<=50) ;%  
    peaks_pos=sort([ind_Metabolite(1) ind_Metabolite(end) ind_Water(1) ind_Water(end)]);
end

peaks_pos_ref=peaks_pos;
Counting=0;
name_peaks={'Lac-NAA','Water'};
peaks_time=ppm(peaks_pos);
N_peaks=max(size(peaks_pos))/2 -1;
ppm_Max=ppm(peaks_pos(end-2)) ;
ppm_zoom=max(find(f<=ppm_Max));
ppm_zoom2=floor(peaks_pos(end-2)*1.3);

%% #Evaluation of SNR of the noisy image/signal
% yf=yf-min(yf)+min(yf0);

name_data=filename(1:end-4); %noisy_file;
fprintf('\n--> Load Data: %s',name_data )

% close all;plot(f,yf)

input_freq_MRS=complex_fid_unsuppressed_FD(:,1);
