function data_load_suppress_init()

% This function is used for loading original data and consequently suppress using HLSVD and SCSA:  

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018

global details_para;
global data_value;
global hsvd_para;

details_para.error=0;
% .mat file loaded
try
    load(data_value.temp_filename);
catch
    details_para.error=1;
    return;
end
org_data_unsuppressed=complex_fid_unsuppressed;
for i=1:size(org_data_unsuppressed,2)
    data_value.FID_FD(:,i)=fftshift(fft((org_data_unsuppressed(:,i))));
end
details_para.Fs= Fs;
details_para.fres= (-(details_para.Fs)/2 + ((details_para.Fs)/(size(data_value.FID_FD,1)*1))*(0:((size(data_value.FID_FD,1)*1-1))));
details_para.t= (1:size(data_value.FID_FD,1))./details_para.Fs;
details_para.PE= 3;
details_para.RE=1;
details_para.Tf = 127776603;
details_para.ref = 4.7;
data_value.file_name{1}=[data_value.spect_name];
data_value.file_name{2}=[data_value.spect_name];
hsvd_para.order=5;

data_value.FID_TD=[org_data_unsuppressed]; % Add your SCSA output signal in time domain here in this matrix (follow the same dimensions), do ifft(ifftshift(your output signal)))
for i=1:size(org_data_unsuppressed,2)
    data_value.FID_FD(:,i)= [fftshift(fft((org_data_unsuppressed(:,i))))]; %Add your SCSA ouutput frquency domain signal here in this matrix (follow the same dimensions)
end
ppm = (details_para.fres)*(1E6/details_para.Tf);
details_para.ppm =  ppm;
details_para.ppm_referenced = details_para.ppm + details_para.ref ;
details_para.ds_flag=zeros(1,size(data_value.FID_TD,2));
details_para.fit_flag=zeros(1,size(data_value.FID_TD,2));
details_para.h_op=zeros(1,size(data_value.FID_TD,2));
details_para.preprocess_done=0;
details_para.fit_done=0;
details_para.ds_count=0;
details_para.ws_count=0;
details_para.preprocess_going=0;
details_para.fitting_going=0;