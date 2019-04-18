function data_load_suppress()

% This function is used for loading original data and consequently suppress using HLSVD and SCSA:  

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018

global details_para;
global data_value;
global hsvd_para;

try
    if (data_value.load_format==1)
        load(data_value.temp_filename);
    elseif(data_value.load_format==2)
        siemens_rda_read(data_value.temp_filename);
    elseif(data_value.load_format==3)
        readSiememsDicom(data_value.temp_filename);
    elseif(data_value.load_format==4)
        [folder, baseFileName, ~] = fileparts(data_value.temp_filename);
         read_philips_file([folder,'\',baseFileName]);
    end 
catch
    details_para.error=1;
    return;
end

preprocessing_SCSA_HSVD();
ppm = (details_para.fres)*(1E6/details_para.Tf);
details_para.ppm =  ppm;
details_para.ppm_referenced = details_para.ppm + details_para.ref ;


