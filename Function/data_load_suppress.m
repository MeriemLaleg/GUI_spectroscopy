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
    load(data_value.temp_filename);
catch
    details_para.error=1;
    return;
end

preprocessing_SCSA_HSVD();
ppm = (details_para.fres)*(1E6/details_para.Tf);
details_para.ppm =  ppm;
details_para.ppm_referenced = details_para.ppm + details_para.ref ;


