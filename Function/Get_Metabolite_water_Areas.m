function [ind_Metabolite,ind_Water]=Get_Metabolite_water_Areas(freq,freq_vector)
freq_desired=freq_vector;
ind_Metabolite= find( freq_desired(1)<= freq & freq<=freq_desired(2));
ind_Water= find( -50<= freq & freq<=50) ;%  

