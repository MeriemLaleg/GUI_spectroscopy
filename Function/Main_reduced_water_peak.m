clear all; close all;
global  show_plot  Loop_Water_reduction Nh_WP_reduction
show_plot=1;              % Show plots
addpath ./Function

%% Load MRS spectrum from Data
load_MRS_data;    

input_freq_MRS=complex_fid_unsuppressed_FD(:,1);  %input_freq_MRS_ref=complex_fid_suppressed_FD(:,1); %% Load one voxel
N=max(size(input_freq_MRS));

%% Reduced water peak - invivo and simulated MRS data
gm=0.5;fs=1;                        % SCSA parameters

Nh_WP=20;                           % the minimum number of eigenfunction. For simulated=20, and for invivo=13;    
Loop_Water_reduction=4;             %  It reduces the water residual. The values range is   from 1 to 4
Nh_WP_reduction= 4;                 % floor(N/100): It make the water residual more flat . The values range is   from 1 to 3

yf_SCSA=SCSA_MRS_Water_Suppression_GUI(ppm,freq,freq_vector, input_freq_MRS,gm,fs,Nh_WP);

