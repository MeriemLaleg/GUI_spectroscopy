clear all; close all;
global  show_plot 
show_plot=1;              % Show plots
addpath ./Function

%% Load MRS spectrum from Data
load_MRS_data;    
input_freq_MRS=complex_fid_unsuppressed_FD(:,1);  %input_freq_MRS_ref=complex_fid_suppressed_FD(:,1); %% Load one voxel
input_freq_MRS_ref=input_freq_MRS; %% Load one voxel


%% ----------------------  Get the GUI Variable   -----------------------------------
% ppm=
% freq=
% freq_vector=
% input_freq_MRS=
% input_freq_MRS_ref=

% --------------  Run the SCSA water suppression algorithm   -----------------------
Nh_list=[14 24]; 
Scale_domain=[0 1];    % the projectin domaine
gm=0.5; fs=1;

[yf_SCSA,h_op,Nh_op]=SCSA_MRS_Water_Suppression_Scaling(Scale_domain, Nh_list, ppm,freq,freq_vector, input_freq_MRS,input_freq_MRS_ref,gm,fs);

%% --------------  The SCSA water suppressed spetrum is  <yf_SCSA> Variable   -----------------------------------------

clearvars -except ppm freq freq_vector  complex_fid_unsuppressed_FD yf_SCSA_Denoised  h_op  Nh
