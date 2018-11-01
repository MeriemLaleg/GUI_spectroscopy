clear all; close all; 
addpath ./Function
load_MRS_data

%  load('Example_MRS_data.mat');input_freq_MRS=complex_fid_suppressed_FD;

global  show_plot 
show_plot=1;

gm=0.5;fs=1;        % SCSA parameters
Nh_list=[15 24];    % defaut range of number of eigenfunctions

[yf_wsUP,h_op,Nh_op]=SCSA_MRS_Water_Suppression_Scaling(Nh_list, ppm,freq,freq_vector, input_freq_MRS,gm,fs);
































% %% Location where there is almost noise
% m1= 368; m2=497;
% 
% %% Location of metabolites 
% a1= 141;b1= 149; % NAA
% a2= 168;b2= 172; % Lac
% a3= 218;b3= 222; % Chr
% gm=0.5; fs=1;
% %% Noise and metabolites peak variable. needed for the optimization
% Noise=[m1,m2];
% Metabolite=[a1  b1; a2  b2; a3  b3];
% 
% [yf_SCSA_Denoised, h_op, Nh]=SCSA_MRS_Denoising(ppm, yf_SCSA, Metabolite, Noise, gm , fs );

clearvars -except ppm freq freq_vector  complex_fid_unsuppressed_FD yf_SCSA_Denoised  h_op  Nh
