function [td_sig] = phase_correction(sig,phi_0,phi_1,fres)
% Performs zero order and first order phase correction
% input: 
% sig - frequency domain signal without phase correction
% phi_0 - zero order phase correction value
% phi_1 - first order phase correction value
% output:
% td_sig - time domain phase corrected signal
a_num = fres;

[N,N_vox] = size(sig); % get the no of samples of each signal and no of signals
% a_num = (0:N-1);
% tot_pha = phi_0 + 2*pi/N*((-phi_1*a_num)); %get the phase correction values
% a_num = [linspace(0,(N-1)/N,N/2) linspace(-(N-1)/N,0,N/2)];

tot_pha = phi_0 + (2*pi*-phi_1*10^-3).*a_num; 
phase_mat = repmat(tot_pha',[1,N_vox]); %generate phase correction matrix
af_pc = (exp(1i*phase_mat).*sig); %do phase correction 
% td_sig = ifft(af_pc);
td_sig = ifft(ifftshift(af_pc));
% fd_sig = fftshift(af_pc,1);


