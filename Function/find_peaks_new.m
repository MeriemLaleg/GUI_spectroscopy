function [No_pks,pk_loc,FWHM,signal_power] = find_peaks_new(seg_sig,fres,msg)
% This function is used for peak selection: Click on the peak with mouse pointer to select the peaks
% Input:
% seg_sig: Segmented real part of the frequency domain signal
% Output:
% No_pks: Number of peaks selected
% pk_loc : Location of selected speaks in samples w.r.t segmented signal
% FWHM : Full Width Half maximum used as initial value for damping

% global error_types;
global details_para;
global fre_seg;
fre_seg=fres;
seg_sig=abs(fftshift(fft(seg_sig)));
sig_len = length(seg_sig);
rng  = details_para.seg;
while(1)
    
    
if(details_para.fla==0)
%plot signal in a full screen figure
h1 = figure('units','normalized','outerposition',[0 0 1 1],'Name','Peak Selection','NumberTitle', 'off'); 
plot(details_para.ppm_referenced',seg_sig);
% xlim([fre_seg(1), fre_seg(end)])
set(gca,'XDir','reverse','XMinorTick','on')

% Get the no of peaks to be selected
if strcmp(msg,'')
    No_pks=1;
else
    
variab = inputdlg(msg);
No_pks = str2double(variab);

end
close(h1)
elseif(details_para.fla==1)
No_pks=details_para.No_pks;
end 
% check whether the entered No_pks is within range
if(details_para.fla==0) 
    %plot signal in a full screen figure
    h2 = figure('units','normalized','outerposition',[0 0 1 1],'Name','Peak Selection','NumberTitle', 'on');
%     plot(ppm_seg,seg_sig)
%     xlim([ppm_seg(1), ppm_seg(end)])
   
    plot(seg_sig)
%     xlim([fre_seg(1), fre_seg(end)])
    set(gca,'XDir','reverse','XMinorTick','on')
    hold;
    plot(zeros(size(seg_sig)),'k')
end
   
    pk_plo = zeros(size(seg_sig)); % pk_plo will have values only at the peak locations and remaining values will be infinity, 
%     initially it will be all infinity, plotting pk_plo will plot only the values at peak location and infinity values are not plotted

    y1 = zeros(1,No_pks);
    for i = 1:No_pks % select peaks locations for the entered No_pks
            if(details_para.fla==0)
                [xi,~] = ginput(1); % get the x - location of mouse click
            elseif(details_para.fla==1)
                xi=details_para.xi(i);
            end
%             fre_pos = ((xi - details_para.ref)*details_para.Tf)/1E6; 
            fre_pos= xi
            vare = round( fre_pos);%round((fre_pos + (details_para.Fs/2))*(details_para.N/details_para.Fs) - rng(1) + 2);
%             vare = round(xi);
            if(vare>sig_len)
                vare = sig_len;
            elseif(vare<1)
                vare = 1;
            end
%             y1(i) = round(xi);
    y1(i) = vare;
    pk_plo(y1(i)) = seg_sig(y1(i)); 
    signal_power(i)=seg_sig(y1(i)); 
    if(details_para.fla==0)
        stem(pk_plo) % plot a stem at the selected peak location 
    end
    end
    if(details_para.fla==0)
    y_n_str = questdlg('Are Peaks Selected Properly');
    if(details_para.fla==0)
    if(strcmp(y_n_str, 'Yes')) 
        break
    elseif(strcmp(y_n_str, 'No')) % redo peak selection
        close(h2)
        continue;
    end
    end
    end
    if(details_para.fla==1)
        break
    end
end
if(details_para.fla==0)
    close(h2)
end

pk_loc = y1;



%% calculate FWHM

lower_hw = zeros(1,No_pks);
higher_hw = zeros(1,No_pks);
max_hwl = 10;
seg_sig=abs(seg_sig);
for i = 1:No_pks % loop for each seleced peak
    
    tmp1 = pk_loc(i);
    pk_val = seg_sig(tmp1);
    half_pk_val = pk_val/2;
    count = 1;
    lower_raising = 0;
    if(seg_sig(tmp1)<seg_sig(tmp1-1)) % Check whether the signal is raising on the lower side
        lower_raising = 1;
    end

    while(1) % find the half value location to the left (lower) side
        cur_pos = tmp1-count;
        if(cur_pos<1)
            lower_hw(i) = 1;
            break
        end
        val1 = seg_sig(cur_pos);
        if(val1<half_pk_val)
            lower_hw(i) = cur_pos;
            break
        end
        if(count == max_hwl)
            lower_hw(i) = cur_pos;
            break;
        end
        count = count + 1;
    end
    
    count = 1;
    upper_raising = 0;
    if(seg_sig(tmp1)<seg_sig(tmp1+1)) % Check whether the signal is raising on the upper side
        upper_raising = 1;
    end

    while(1) % find the half value location to the right (raising) side
        cur_pos = tmp1+count;
        if(cur_pos>sig_len)
            higher_hw(i) = sig_len;
            break
        end
        val1 = seg_sig(cur_pos);
        if(val1 < half_pk_val)
            higher_hw(i) = cur_pos;
            break
        end
        if(count == max_hwl)
            higher_hw(i) = cur_pos;
            break;
        end
        count = count + 1;
    end
    
    if((lower_raising == 1)&&(upper_raising == 0)) 
        lower_hw(i) = tmp1 - (higher_hw(i)-tmp1);
    elseif((lower_raising == 0)&&(upper_raising == 1))
        higher_hw(i) = tmp1 + (tmp1-lower_hw(i));
    end
end

fres_dif = fres(2)- fres(1); 
FWHM = (higher_hw-lower_hw)*fres_dif ;%convert sample difference to frequency difference
% FWHM
