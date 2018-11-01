    function [FitParams, rejectframe, residCr]  = FitPeaks(freq, FrameData, initx)
    
    
    nlinopts = statset('nlinfit');
    nlinopts = statset(nlinopts, 'MaxIter', 1e4, 'Display','Off');
    nframes = size(FrameData,2);
    FrameData=real(FrameData);
    for jj = 1:nframes
        [fit_param, residCr] = nlinfit(freq', (FrameData(:,jj)), ...
            @(xdummy, ydummy) LorentzModel(xdummy, ydummy), ...
            initx, nlinopts);
        FitParams(jj,:) = fit_param;
        FitParams2(jj,:) = fit_param;
        fit_plot = LorentzModel(fit_param, freq);
    end
    
    for kk=1:size(FitParams,1)
        if FitParams(kk,1)<0
            FitParams(kk,4)= FitParams(kk,4)+pi;
        end
    end
    
    % Need to deal with phase wrap:
    % Convert to complex number then recalculate phase within 2*pi range
    phase_wrapped = FitParams(:,4);
    cmplx = cos(phase_wrapped) + 1i * sin(phase_wrapped);
    phase_unwrapped = angle(cmplx);
    

    % then fix to be within -pi..pi
    offsetpos =  2*pi*lt(phase_unwrapped, -pi);
    offsetneg = -2*pi*gt(phase_unwrapped,  pi);
    phase_unwrapped = phase_unwrapped + offsetpos + offsetneg;
    FitParams(:,4) = phase_unwrapped;
    
    % Fix area and linewidth to be positive

    FitParams(:,1) = abs(FitParams(:,1));
    FitParams(:,2) = abs(FitParams(:,2));
    
    end
