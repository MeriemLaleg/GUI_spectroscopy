function [amp,phas]=phas_correct(amp_phs,order)

x=real(amp_phs); y=imag(amp_phs);
for i=1:order,
  if x(i)==0,
    if y(i)>=0,
      phas(i)=90;
    else
      phas(i)=270;
    end
  else,
    phas(i)=atan(y(i)/x(i))*180/pi;
  end
end

% the following makes phas correct and positive : from 0 to 360- 
for i=1:order,
  if y(i)<0 & x(i)>0;
    phas(i)=phas(i)+360;
  end
  if x(i)<0;                      % necessary when x<0 and y>0.
    phas(i)=phas(i)+180;
  end
end
% to change phas from 0 -- 360 to -180 -- 180-, use phaseadj or 
% delete the above lines and change them to if x<0 & y>0, phas+180, end
phas=phas*pi/180;
amp=abs(amp_phs).*exp(1i*phas');