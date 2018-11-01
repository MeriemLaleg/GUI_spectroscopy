function [yd]=scale_signal(y,y0)
y=y-min(y);
y=y/max(y);

% min(y)
% max(y)


yd=y*(max(y0)-min(y0));
yd=yd+min(y0);


% min(yd)
% max(yd)

d=1;


