function [h_list,Nh_list]= SCSA_1D_Nh_increase(y,fs,Nh_list,gm)

global  Plot_fig

h=max(y)-min(y);
[h, yscsa,Nh]= SCSA_1D(y,fs,h,gm);

Nh0=Nh_list(1);
while Nh/Nh0~=1

    if Nh==0
        h=h/2;
    else
    h=h*Nh/Nh0;
    end
    [h, yscsa,Nh]= SCSA_1D(y,fs,h,gm);
    d=1;
end

h0=h*Nh/Nh_list(2);

[h0, yscsa0,Nh0,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h0,gm);

Dh=(h-h0)/(Nh-Nh0);

h_list=h:Dh/4:h0;
Nh_list=[Nh:Nh0];



