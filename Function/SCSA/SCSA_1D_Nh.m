function [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D_Nh(y,fs,Nh0,gm)

h=max(y)-min(y);

[h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm);

Dh=0;Dh2=0;

while Nh/Nh0~=1
    h_old=h;
    Dh_old=Dh;Dh2_old=Dh2;

    
    if Nh==0
        h=h/2;
    else
    h=h*Nh/Nh0;
    end
    
    
    [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm);
    
    Dh=h_old-h;Dh2= abs(Dh-Dh_old);
    Dh_loop= abs(Dh2-Dh2_old);
    
    if Dh_loop==0 , h=h*0.923; end
    
    d=1;
end

