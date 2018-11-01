function [sm_signal,function_weight] = apodize(FID,type,para)
% apodization
% Input:
% FID= Time domain signal
% type= Apodization type: 1-Lorentz; 2-Gauss; 3-Lorentz-Gauss; 4-Sigmoid
% para= Apodization parameters: para(1)-sampling frequency; para(2)-lorentz para; para(3)-Gauss para(4)- sigmoid para
% Output:
% sm_signal= Signal after apodization in time domain
% w = weight vector used for apodisation 


[N, N_vox] = size(FID);
n_t = ((0:N-1))/para(1);
n = 0:N-1;
function_weight= ones(N,N_vox);

switch type
    case 1 % Exponential
        if(para(2) == 0)
            sm_signal = FID;
            return
        else
            a= (pi*para(2));
            w = exp((-a)*n_t);
%           w = exp(-para(2)*n_t);

        end
    case 2 % Gauss
        if(para(3) == 0)
            sm_signal = FID;
            return
        else
            a= (pi*para(3))/2*sqrt(log(2));
            w = exp(-((n_t).*(n_t))*(a^2));
%             w = exp(-((n_t).*(n_t))*para(3));
%             a = 4*log(2)/(para(3)^2);
%             w = sqrt(a/pi)*exp(-a*n_t.^2)

        end
    case 3 % Gauss-Exponential
        if((para(2) == 0) && (para(3) == 0))
            sm_signal = FID;
            return
        else
            a1 = (pi*para(2)); a2 = (pi*para(3))/2*sqrt(log(2));
            w = exp(-(a1*n_t + (a2^2)*(n_t.*n_t)));
        end
        
    case 4 % Sigmoid
        if(para(4) == 0)
            sm_signal = FID;
            return
        else
            w = 1./(1+exp(0.05*(n - para(4))));
        end
        
    otherwise % do Exponential 
        w = exp(-para(2)*n/para(1));    
end

 weight_mat = repmat(w',[1,N_vox]);
 function_weight = weight_mat;
 sm_signal = weight_mat.*FID; % multiply the weight matrix with the signal matrix to perform apodization
 