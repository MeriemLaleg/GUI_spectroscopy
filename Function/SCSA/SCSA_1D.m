

function [h, yscsa,Nh,psinnor,kappa,Ymin,squaredEIGF0]= SCSA_1D(y,fs,h,gm)

Lcl = (1/(2*sqrt(pi)))*(gamma(gm+1)/gamma(gm+(3/2)));
N=max(size(y));
%% remove the negative part
Ymin=min(y);

 y_scsa = y -Ymin;
%% Build Delta metrix for the SC_hSA
feh = 2*pi/N;
D=delta(N,fs,feh);

%% start the SC_hSA
Y = diag(y_scsa);
SC_h = -h*h*D-Y; % The Schrodinger operaor

% = = = = = = Begin : The eigenvalues and eigenfunctions
[psi,lamda] = eig(SC_h); % All eigenvalues and associated eigenfunction of the schrodinger operator
% Choosing  eigenvalues
All_lamda = diag(lamda);
ind = find(All_lamda<0);


%  negative eigenvalues
Neg_lamda = All_lamda(ind);
kappa = diag((abs(Neg_lamda)).^gm); 
Nh = size(kappa,1); %%#ok<NASGU> % number of negative eigenvalues



if Nh~=0
    
% Associated eigenfunction and normalization
psin = psi(:,ind(:,1)); % The associated eigenfunction of the negarive eigenvalues
I = simp(psin.^2,fs); % Normalization of the eigenfunction 
psinnor = psin./sqrt(I);  % The L^2 normalized eigenfunction 

squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;

%yscsa =4*h*sum((psinnor.^2)*kappa,2); % The 1D SC_hSA formula
yscsa1 =((h/Lcl)*sum((psinnor.^2)*kappa,2)).^(2/(1+2*gm));
else
    
  psinnor = 0*y;  % The L^2 normalized eigenfunction 
  kappa=0;
  squaredEIGF0=(h/Lcl)*(psinnor.^2)*kappa;

  yscsa1=0*y;
  yscsa1=yscsa1-10*abs(max(y));
%   disp('There are no negative eigenvalues. Please change the SCSA parameters: h, gm ')
  fprintf('\0')

end


if size(y_scsa) ~= size(yscsa1)
yscsa1 = yscsa1';
end
 
 %% add the removed negative part
 yscsa = yscsa1 + Ymin;






    %**********************************************************************
    %*********              Numerical integration                 *********
    %**********************************************************************

    % Author: Taous Meriem Laleg

    function y = simp(f,dx);
    %  This function returns the numerical integration of a function f
    %  using the Simpson method

    n=length(f);
    I(1)=1/3*f(1)*dx;
    I(2)=1/3*(f(1)+f(2))*dx;

    for i=3:n
        if(mod(i,2)==0)
            I(i)=I(i-1)+(1/3*f(i)+1/3*f(i-1))*dx;
        else
            I(i)=I(i-1)+(1/3*f(i)+f(i-1))*dx;
        end
    end
    y=I(n);
    

    %**********************************************************************
    %*********             Delata Metrix discretization           *********
    %**********************************************************************
    
    
%Author: Zineb Kaisserli

function [Dx]=delta(n,fex,feh)
    ex = kron([(n-1):-1:1],ones(n,1));
    if mod(n,2)==0
        dx = -pi^2/(3*feh^2)-(1/6)*ones(n,1);
        test_bx = -(-1).^ex*(0.5)./(sin(ex*feh*0.5).^2);
        test_tx =  -(-1).^(-ex)*(0.5)./(sin((-ex)*feh*0.5).^2);
    else
        dx = -pi^2/(3*feh^2)-(1/12)*ones(n,1);
        test_bx = -0.5*((-1).^ex).*cot(ex*feh*0.5)./(sin(ex*feh*0.5));
        test_tx = -0.5*((-1).^(-ex)).*cot((-ex)*feh*0.5)./(sin((-ex)*feh*0.5));
    end
    Ex = full(spdiags([test_bx dx test_tx],[-(n-1):0 (n-1):-1:1],n,n));
    
    Dx=(feh/fex)^2*Ex;
