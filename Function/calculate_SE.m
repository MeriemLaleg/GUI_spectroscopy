%% Calculation of standard error

% Input:
% beta: estimated parameters
% resid: Residual signal
% J: Jacobian matrix computed form NLLS
% Output:
% se: Standard error

function se = calculate_SE(beta,resid,J)

n = length(resid);
p = numel(beta);
v = n-p;

[~,R] = qr(J,0);   % upper triangualr matrix of the Jacobian computed from the residual 
Rinv = R\eye(size(R)); 
diag_info = sum(Rinv.*Rinv,2); % Computation of Cramer Rao lower bound
% diag_info=diag(pinv(full(J'*J)));

rmse = (norm(resid) / sqrt(v)); %s.d. of residual matrix
se = sqrt(diag_info) * rmse; % computation of standard error 
% se= se.^2; 
