function [theta_s] = admm_knnfl_varying_rho(X, Y,W, lambda, k, rho, max_iter,mu,tau_inc,tau_dec)
% Inputs:
% X: n x d feature matrix
% Y: n x 1 target variable vector
% lambda: regularization parameter
% k: number of nearest neighbors for the knn graph
% rho: ADMM penalty parameter
% max_iter: maximum number of iterations
% 
% Outputs:
% beta: d x 1 weight vector

n = size(X, 1);
d = size(X, 2);

% Initialize variables
theta_s = zeros(n, 1);
z = zeros(n, 1);
u = zeros(n, 1);


%initialize the residuals
s=0;
r=0;

% ADMM iterations
for iter = 1:max_iter
    % Update beta
    for i=1:n
        theta_s(i)=(-u(i)-rho*(Y(i)-z(i)))/(rho+2*(W(i)))+Y(i);
    end
    % Update z
    z_old=z;
    param_knfl=theta_s+u/rho;
    new_lambda=2*lambda/rho;
    z=knnfl(X.',param_knfl.',k,new_lambda);
    z=z.';
    
    %update rho
    s=rho*(z-z_old);
    r=theta_s-z;
    if norm(r) > mu*norm(s)
        rho=rho*tau_inc;
    elseif norm(s) > mu*norm(r)
        rho=rho/tau_dec;
    else
        rho=rho;
    end

    % Update u
    u = u + rho*(theta_s - z);
    % Check convergence
    %if norm(z - z_old) / norm(z_old) < 1e-3
     %   break;
   % end
end

end

