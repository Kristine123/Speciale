function [logL,lambda] = logL_PARX_saturday(theta,Y,X,p,q)

[d,T] = size(X);

omega = theta(1)^2;
alpha = theta(2:p+1).^2;
beta = theta(p+2:p+q+1).^2;
gamma = theta(p+q+2:p+q+1+d).^2;

if (d == 0)
    T = length(Y);
    X = zeros(1,T);
    gamma = 0;
end

lambda = zeros(1,T);
lambda(1:max(p,q)) = mean(Y)*ones(1,max(p,q));

for t = max(p,q)+1:T
    lambda(t) = omega + alpha'*Y(t-1:-1:t-p)' + beta'*lambda(t-1:-1:t-q)' + gamma'*X(:,t-1);
end

logf = Y.*log(lambda) - lambda;
logL = -sum(logf(p+1:T))/T; 