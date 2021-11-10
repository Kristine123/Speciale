function lambda_forecast = PARX_forecast(theta,Y,X,p,q)
% Hvad gør funktionen?

[logL,lambda] = logL_PARX(theta,Y,X,p,q);

[d,T] = size(X);
if (d == 0)
    T = length(Y);
    X = zeros(1,T);
end

omega = theta(1)^2;
alpha = theta(2:p+1).^2;
beta = theta(p+2:p+q+1).^2;
gamma = theta(p+q+2:p+q+1+d).^2;

if (d == 0)
    gamma = 0;
end

lambda_forecast = omega + alpha'*Y(T:-1:T+1-p)' + beta'*lambda(T:-1:T+1-q)' + gamma'*X(:,T);