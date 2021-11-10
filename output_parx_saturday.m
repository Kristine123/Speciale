%function[OUT] = output_parx_robb(y,X,theta0,p,q,s,filename)
%function[OUT,pearson_residuals,pred,zero_prob,meany,rpit,upper_pred]=output_parx_rob(y,X,theta0,p,q,s,filename)
function[OUT,pearson_residuals]=output_parx_saturday(y,X,theta0,p,q,s,filename)
% data = sample
% cov = covariate matrix (each covariate is a column)
% theta0= imitial values for optimization
% lb = lower bound for optimization
% opt = optimization algorithm
% p = number of response lags
% q = number of intensity lags
% s = xls file sheet


theta_sqr0 = sqrt(theta0);

options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun', 1e-8, 'TolX', 1e-8,'Display', 'iter','PlotFcns',@optimplotfval)
[theta_sqr,logLmin] = fminsearch('logL_PARX_saturday',theta_sqr0,options,y,X,p,q);

theta = theta_sqr.^2;

[ncov,T] =size(X);

omega=theta(1);
alpha=theta(2:p+1);
beta=theta(p+2:p+q+1);
gamma=theta(p+q+2:p+q+1+ncov);

if (ncov == 0)
    T = length(y);
    X = zeros(0,T);
    gamma = zeros(0,1);
end

global lambda
lambda = zeros(1,T);
score = zeros(1+p+q+ncov,1);
hessian = zeros(1+p+q+ncov,1+p+q+ncov);

Dlambda = zeros(1+p+q+ncov,1);
Dlambda_oa = zeros(1+p,1);              %1st order deriv of lambda w.r.t. (omega,alpha)
Dlambda_g = zeros(ncov,1);              %1st order deriv of lambda w.r.t. gamma
Dlambda_beta = zeros(q,1);               %1st order deriv of lambda w.r.t. beta
D2lambda = zeros(p+q+1+ncov,p+q+1+ncov);
D2lambda_oa_oa = zeros(1+p,1+p);         %2nd order deriv of lambda w.r.t. (omega,alpha)
D2lambda_oa_g = zeros(1+p,ncov);        %2nd order deriv of lambda w.r.t. (omega,alpha)|gamma
D2lambda_oa_b = zeros(1+p,q);           %2nd order deriv of lambda w.r.t. (omega,alpha)|beta
D2lambda_g_g = zeros(ncov,ncov);        %2nd order deriv of lambda w.r.t. gamma
D2lambda_g_b = zeros(ncov,q);           %2nd order deriv of lambda w.r.t. gamma|beta
D2lambda_b_b = zeros(q,q);               %2nd order deriv of lambda w.r.t. beta

Omega = zeros(p+q+1+ncov,p+q+1+ncov);
H = zeros(p+q+1+ncov,p+q+1+ncov);
                    
lambda(1:max(p,q)) = mean(y)*ones(1,max(p,q));
if (q > 0)
    for i = max(p,q)+1:T    
        lambda(i) = omega + alpha'*y(i-1:-1:i-p)' + beta'*lambda(i-1:-1:i-q)' + gamma'*X(:,i-1);

        D2lambda_oa_b = Dlambda_oa + beta*D2lambda_oa_b;
        D2lambda_g_b =  Dlambda_g + beta*D2lambda_g_b;
        D2lambda_b_b = 2*Dlambda_beta + beta*D2lambda_b_b ;     

        D2lambda = [D2lambda_oa_oa,D2lambda_oa_b,D2lambda_oa_g; ...
                    D2lambda_oa_b',D2lambda_b_b,D2lambda_g_b'; ...
                    D2lambda_oa_g',D2lambda_g_b,D2lambda_g_g];

        Dlambda_oa = [1;y(i-1:-1:i-p)'] + beta'*Dlambda_oa;
        Dlambda_g = X(:,i-1) + beta*Dlambda_g;
        Dlambda_beta = lambda(i-1) + beta*Dlambda_beta;              
        Dlambda = [Dlambda_oa;Dlambda_beta;Dlambda_g];

        score =  (y(i)/lambda(i) - 1)*Dlambda;
        hessian = (y(i)/lambda(i) - 1)*D2lambda - y(i)/lambda(i)^2*Dlambda*Dlambda';

        Omega = Omega + score*score';
        H = H + hessian;
    end
else
    for i = max(p,q)+1:T    
        lambda(i) = omega + alpha'*y(i-1:-1:i-p)' + gamma'*X(:,i-1);

        Dlambda = [1;y(i-1:-1:i-p)';X(:,i-1)];

        score = (y(i)/lambda(i) - 1)*Dlambda;
        hessian = -y(i)/lambda(i)^2*Dlambda*Dlambda';

        Omega = Omega + score*score';
        H = H + hessian;
    end
end
Omega = Omega/T;
H = H/T;

logl = -T*logLmin;
[aic,bic]=aicbic(logl,p+q+1+ncov,T);
global VAR
VAR = Omega^(-1)/T;
global VAR_ROB
VAR_ROB = H^(-1)*Omega*H^(-1)/T;

OUT=zeros(p+q+1+ncov+1,5);
OUT(p+q+1+ncov+1,1)=logl;
OUT(p+q+1+ncov+1,2:3)=[aic,bic];


for k=1:p+q+1+ncov
    OUT(k,1)=theta(k);
    OUT(k,2)=sqrt(VAR(k,k));
    OUT(k,3)=OUT(k,1)/OUT(k,2);   
    OUT(k,4)=sqrt(VAR_ROB(k,k));
    OUT(k,5)=OUT(k,1)/OUT(k,4);   

end

global pearson_residuals
for i=1:T
    pears(i)=(y(i)-lambda(i))/sqrt(lambda(i));
end
pearson_residuals = pears(p+1:T)';

writematrix(OUT,filename,'Sheet',s)
end





