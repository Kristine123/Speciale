function[OUT,pearson_residuals,pred,probability,meany,rpit,confidence_interval,theta]=output_parx_rob(y,X,theta0,p,q,s,filename)

%Numerical approximatation of maximum likelihood and restriction of outcome
%been none-negative

 theta_sqr0 = sqrt(theta0);
 %theta_sqr0 = theta0;   %to remove none negative restriction activate this line

options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun', 1e-8, 'TolX', 1e-8)
[theta_sqr,logLmin] = fminsearch('logL_PARX',theta_sqr0,options,y,X,p,q);
theta = theta_sqr.^2;
%theta = theta_sqr; %to remove none negative restriction activate this line

%ncov is specified as amount of rows and T is amount of columns
[ncov,T] =size(X);

omega = theta(1);
alpha = theta(2:p+1);
beta = theta(p+2:p+q+1);
gamma = theta(p+q+2:p+q+1+ncov);

%If X is empty zeros are filled in on empty spaces
if (ncov == 0)
    T = length(y);
    X = zeros(0,T);
    gamma = zeros(0,1);
end

%To derived inferrence the score and the hessian matrix is derived.
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

%this loop finish the calculation of score and hessian for the model including beta
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
    %This loop makes it possible to calculate with out including a beta
    %in the model
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

%log-likelihood and information critirias AIC and BIC is calculated
logl = -T*logLmin;
[aic,bic]=aicbic(logl,p+q+1+ncov,T);

%Variance and robust variance is calculated
global VAR
VAR = Omega^(-1)/T;
global VAR_ROB
VAR_ROB = H^(-1)*Omega*H^(-1)/T;

%Results are precented in matrix OUT
OUT=zeros(p+q+1+ncov+1,5);
OUT(p+q+1+ncov+1,1)=logl;
OUT(p+q+1+ncov+1,2:3)=[aic,bic];

%Variance is squared to obtain the standard deviation
for k=1:p+q+1+ncov
    OUT(k,1)=theta(k);
    OUT(k,2)=sqrt(VAR(k,k));
    OUT(k,3)=OUT(k,1)/OUT(k,2);   
    OUT(k,4)=sqrt(VAR_ROB(k,k));
    OUT(k,5)=OUT(k,1)/OUT(k,4);   

end

%Pearson residuals are calculated 
global pearson_residuals
for i=1:T
    pears(i)=(y(i)-lambda(i))/sqrt(lambda(i));
end
pearson_residuals = pears(p+1:T)';

%p is the amout of autoregressive lags 
%T is the length of time series

%Calculated lambda values in each point of the time series is stored as pred
global pred
for i=1:T
    pred(i)=lambda(i);
end
pred = pred(p+1:T)';

%Mean of y calculated at each point of the time series is stored as meany
global meany
    for i = 1:T
    meany(i) = mean(y(i));
    end
meany = meany(p+1:T)';    

%Calculates the probability of an specific integer beeing
%observed at a specific point in the time series
global probability
    for i=1:T-p 
    probability(i)=poisscdf(4,pred(i));
    end
probability = probability(1:T-p)';    
    
%Calculation of the Randomized Probabilti Integral Transform (PIT)
global rpit
    for i=1:T-p
    rpit(i)=poisscdf((y(i)-1),pred(i))+rand*(poisscdf(y(i),pred(i))-poisscdf((y(i)-1),pred(i)));
    end     
  rpit = rpit(1:T-p)';    


global confidence_interval
confidence_interval = confidence_interval';

%loops the specific time series
for i=1:size(pred)
    
    %Formular from
    %https://stats.stackexchange.com/questions/15371/how-to-calculate-a-confidence-level-for-a-poisson-distribution 
    confidence_interval(i) = pred(i) - 1.64 * (pred(i)/sqrt(i));
    %Ændret 01-09-2022 fra:
    %confidence_interval(i) = pred(i) - 1.64 * sqrt((pred(i)/i));
end

% %loops the specific time series
% for i=1:T-p 
%     hos = y(1:i)';
%     
%     %when the length of timeseries is only one,
%     %the next observation is considered additionally
%     %to be able to calculate a result
%     if(size(hos)==1)
%         hos = y(i:i+1)';
%     end
%     
%     %the poisson distribution is fitted
%     pd = fitdist(hos,'Poisson');
%     %the 95% confidence interval is calculated
%     ci = paramci(pd,'Alpha',.05);
%     %The calculated lower bound of the confidence interval is stored in
%     %cofidence_interval
%     confidence_interval(end+1) = ci(1,1);
% end


writematrix(OUT,filename,'Sheet',s)
end





