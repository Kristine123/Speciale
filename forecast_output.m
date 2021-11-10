function[OUT1]=forecast_output(data,cov,theta0,lb,opt,p,q,k,s,filename)
%function[OUT]=output(data,cov,theta0,lb,opt,p,q,s,filename)

%data = count data (length T = sample size + length of forecasting horizon) 
%cov = covariate matrix (each covariate is a column)
%theta0 = initial values for optimization
%lb = lower bound for optimization
%opt = optimization algorithm
%p = number of response lags
%q = number of intensity lags
%k = length of forecasting horizon
%s = xls file sheet 

T=length(data)-k;
dim=size(cov);
ncov=dim(2);
y(1:max(p,q))=0;
y(max(p,q)+1:T+max(p,q)+k)=data;
X(1:max(p,q),:)=0;
X(max(p,q)+1:T+max(p,q)+k,1:ncov)=cov;
TH=zeros(p+q+1+ncov,k);

y_hat=zeros(k+max(p,q));
y_hat(1:max(p,q))=y(T+1:T+max(p,q));


for i=1:k

y_sample=[y,y_hat(max(p,q)+i)];
TH(:,i)=forecast_estimate(y_sample,X,theta0,lb,opt,p,q,ncov,T-1+i);

omega=TH(1,i);
alpha=TH(2:p+1,i);
beta=TH(p+2:p+q+1,i);
gamma=TH(p+q+2:p+q+1+ncov,i);

y_hat(max(p,q)+i)=omega+sum(alpha'.*y(T+max(p,q)+i-1:-1:T+max(p,q)+i-p))+sum(beta'.*y_hat(max(p,q)+i-1:-1:max(p,q)+i-q))+sum(gamma'.*X(T+max(p,q)+i-1,:));
%y_hat(max(p,q)+i)=omega+sum(alpha'.*y(T+max(p,q)+i-1:-1:T+max(p,q)+i-p))+sum(beta'.*y_hat(max(p,q)+i-1:-1:max(p,q)+i-q));
end

OUT1(:,1)=y(T+max(p,q)+1:T+max(p,q)+k);
OUT1(:,2)=y_hat(max(p,q)+1:max(p,q)+k);

%xlswrite(filename,OUT1)
writematrix(OUT1,filename,'Sheet',s)
end




%xlswrite(filename,OUT,s)


