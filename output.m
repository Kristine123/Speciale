function[OUT]=output(data,cov,theta0,lb,opt,p,q,s,filename)

% data = sample
% cov = covariate matrix (each covariate is a column)
% theta0= imitial values for optimization
% lb = lower bound for optimization
% opt = optimization algorithm
% p = number of response lags
% q = number of intensity lags
% s = xls file sheet



T=length(data);
dim=size(cov);
ncov=dim(2);
y(1:max(p,q))=0;
y(max(p,q)+1:T+max(p,q))=data;
X(1:max(p,q),:)=0;
X(max(p,q)+1:T+max(p,q),1:ncov)=cov;
theta=estimate(y,X,theta0,lb,opt,p,q,ncov,T);


omega=theta(1);
alpha=theta(2:p+1);
beta=theta(p+2:p+q+1);
gamma=theta(p+q+2:p+q+1+ncov);

global lambda
lambda=zeros(T+max(p,q),1);
lambda(1:max(p,q))=0;
s_omega=zeros(T+max(p,q));
s_omega(1:max(p,q))=0;
s_alpha=zeros(T+max(p,q),p);
s_alpha(1:max(p,q),:)=0;
s_beta=zeros(T+max(p,q),q);
s_beta(1:max(p,q),:)=0;
s_gamma=zeros(T+max(p,q),ncov);
s_gamma(1:max(p,q),:)=0;
l=zeros(T+max(p,q),1);


A=zeros(p+q+1+ncov,p+q+1+ncov);

%Der loopes gennem hele modellen, 
for i=max(p,q)+1:T+max(p,q)
    lambda(i)=omega+sum(alpha'.*y(i-1:-1:i-p))+sum(beta.*lambda(i-1:-1:i-q))+sum(gamma'.*abs(X(i-1,:)));
    %lambda(i)=omega+sum(alpha'.*y(i-1:-1:i-p))+sum(beta.*lambda(i-1:-1:i-q));
    l(i)=y(i)*log(lambda(i))-lambda(i);    
    
    s_omega(i)=1+sum(beta'.*s_omega(i-1:-1:i-q));
    for j=1:p
    s_alpha(i,j)=y(i-j)+sum(beta.*s_alpha(i-1:-1:i-q,j));
    end
    for h=1:q
    s_beta(i,h)=lambda(i-h)+sum(beta.*s_beta(i-1:-1:i-q,h));
    end
    for m=1:ncov
    s_gamma(i,m)=abs(X(i-1,m))+sum(beta.*s_gamma(i-1:-1:i-q,m));
    end
    
    %Hvad foregår der på den her linje? 
    G_=1/lambda(i)*[s_omega(i);s_alpha(i,:)';s_beta(i,:)';s_gamma(i,:)']*[s_omega(i);s_alpha(i,:)';s_beta(i,:)';s_gamma(i,:)']';
    G=A+G_;
    A=G;
    
   %det har har jo noget og gøre med, at det er et loop ... inverterer y(i)
   %og kalder det G_ ... summer så til G... som de så "fylder ind" på A's pladser? .. A er en predefineret nulmatrice %omdøber til A
    
    
end

logl=sum(l);
global S
S=G^-1;
%Så vi inverterer G tilbage igen? og kalder den S? noget trick med positivt
%og negativt?

OUT=zeros(p+q+1+ncov+1,3);
OUT(p+q+1+ncov+1,1)=logl;
OUT(p+q+1+ncov+1,2:3)=aicbic(logl,p+q+1+ncov,T);


for k=1:p+q+1+ncov
    OUT(k,1)=theta(k);
    OUT(k,2)=sqrt(S(k,k));
    OUT(k,3)=OUT(k,1)/OUT(k,2);
   
end

global pears
for i=max(p,q)+1:T+max(p,q)
    pears(i)=y(i)-lambda(i)/sqrt(lambda(i));

end
%temp = array2table(OUT);
%xlswrite(filename,OUT,s)
writematrix(OUT,filename,'Sheet',s)
end





