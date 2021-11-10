function[theta]=estimate(y,X,theta0,lb,opt,p,q,ncov,T)
theta=fmincon(@nested,theta0,[],[],[],[],lb,[],[],opt);


function z=nested(theta)

omega=theta(1);
alpha=theta(2:p+1);
beta=theta(p+2:p+q+1);
gamma=theta(p+q+2:p+q+1+ncov);


global lambda
lambda=zeros(T+max(p,q),1);
lambda(1:max(p,q))=0;
l=zeros(T+max(p,q),1);

for i=max(p,q)+1:T+max(p,q)
    lambda(i)=omega+sum(alpha'.*y(i-1:-1:i-p))+sum(beta.*lambda(i-1:-1:i-q))+sum(gamma'.*abs(X(i-1,:)));
    %lambda(i)=omega+sum(alpha'.*y(i-1:-1:i-p))+sum(beta.*lambda(i-1:-1:i-q));
    l(i)=y(i)*log(lambda(i))-lambda(i);
end

z=-sum(l);
end
end