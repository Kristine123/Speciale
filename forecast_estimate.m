function[theta_hat]=forecast_estimate(y,X,theta0,lb,opt,p,q,ncov,T)
%theta_hat=fmincon(@nested,theta0,[],[],[],[],lb,[],[],opt);
beta_ = [];
lambda = [];

    function z=nested(theta)

    omega=theta(1);
    alpha=theta(2:p+1);
    beta_=theta(p+2:max(p+2,p+q+1));
    gamma=theta(p+q+2:p+q+1+ncov);

    if q<=0
      beta_=0;
    end

    lambda=zeros(T+max(p,q),1);
    lambda(1:max(p,q))=0;
    l=zeros(T+max(p,q),1);

    for i=max(p,q)+1:T+max(p,q)
        lambda(i)=omega+sum(alpha'.*y(i-1:-1:i-p))+sum(beta_.*lambda(i-1:-1:i-q))+sum(gamma'.*(X(i-1,:)));
        %lambda(i)=omega+sum(alpha'.*y(i-1:-1:i-p))+sum(beta_.*lambda(i-1:-1:i-q));
        l(i)=y(i)*log(lambda(i))-lambda(i);
    end

    z=-sum(l);
    end

    theta=fmincon(@nested,theta0,[],[],[],[],lb,[],[],opt);
    %options = optimoptions(@fmincon,'MaxFunctionEvaluations',10000)


    SE=zeros(p+q+1+ncov,1);
    s_omega=zeros(T+max(p,q));
    s_omega(1:max(p,q))=0;
    s_alpha=zeros(T+max(p,q),p);
    s_alpha(1:max(p,q),:)=0;
    s_beta=zeros(T+max(p,q),q);
    s_beta(1:max(p,q),:)=0;
    s_gamma=zeros(T+max(p,q),ncov);
    s_gamma(1:max(p,q),:)=0;

    A=zeros(p+q+1+ncov,p+q+1+ncov);

    for m=max(p,q)+1:T+max(p,q)

        s_omega(m)=1+sum(beta_'.*s_omega(m-1:-1:m-q));

        for j=1:p
        s_alpha(m,j)=y(m-j)+sum(beta_.*s_alpha(m-1:-1:m-q,j));
        end

        for h=1:q
        s_beta(m,h)=lambda(m-h)+sum(beta_.*s_beta(m-1:-1:m-q,h));
        end

        for o=1:ncov
        s_gamma(m,o)=X(m-1,o)+sum(beta_.*s_gamma(m-1:-1:m-q,o));
        end


        G_=1/lambda(m)*[s_omega(m);s_alpha(m,:)';s_beta(m,:)';s_gamma(m,:)']*[s_omega(m);s_alpha(m,:)';s_beta(m,:)';s_gamma(m,:)']';
        G=A+G_;
        A=G;

    end

    SE=sqrt(diag(G^-1));
    theta_hat=zeros(p+q+1+ncov,1);


    for n=1:p+q+1+ncov
        if abs(theta(n)/SE(n))>1.65
        theta_hat(n)=theta(n);
        else
        theta_hat(n)=0;
        end
    end
end








   
    
    
    



