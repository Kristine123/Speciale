%Model with one covariate, p=1, q=2:

def=importdata('/Users/krmmm/Documents/Speciale/PARX_codes/PARX_code_October2014 2/code&script/defcountUS.txt');
y=def;

rv=importdata('/Users/krmmm/Documents/Speciale/PARX_codes/PARX_code_October2014 2/code&script/rv.txt');
X=rv;

theta0=[0.1;0.2;0.2;0.2;0.2];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');

output(y,X,theta0,lb,opt,2,1,1)   %%the final 1 means it writes in the first sheet of the excel file%%
%output(data,cov,theta0,lb,opt,p,q,s)
% output = LogL, AIC, BIC  x 6

%% Model with two covariates:

def = importdata('/Users/krmmm/Documents/Speciale/PARX_codes/PARX_code_October2014 2/code&script/defcountUS.txt');
y   = def;
rv  = importdata('/Users/krmmm/Documents/Speciale/PARX_codes/PARX_code_October2014 2/code&script/rv.txt');
rvlag =importdata('/Users/krmmm/Documents/Speciale/rvlag2.txt');
spread=importdata('/Users/krmmm/Documents/Speciale/PARX_codes/PARX_code_October2014 2/code&script/credit_spread.txt');
X(:,1)=rv;
X(:,2)=rvlag;
theta0=[0.1;0.2;0.2;0.2;0.2;0.2];
lb=-Inf*ones(6,1);

output(y,X,theta0,lb,opt,2,1,1)

%%




%Model when you use the dummy for taking only negative values:

def=importdata('defcountUS.txt');
y=def;
ip=importdata('IP_index_change.txt');
ip=importdata('ipindex.txt');

for i=1:length(ip)
    
    if ip(i)>=0
        ip_(i)=0;    
    else ip_(i)=ip(i);    
    end
    
end

X(:,1)=ip_';
theta0=[0.1;0.2;0.2;0.2;0.2];     
lb=-Inf*ones(5,1);
output(y,X,theta0,lb,opt,2,1,1)



