%Forsøg på at køre modellen med eget data

%Model with one covariate, p=1, q=2:

fed=importdata('/Users/krmmm/Documents/Speciale/Model_data_hospitalization.txt');
y=fed;

con=importdata('/Users/krmmm/Documents/Speciale/contact_numbers.txt');
X=con;

%%

theta0=[0.1;0.2;0.2;0.2;0.2];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');

%%

filename = 'estimates_output_test2.xlsx';

out1 = output(y,X,theta0,lb,opt,2,1,1,filename);   %%the final 1 means it writes in the first sheet of the excel file%%
%output(data,cov,theta0,lb,opt,p,q,s)
% output = LogL, AIC, BIC  x 6

%%

%Model with two covariates:

def=importdata('/Users/krmmm/Documents/Speciale/hospitalization.txt');
y=def;
rv=importdata('/Users/krmmm/Documents/Speciale/contact_numbers.txt');
spread=importdata('/Users/krmmm/Documents/Speciale/Tested_dk_382.txt');
X(:,1)=rv;
X(:,2)=spread;
theta0=[0.1;0.2;0.2;0.2;0.2;0.2];
lb=-Inf*ones(6,1);
%output(y,X,theta0,lb,opt,2,1,1)

filename = 'estimates_output_test3.xlsx';

out2 = output(y,X,theta0,lb,opt,2,1,1,filename);