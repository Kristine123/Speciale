clc

%attempts with different p and q
%RAR & PARX(2,1) med covariater antal testede og kontakttal 

%PAR p=2, q=1:
%Der skal ændres til de to grønne linjer i output og estimate for at køre
%koden. 

hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
y=hosp;


ones=importdata('/Users/krmmm/Documents/Speciale/Model_data/1.txt');
X=ones;

theta0=[0.1;0.2;0.2;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
lb=-Inf*ones(length(theta0),1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PAR/estimates_output_par_q2.xlsx';

out1 = output_parx_rob(y,X,theta0,lb,opt,2,2,1,filename);   %%the final 1 means it writes in the first sheet of the excel file%%
%output(data,cov,theta0,lb,opt,p,q,s,output_filename)
% output = LogL, AIC, BIC  x 6

%%

%attempts with different p and q
%RAR & PARX(2,1) med covariater antal testede og kontakttal 

%PAR p=2, q=1:
%Der skal ændres til de to grønne linjer i output og estimate for at køre
%koden. 

hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
y=hosp;


ones=importdata('/Users/krmmm/Documents/Speciale/Model_data/1.txt');
X=ones;


%%
theta0=[0.1;0.2;0.2;0.2;0.2;0.2;0.1;0.2;0.2;0.2;0.2;0.2;0.2;0.1;0.2;0.2;0.2;0.2;0.2;0.1;0.2;0.2;0.2;0.2;0.2;0.2;;0.2;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
lb=-Inf*ones(length(theta0),1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PAR/estimates_output_par_q21p21.xlsx';

out1 = output(y,X,theta0,lb,opt,21,21,filename);   %%the final 1 means it writes in the first sheet of the excel file%%
%output(data,cov,theta0,lb,opt,p,q,s,output_filename)
% output = LogL, AIC, BIC  x 6



%%

%attempts with different p and q
%RAR & PARX(2,1) med covariater antal testede og kontakttal 

%PAR p=2, q=1:
%Der skal ændres til de to grønne linjer i output og estimate for at køre
%koden. 

hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
y=hosp;


ones=importdata('/Users/krmmm/Documents/Speciale/Model_data/1.txt');
X=ones;

theta0=[0.1;0.2;0.2;0.2;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
lb=-Inf*ones(length(theta0),1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PAR/estimates_output_par_q2p4.xlsx';

out1 = output(y,X,theta0,lb,opt,3,2,1,filename);   %%the final 1 means it writes in the first sheet of the excel file%%
%output(data,cov,theta0,lb,opt,p,q,s,output_filename)
% output = LogL, AIC, BIC  x 6
%%


%Model with one covariate, p=2, q=1:


test=importdata('/Users/krmmm/Documents/Speciale/Model_data/Tested_366.txt');
X=test;


theta0=[0.2;0.2;0.2;0.2;0.2];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PARX_test_contact/estimates_output_one_covariate_test_theta.xlsx';

out2 = output(y,X,theta0,lb,opt,2,1,1,filename);   %%the final 1 means it writes in the first sheet of the excel file%%
%output(data,cov,theta0,lb,opt,p,q,s)
% output = LogL, AIC, BIC  x 6

%%

%Model with two covariates:

%hosp=importdata('/Users/krmmm/Documents/Speciale/hospitalization.txt');
%y=hosp;

%test=importdata('/Users/krmmm/Documents/Speciale/contact_numbers.txt');

contact=importdata('/Users/krmmm/Documents/Speciale/Model_data/kontakttal_366.txt');

X(:,1)=test;
X(:,2)=contact;
theta0=[0.1;0.2;0.2;0.2;0.2;0.2];
lb=-Inf*ones(6,1);
%output(y,X,theta0,lb,opt,2,1,1)

filename = 'Outputs/PARX_test_contact/estimates_output_two_covariates_test_contact.xlsx';

out3 = output(y,X,theta0,lb,opt,2,1,1,filename);