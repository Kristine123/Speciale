clc
%%
%PAR model with different lag lengths
%Der skal ændres til de to grønne linjer i output og estimate for at køre
%koden. 

hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
y=hosp;

ones=importdata('/Users/krmmm/Documents/Speciale/Model_data/1.txt');
X=ones;
%Denne er ligegyldig - køres ikke med - det er bare startværdien den
%tager med

%Disse varieres
p = 2; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
%k = p+q+1;

theta0=[0.1;0.2;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
%theta0 = [0.1*ones(k,1);0.0];
%Denne bid kører ikke, der er noget galt med dimensionerne
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*ones(length(theta0),1);
%lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PAR/estimates_output_par_p2.xlsx';

out1 = output(y,X,theta0,lb,opt,p,q,s,filename);   
%output(y,X,theta0,lb,opt,2,1,1,filename);
%%output(data,cov,theta0,lb,opt,p,q,s,output_filename)
%output(data,cov,theta0,lb,opt,p,q,s,output_filename)
% output = LogL, AIC, BIC  x 6
%output_parx_rob
%y,X,theta0,p,q,s,filename

%%
hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
y=hosp;

ones=importdata('/Users/krmmm/Documents/Speciale/Model_data/1.txt');
X=ones;
%Denne er ligegyldig - køres ikke med - det er bare startværdien den
%tager med

%Disse varieres
p = 7; % autoreg led
q = 1; % intensity lagged
s = 2; % Excel sheet 

%theta0=[0.1;0.2;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
theta0 = [0.1*ones(p+q+1,1);0.0];
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*ones(length(theta0),1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PAR/estimates_output_par_p7.xlsx';

out1 = output_parx_rob(y,X,theta0,p,q,s,filename);   
%output(data,cov,theta0,lb,opt,p,q,s,output_filename)
% output = LogL, AIC, BIC  x 6

%%



%Model with one covariate, p=2, q=1:


test=importdata('/Users/krmmm/Documents/Speciale/Model_data/Tested_366.txt');
X=test;


theta0=[0.1;0.2;0.2;0.2;0.2];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PARX_test_contact/estimates_output_one_covariate_test.xlsx';

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