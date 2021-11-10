%RAR & PARX(2,1) med covariater antal testede og kontakttal 

%PAR p=2, q=1:
%Der skal ændres til de to grønne linjer i output og estimate for at køre
%koden. 
clc
          
hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
y=hosp;


ones=importdata('/Users/krmmm/Documents/Speciale/Model_data/1.txt');
X=ones;


theta0=[0.1;0.2;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
lb=-Inf*ones(length(theta0),1);
opt=optimset('Algorithm','interior-point');
 
filename = '/Users/krmmm/Documents/MATLAB/PARX_1/Outputs/PAR/estimates_output_TRY_forecast3.xlsx';

out1 = forecast_output(y,X,theta0,lb,opt,2,1,2,1,filename); %%the final 1 means it writes in the first sheet of the excel file%%
%function[OUT]=output(data,cov,theta0,lb,opt,p,q,s,filename)
%function[OUT]=output(data,cov,theta0,lb,opt,p,q,s,filename)
%output(data,cov,theta0,lb,opt,p,q,s,output_filename)
% output = LogL, AIC, BIC  x 6