%Forecast per 26-05-2021

%Der skal ændres til de to grønne linjer i output og estimate for at køre
%koden for PAR

%%
clc
clear

hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
temp=importdata('/Users/krmmm/Documents/Speciale/Model_data/Temperatur_366.txt');
contact=importdata('/Users/krmmm/Documents/Speciale/Model_data/kontakttal_366.txt');
tested=importdata('/Users/krmmm/Documents/Speciale/Model_data/Tested_366.txt');
vaccination=importdata('/Users/krmmm/Documents/Speciale/Model_data/vaccinationer_pctafbefolkning_366.txt');

%Bruges kun til PAR
%ones=importdata('/Users/krmmm/Documents/Speciale/Model_data/1.txt');
%X=ones;

%%

y=hosp;
X(:,1)=temp;
X(:,2)=contact;
X(:,3)=tested;
X(:,4)=vaccination;

%X=temp; , for kun en variabel
%%

%Disse varieres - længden af theta0 justeres pt manuelt
p = 14; % autoreg led
q = 3; % intensity lagged
s = 2; % Excel sheet 
k = 7;
%o = p+q+1;

theta0=[0.2;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.1;0.2;0.2;0.2;0.1;0.2]; 
% Denne kan automatiseres
lb=-Inf*ones(length(theta0),1);
opt=optimset('Algorithm','interior-point');
 
filename = '/Users/krmmm/Documents/MATLAB/PARX_1/Outputs/Forecasts/Forecasts_26.xlsx';

out1 = forecast_output(y,X,theta0,lb,opt,p,q,k,s,filename); %%the final 1 means it writes in the first sheet of the excel file%%
%forecast_output(data,cov,theta0,lb,opt,p,q,k,s,filename)
% output = LogL, AIC, BIC  x 6

%%

%Disse varieres - længden af theta0 justeres pt manuelt
p = 14; % autoreg led
q = 3; % intensity lagged
s = 2; % Excel sheet 
k = 21;
%o = p+q+1;

theta0=[0.2;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.1;0.2;0.2;0.2;0.1;0.2]; 
% Denne kan automatiseres
lb=-Inf*ones(length(theta0),1);
opt=optimset('Algorithm','interior-point');
 
filename = '/Users/krmmm/Documents/MATLAB/PARX_1/Outputs/Forecasts/Forecasts_26.xlsx';

out1 = forecast_output(y,X,theta0,lb,opt,p,q,k,s,filename); %%the final 1 means it writes in the first sheet of the excel file%%
%forecast_output(data,cov,theta0,lb,opt,p,q,k,s,filename)
% output = LogL, AIC, BIC  x 6