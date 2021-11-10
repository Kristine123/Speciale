
clc
clear
%%
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');

%%
y=hosp';
%X(1,:)=temp;
X= [];
%%
%PAR(1,1)
%Benchmark model

p = 1; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
%g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[0.2;0.1;0.2;0.2]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PAR_22/output_par_22.xlsx';

%function[OUT,pears,pred]=output_parx_rob(y,X,theta0,p,q,s,filename)
[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  

%%
figure; autocorr(pearson_residuals)

%%
%PAR(2,1)
%Benchmark model

p = 2; % autoreg led
q = 1; % intensity lagged
s = 2; % Excel sheet 
%g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[0.2;0.1;0.2;0.2;0.2;0.2]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PAR_22/output_par_22.xlsx';

%function[OUT,pears,pred]=output_parx_rob(y,X,theta0,p,q,s,filename)
[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  
 
