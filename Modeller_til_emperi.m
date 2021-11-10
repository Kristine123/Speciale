%% Alle model estimater til den empiriske del
clc
clear
%%
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');  
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_16.txt'); 
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_16.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_16.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
%%

y=hosp';
%X(1,:)=temp;
X(1,:)=temp;
X(2,:)=contact;
X(3,:)=tested;
X(4,:)=vacc;

%OBS y og X er her transponeret ifht den originale kode for at passe ind i
%output_parx_rob funktionen


%%
%PARX(1,1)
%Benchmark model

p = 1; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[0.2;0.1;0.2;0.2;0.54;0.01;0.81]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX_22/output_parx_22.xlsx';

%function[OUT,pears,pred]=output_parx_rob(y,X,theta0,p,q,s,filename)
[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  

%%
%PARX(2,1)
%Benchmark model

p = 2; % autoreg led
q = 1; % intensity lagged
s = 2; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX_22/output_parx_22.xlsx';

%function[OUT,pears,pred]=output_parx_rob(y,X,theta0,p,q,s,filename)
[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  

%%
%%PARX med vacc og temp
y=hosp';
%X(1,:)=temp;
X(1,:)=temp;
X(2,:)=vacc;
%X(3,:)=tested;
%X(4,:)=vacc;

%%
%PARX(1,1)
%Benchmark model

p = 1; % autoreg led
q = 1; % intensity lagged
s = 3; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[0.2;0.1;0.2;0.2;0.54]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX_22/output_parx_22.xlsx';

%function[OUT,pears,pred]=output_parx_rob(y,X,theta0,p,q,s,filename)
[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  
%%PARX med vacc og temp
y=hosp';
%X(1,:)=temp;
X(1,:)=temp;
X(2,:)=vacc;
%X(3,:)=tested;
%X(4,:)=vacc;

%%
%%PARX med temp
y=hosp';
%X(1,:)=temp;
X(1,:)=temp;
%X(2,:)=vacc;
%X(3,:)=tested;
%X(4,:)=vacc;

%%
%PARX(1,1)
%Benchmark model

p = 1; % autoreg led
q = 1; % intensity lagged
s = 4; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[0.2;0.1;0.2;0.2]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX_22/output_parx_22.xlsx';

%function[OUT,pears,pred]=output_parx_rob(y,X,theta0,p,q,s,filename)
[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  
