%Forsøg på at kører output_parx_rob
clc
clear
%%

hosp=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/newlyadmitted_366.txt');
%temp=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/Temperatur_366.txt');
%X=temp;
fahrenheit = importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/fahrenheit_366.txt'); 
contact=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/kontakttal_366.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/Tested_366.txt');
vaccination=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/vaccinationer_pctafbefolkning_366.txt');

%%

y=hosp';
%X(1,:)=temp;
X(1,:)=fahrenheit;
X(2,:)=contact;
X(3,:)=tested;
X(4,:)=vaccination;

%OBS y og X er her transponeret ifht den originale kode for at passe ind i
%output_parx_rob funktionen


%%
%PARX(14,1)
%Benchmark model

p = 1; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[1.4;0.18;0.11;0.08;0.14;0.02;0.01]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX/output_parx_rob_05.xlsx';

out1 = output_parx_rob(y,X,theta0,p,q,s,filename);  
%function[OUT]=output_parx_rob(y,X,theta0,pred,pears,p,q,s,filename)

%%
%PARX(14,2)
%Benchmark model

p = 14; % autoreg led
q = 2; % intensity lagged
s = 2; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[1.4;0.18;0.11;0.08;0.14;0.02;0.01;0.05;0.00;0.05;0.04;0.06;0.23;0.3;0.07;0.01;0.01;0.5;1.9;1.4;0.5]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX/output_parx_rob_juli.xlsx';

out1 = output_parx_rob(y,X,theta0,p,q,s,filename);  

%%
%PARX(7,1)
%Benchmark model

p = 7; % autoreg led
q = 1; % intensity lagged
s = 3; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[1.4;0.18;0.11;0.08;0.14;0.02;0.01;0.05;0.60;0.05;0.04;0.06;0.23;0.3;0.07;0.01;0.01;0.5;1.9;1.4]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX/output_parx_rob_juli.xlsx';

out1 = output_parx_rob(y,X,theta0,p,q,s,filename);  
