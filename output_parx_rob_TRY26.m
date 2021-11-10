%Forsøg på at kører output_parx_rob
clc
clear
%%

hosp=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/newlyadmitted_366.txt');
temp=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/Temperatur_366.txt');
%X=temp;
contact=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/kontakttal_366.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/Tested_366.txt');
vaccination=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/vaccinationer_pctafbefolkning_366.txt');

%%

y=hosp';
X(1,:)=temp;
X(2,:)=contact;
X(3,:)=tested;
X(4,:)=vaccination;

%OBS y og X er her transponeret ifht den originale kode for at passe ind i
%output_parx_rob funktionen


%%
%PARX(14,1)
%Benchmark model

%Disse varieres måske senere 
p = 14; % autoreg led
q = 1; % intensity lagged
s = 6; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[1.4;0.18;0.11;0.08;0.14;0.02;0.01;0.05;0.00;0.05;0.04;0.06;0.23;0.3;0.07;0.01;0.01;0.5;1.9;1.4]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX/output_parx_rob_juni.xlsx';

% [OUT]=output_parx_rob(y,X,theta0,lb,opt,p,q,s,filename)
out1 = output_parx_rob(y,X,theta0,p,q,s,filename);  

% output_parx_rob(y,X,theta0,lb,opt,p,q,s,filename) 
% output(y,X,theta0,lb,opt,2,1,1,filename);
% output = LogL, AIC, BIC  x 6

%%
%Problem: Vi får forskellige resultater med den nye og den gamle kode: det
%er nok fordi der ikke er nok itterationer på den gamle: letteste løsning
%er at kunne sætte disse op, næste er manuelt at justere theta værdierne.
%Du har lige konstateret det er det der er problemet, så det behøver du
%ikek konstatere en gang til.

