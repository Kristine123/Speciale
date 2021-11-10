% Estimationer til vejledermøde d. 25-05-2021

%Parx modeller

%Der skal ændres til de to grønne linjer i output og estimate for at køre
%koden. 

%%
clc
clear

hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
temp=importdata('/Users/krmmm/Documents/Speciale/Model_data/Temperatur_366.txt');
%X=temp;
contact=importdata('/Users/krmmm/Documents/Speciale/Model_data/kontakttal_366.txt');
tested=importdata('/Users/krmmm/Documents/Speciale/Model_data/Tested_366.txt');
vaccination=importdata('/Users/krmmm/Documents/Speciale/Model_data/vaccinationer_pctafbefolkning_366.txt');

%%

y=hosp;
X(:,1)=temp;
X(:,2)=contact;
X(:,3)=tested;
X(:,4)=vaccination;

%%

%Disse varieres måske senere 
p = 2; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
%k = p+q+1;

theta0=[0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.1];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
%theta0 = [0.1*ones(k,1);0.0];
%Denne bid kører ikke, der er noget galt med dimensionerne
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*temp(length(theta0),1);
%lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/estimates_output_parx_25.xlsx';

out1 = output(y,X,theta0,lb,opt,p,q,s,filename);   
%output(y,X,theta0,lb,opt,2,1,1,filename);
% output = LogL, AIC, BIC  x 6


%%
%PARX(7,1)

%Disse varieres måske senere 
p = 7; % autoreg led
q = 1; % intensity lagged
s = 2; % Excel sheet 
%k = p+q+1;

theta0=[0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
%theta0 = [0.1*ones(k,1);0.0];
%Denne bid kører ikke, der er noget galt med dimensionerne
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*temp(length(theta0),1);
%lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/estimates_output_parx_25.xlsx';

out1 = output(y,X,theta0,lb,opt,p,q,s,filename);   
%output(y,X,theta0,lb,opt,2,1,1,filename);
% output = LogL, AIC, BIC  x 6

%%
%PARX(14,1)

%Disse varieres måske senere 
p = 14; % autoreg led
q = 1; % intensity lagged
s = 3; % Excel sheet 
%k = p+q+1;

theta0=[0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1;0.2;0.2];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
%theta0 = [0.1*ones(k,1);0.0];
%Denne bid kører ikke, der er noget galt med dimensionerne
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*temp(length(theta0),1);
%lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/estimates_output_parx_25.xlsx';

out1 = output(y,X,theta0,lb,opt,p,q,s,filename);   
%output(y,X,theta0,lb,opt,2,1,1,filename);
% output = LogL, AIC, BIC  x 6

%%
%PARX(21,1)

%Disse varieres måske senere 
p = 21; % autoreg led
q = 1; % intensity lagged
s = 4; % Excel sheet 
%k = p+q+1;

theta0=[0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1;0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.1;0.2;0.2;0.2;0.1;0.2;0.2];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
%theta0 = [0.1*ones(k,1);0.0];
%Denne bid kører ikke, der er noget galt med dimensionerne
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*temp(length(theta0),1);
%lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/estimates_output_parx_25.xlsx';

out1 = output(y,X,theta0,lb,opt,p,q,s,filename);   
%output(y,X,theta0,lb,opt,2,1,1,filename);
% output = LogL, AIC, BIC  x 6