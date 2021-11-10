% Estimationer til vejledermøde d. 25-05-2021

%Par modeller
%lag 2, 7, 14, 21

%Der skal ændres til de to grønne linjer i output og estimate for at køre
%koden. 

hosp=importdata('/Users/krmmm/Documents/Speciale/Model_data/newlyadmitted_366.txt');
y=hosp;

ones=importdata('/Users/krmmm/Documents/Speciale/Model_data/1.txt');
X=ones;
%Denne er ligegyldig - køres ikke med - det er bare startværdien den
%tager med

%%

%Par(2,1)

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
% output = LogL, AIC, BIC  x 6


%%

%Par(7,1)

%Disse varieres
p = 7; % autoreg led
q = 1; % intensity lagged
s = 2; % Excel sheet 
%k = p+q+1;

theta0=[0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
%theta0 = [0.1*ones(k,1);0.0];
%Denne bid kører ikke, der er noget galt med dimensionerne
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*ones(length(theta0),1);
%lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PAR/estimates_output_par_p2.xlsx';

out2 = output(y,X,theta0,lb,opt,p,q,s,filename);   
%output(y,X,theta0,lb,opt,2,1,1,filename);
% output = LogL, AIC, BIC  x 6


%%

%Par(14,1)

%Disse varieres
p = 14; % autoreg led
q = 1; % intensity lagged
s = 3; % Excel sheet 
%k = p+q+1;

theta0=[0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.2;0.2;0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
%theta0 = [0.1*ones(k,1);0.0];
%Denne bid kører ikke, der er noget galt med dimensionerne
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*ones(length(theta0),1);
%lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PAR/estimates_output_par_p2.xlsx';

out2 = output(y,X,theta0,lb,opt,p,q,s,filename);   
%output(y,X,theta0,lb,opt,2,1,1,filename);
% output = LogL, AIC, BIC  x 6


%%

%Par(21,1)

%Disse varieres
p = 21; % autoreg led
q = 1; % intensity lagged
s = 4; % Excel sheet 
%k = p+q+1;

theta0=[0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.2;0.2;0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.1;0.2;0.2;0.2;0.1;0.2;0.2;0.0];     %%I USE ARBITRARY INITIAL VALUES AS I CHECKED THE ESTIMATES ARE NOT SENSITIVE %%
%theta0 = [0.1*ones(k,1);0.0];
%Denne bid kører ikke, der er noget galt med dimensionerne
%Her er 0.0 tilføjet fordi det er en gypsy løsning PAR model
lb=-Inf*ones(length(theta0),1);
%lb=-Inf*ones(5,1);
opt=optimset('Algorithm','interior-point');


filename = 'Outputs/PAR/estimates_output_par_p2.xlsx';

out2 = output(y,X,theta0,lb,opt,p,q,s,filename);   
%output(y,X,theta0,lb,opt,2,1,1,filename);
% output = LogL, AIC, BIC  x 6
% PARX modeller
%lag 2, 7, 14, 21 alle covariater
