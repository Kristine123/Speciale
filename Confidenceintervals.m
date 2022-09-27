%% Experimenting with confidence intervals

%% Estimate PARX
clc 
clear all
close all

%% Import data
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');     
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_18.txt'); % Ikke transformerede. Trans er i temp_16
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_18.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');

%%
%Covariates
X(1,:)=max(0,diff(tested));%Only significant variable in empirical analysis

%Length of time series
 startDate = 1;
 endDate = 464;
 
%Length of timeseries continued
y = hosp';
y = y(:,1:endDate); 

   y = y(:,startDate:endDate); %To shorten time series further use this format 
   X = X(:,startDate:endDate);

% PARX(1,1)

%Benchmark model

p = 1; % autoreg lags
q = 1; % intensity lagged
s = 1; % Excel sheet 
[r,T] = size(X);
g = p+q+r+1; % determine number of start values used

% Set initial values
% theta0 = [0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81];
theta0 = [0.5;0.3;0.1;0.26;0.64;0.05;0.71;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81];
theta0 = theta0(1:g);

% Set output save location
filename = 'Outputs/PARX_22/0112_rerun_the_results.xlsx';

% Use parx_rob function to get outputs

[OUT, pearson_residuals, pred, probability, meany, rpit,confidence_interval] = output_parx_rob(y(:,1:end),X(:,1:end),theta0,p,q,s,filename); 

% If a PAR is estimated: output_parx_rob(y(:,1:T0),X,theta0,p,q,s,filename);
% If a PARX is estimated: output_parx_rob(y(:,1:T0),X(:,1:T0),theta0,p,q,s,filename);

%% Calculating upper confindence interval lambda IKKE = Var
%start at 2 because need at least 1 degree of freedom to calculate
%confidence level:
p = .95;
%confidence_interval_upper = zeros((size(pred)-5));

for i = 2:size(pred)
%at leat one degree of fredoom needed to calculate thats why we start on 2
nu = i-1 ; 
t = tinv(p,nu);
%Her har du antaget lambda IKKE lig mean
upper = pred(i) + t * (std(pred(1:i))/sqrt(size(pred));  
%To adjust 
confidence_interval_upper(i-1) = upper;
end
% The result is a way to narrow interval

%% Calculating upper confindence interval lambda IKKE = Var OG der er ikke dobbelt normaliseret 
%start at 2 because need at least 1 degree of freedom to calculate
%confidence level:
p = .95;
%confidence_interval_upper = zeros((size(pred)-5));

for i = 2:size(pred)
%at leat one degree of fredoom needed to calculate thats why we start on 2
nu = i-1 ; 
t = tinv(p,nu);
%Her har du antaget lambda IKKE lig mean
upper = pred(i) + t * std(pred(1:i);  
%To adjust 
confidence_interval_upper(i-1) = upper;
end
% The result is a way to narrow interval

%% Calculating upper confindence interval - t test where lambda = var

for i = 2:size(pred)
nu = i-1 ; 
t = tinv(p,nu);
%Her har du antaget lambda = mean = var
upper = pred(i) + t * (sqrt(pred(i))/sqrt(i));  
confidence_interval_upper1(i) = upper;
end
% The result is a way to narrow interval - the adjustment makes no
% difference

%% Attempt with chi^2 

for i = 2:size(pred)
%at leat one degree of fredoom needed to calculate thats why we start on 2
nu = i-1 ; 
t = chi2inv(p,nu);
upper = (nu*var(pred(1:i)))/t;
%pred(i) + t * (pred(i)/sqrt(i));  
%To adjust 
confidence_interval_upper(i-1) = upper;
end


%% PLOT it Predicted defaults per month
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2-3);
%Date = (t1:t2-2);

figure; 
u = y(p+1:end-1)
%u = y(p+1:end)
%plot(Date,u,Date,pred,Date, confidence_interval, 'LineWidth',1.5)
plot(Date,u,Date,pred(1:end-1),Date, confidence_interval_upper, 'LineWidth',1.5)
%plot(Date,g,'LineWidth',0.7)
ylabel('Daily Hospitalizations','FontSize',  14)
xlabel('Date','FontSize',  14)
%plot(Date,confidence_interval)

xlim(datetime([2020 2021],[03 08],1));
%title('PAR(1,1) Estimates  vs Actual Observations','FontSize', 14);
legend('Actual Hospitalizations','Estimated Hospitalizations','Confidence interval lower bound','FontSize', 14);
ax = gca;
ax.FontSize = 14;
grid minor

saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/TestedP/predictions2.jpg')


%pred = [1:100];
