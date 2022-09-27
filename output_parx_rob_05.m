%% PARX model estimation
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

%% Select data to use

% X=[]; %to run a PAR model activate this line

% Transformed variables for PARX model

% test = ["a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11"];
% for c = 1:11
%     test(1,c);
% end
% 
% a3 = 1;
% a4 = 2;

% X(a1,:)=abs(min(0,diff(temp)));
% X(a2,:)=max(0,diff(temp));
% 
% X(1,:)=abs(min(0,diff(contact)));
% X(2,:)=max(0,diff(contact));
% 
% X(a5,:)=abs(min(0,diff(tested))); 
 X(1,:)=max(0,diff(tested));%Only significant variable in empirical analysis
% 
% X(a7,:)=abs(diff(vacc)); 
% 
% X(a8,:)=abs(min(0,diff(varbrit)));
% X(a9,:)=max(0,diff(varbrit));
% 
% X(a10,:)=abs(min(0,diff(vardelta)));
% X(a11,:)=max(0,diff(vardelta));
 %% Length of time series
 
 startDate = 1;
 endDate = 464;
% endDate observations % observations on 55, (9/10-20)197, (15/02-21)326, (02/07-21)464
% endDate is set to lenght of time series -1 when including the exogenes variables. This is necessary as taking first difference of exogenes variables results
% in a shorter time series
 %% Dummy
% Creates dummy for robustness analysis
% dummy = zeros(endDate,1);
% 
% if (endDate == 464)
%     for i = 1:endDate
%         if(i >= 197)
%         dummy(i) = 1;
%         end
%     end
% else
%     for i = 1:endDate
%         if(i >= 197)
%         dummy(i) = 1;
%         end
%     end
%     for i = endDate:464
%         dummy(i) = 0;
%     end
% end
% 
% X(1,:)=dummy; % To include dummy activate this line

%% Length of timeseries continued
y = hosp';
y = y(:,1:endDate); 

   y = y(:,startDate:endDate); %To shorten time series further use this format 
   X = X(:,startDate:endDate);
%% PARX(1,1)

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

%%
% writematrix(pred, '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Outputs/PARX/pred.xlsx')
% writematrix(, '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Outputs/PARX/lambda_f_2para.xlsx')
%% Settings dates for misspecification tests

% Set dates for plotting
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
%k = p+1; %(Not first difference)
k = p+2; % (First difference)
%g = y(k:end)'; % (Not first difference)
g = y(k-1:end)'; % (First difference)
Date = Date(k:end);

%% Plot predicted values vs. actual observations

figure; 
plot(Date,pred, 'LineWidth',1.5)
hold on
plot(Date,g, 'LineWidth',0.7)
ylabel('Hospitalizations','FontSize',  14)
xlabel('Date','FontSize',  14)
hold off

xlim(datetime([2020 2021],[03 08],1));
%title('PAR(1,1) Estimates  vs Actual Observations','FontSize', 14);
legend('Estimated Hospitalizations','Actual Hospitalizations', 'FontSize', 14);
ax.FontSize = 14;
%saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/TestedP/predictions.jpg')

%% Pearson Residuals

% k = p+1;
% g = y(k:465);

figure; 
    plot(Date, pearson_residuals)
    xlim(datetime([2020 2021],[03 08],1));
    %title('Pearson Residuals', 'FontSize', 20);
    legend('Predictions','True', 'FontSize', 14);
    xlabel('Date', 'FontSize', 14);     
    ylabel('Pearson Residuals', 'FontSize', 14);
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/predictions.jpg')

%% INDSAT 02-09-2022 - confidence interval

%% Calculating upper confindence interval
%start at 2 because need at least 1 degree of freedom to calculate
%confidence level:
p = .95;
%pre-allocate ? two observations = one degree of freedom needed to
%calculate confidence interval
%confidence_interval_upper = zeros((size(pred)-5));

for i = 2:size(pred)
%at leat one degree of fredoom needed to calculate thats why we start on 2
nu = i-1 ; 
t = tinv(p,nu);
%Her har du antaget lambda = mean = var
upper = pred(i) + t * (std(pred(1:i))/sqrt(i));  
%To adjust 
confidence_interval_upper(i-1) = upper;
end

% The result is a way to narrow interval

%% Attempt with chi^2 

for i = 2:size(pred)
%at leat one degree of fredoom needed to calculate thats why we start on 2
nu = i-1 ; 
t = chi2inv(p,nu)
upper = pred(i) + t * (pred(i)/sqrt(i));  
%To adjust 
confidence_interval_upper(i-1) = upper;
end



%% Predicted defaults per month
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
%% Residual Autocorrelation of pearson residuals

figure; 
    autocorr(pearson_residuals)
    title('', 'FontSize', 20);
    legend('Predictions','True', 'FontSize', 14);
    xlabel('Lag', 'FontSize', 14);     
    ylabel('Sample Autocorrelation', 'FontSize', 14);
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/TestedP/ACF.jpg')
%% Ljungbox test

%e
[h,pValue] = lbqtest(pearson_residuals,'lags',[4,7,14])
%e^2 - ARCH effects
[h,pValue] = lbqtest(pearson_residuals.^2,'lags',[4,7,14])

%% Zero counts and probabilities

% Find supsamples of observations equal to three
hospp = hosp(1:464);
z = (hospp==3);

% Plot observations equal to three and corresponding probabilities assigned
% by the model



fig1 = figure;
    yyaxis left
    plot(Date, probability)
    ylabel('Estimated Probability', 'Color', 'k', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    
    yyaxis right
    plot(Date(z),hospp(z),'linestyle','none','marker','*','LineWidth',0.7)
    ylim([3,10])
    
    ax = gca;
    ax.FontSize = 14; 
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    set(gca,'ytick',[])
    %title('Observations equal to the mode and corresponding estimated probabilities', 'FontSize', 14)
    legend('Estimated probability of mode','Observations equal to the mode', 'FontSize', 14)
    grid minor
    
    %set(fig1,'Position',[100 100 1000 500])
    
    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/TestedP/zeroCount.jpg')
%% PIT

% Number of bins in histogram
nbins=25;

figure;
    histogram(rpit,nbins)
    title('', 'FontSize', 14)
    %xlabel('Her', 'FontSize', 14)
    %ylabel('Her også', 'FontSize', 14)
    
    ax = gca;
    ax.FontSize = 14; 
    grid minor

% Kolmogorov-Sminov test
[h,p] = kstest(norminv(rpit))

    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/TestedP/PIT3.jpg')
%% 2208testafbetydningenaftransformeretellerikketransformerettest
kk=abs(min(0,diff(tested))); %Bliver signifikant
[h,pValue] = adftest(kk)
%%
testednonnegative=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
%%
A5=abs(min(0,diff(testednonnegative)));
B6=max(0,diff(testednonnegative));
%% Initial plot
figure
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
%x = Date;
y = B6;
plot(y)
%xlim(datetime([2020 2021],[03 08],1));
title('Initial plot');
%%
kk=B6; %Bliver signifikant
[h,pValue] = adftest(kk)

%% Forsøg på at lave diagram over skewness

% Number of bins in histogram
nbins=50;

figure;
    histogram(hosp,nbins)
    title('', 'FontSize', 14)
    xlabel('Number of hospitalizations', 'FontSize', 14)
    ylabel('Counts', 'FontSize', 14)
    
    ax = gca;
    ax.FontSize = 14; 
    grid minor

% Kolmogorov-Sminov test
[h,p] = kstest(norminv(hosp))

saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/2508_endelige_misspecifikation/screwnessHistogram.jpg')
%%
x = 0:15;
y = poisscdf(x, 4.946845228258169);

figure
stairs(x,y)
xlabel('Observation')
ylabel('Cumulative Probability')

%% Attempt to calculate confidence intervals useing t-distribution : should maybe be Z?
clear;
clc;
pred = [1:100];
%goal is to do a loop through i and calculate upper an lower bound in each
%step

%zeros(size(pred)-1);

%% Calculating upper confindence interval - t test where lambda = var
%start at 2 because need at least 1 degree of freedom to calculate
%confidence level:
p = .95;
%pre-allocate ? two observations = one degree of freedom needed to
%calculate confidence interval
%confidence_interval_upper = zeros((size(pred)-1));

for i = 2:size(pred)
nu = i-1 ; 
t = tinv(p,nu);
upper = pred(i) + t * (sqrt(pred(i))/sqrt(i));  
confidence_interval_upper1(i) = upper;
end

%% Calculating upper confindence interval - t test where lambda is not equal to var
%start at 2 because need at least 1 degree of freedom to calculate
%confidence level:
p = .95;
%pre-allocate ? two observations = one degree of freedom needed to
%calculate confidence interval
%confidence_interval_upper = zeros((size(pred)-1));

for i = 2:size(pred)
nu = i-1 ; 
t = tinv(p,nu);
upper = pred(i) + t * (std(pred(1:i))/sqrt(i));  
confidence_interval_upper2(i) = upper;
end


%%
%% Calculating upper confindence interval - CHI^2
%start at 2 because need at least 1 degree of freedom to calculate
%confidence level:
p = .95;
%pre-allocate ? two observations = one degree of freedom needed to
%calculate confidence interval
%confidence_interval_upper = zeros((size(pred)-1));

for i = 2:size(pred)
nu = i-1 ; 
t = tinv(p,nu);
upper = pred(i) + t * (pred(i)/sqrt(i));  
confidence_interval_upper(i) = upper;
end



   

