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

%% ikketransformerettest
testednonnegative=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
%% Length of time series
T0 = 204;
%T0 = 236;
y=hosp';
y = y(:,62:T0); % y is shortend in order to allow us to run the model. 
%X = X(:,1:T0);
%%
dummy = zeros(T0,1);
for i = 1:T0
    if(i >= 236)
        dummy(i) = 1;
    end
end
%% Select data to use
% X=[];
% X(1,:)=temp;
% X(2,:)=contact;
% X(3,:)=tested;
% X(4,:)=vacc;
% X(5,:)=varbrit;
% X(6,:)=vardelta;
% 
%New transformed variables
% 
% X(1,:)=abs(min(0,diff(temp)));
% X(1,:)=max(0,diff(temp));
% 
% X(3,:)=abs(min(0,diff(contact)));
% X(4,:)=max(0,diff(contact));
% 
% X(5,:)=abs(min(0,diff(tested))); 
X(1,:)=max(0,diff(tested));%Bliver signifikant
X = X(:,62:T0);

% X(1,:)=abs(diff(vacc)); 

% X(1,:)=abs(min(0,diff(varbrit)));
%  X(2,:)=max(0,diff(varbrit));

% X(1,:)=abs(min(0,diff(vardelta)));
% X(2,:)=max(0,diff(vardelta));

% X(1,:)=dummy;
%% PARX(1,1)
%Benchmark model

p = 1; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
[r,T] = size(X);
g = p+q+r+1; % determine number of start values used

% Set initial values
%theta0 = [0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81];
theta0 = [0.5;0.3;0.1;0.26;0.64;0.05;0.71;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81];
theta0 = theta0(1:g);

% Set output save location
filename = 'Outputs/PARX_22/2709_lockdown_62_204_short.xlsx';

% Use parx_rob function to get outputs
[OUT, pearson_residuals, pred, probability, meany, rpit,confidence_interval] = output_parx_rob(y,X,theta0,p,q,s,filename); 

%Y(1:T0),X(:,1:T0)
%%
writematrix(pred, '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Outputs/PARX/pred.xlsx')
%writematrix(, '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Outputs/PARX/lambda_f_2para.xlsx')
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
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/TestedP/predictions.jpg')

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
%% Predicted defaults per month
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2-2);

figure; 
u = y(p+1:end)
plot(Date,u,Date,pred,Date, confidence_interval, 'LineWidth',1.5)
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

