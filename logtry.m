%% PARX model estimation
clc 
clear all
close all

%% Import data
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');     
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_16.txt'); 
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_16.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_16.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');

%%
logcontact=reallog(contact);

%%
x=tested;
%x1=reallog(x);
x1=diff(x,2);
[h,pValue] = adftest(x1)
%%

temp=diff(temp);
contact=diff(contact);
tested=diff(tested);
vacc=diff(vacc);
varbrit=diff(varbrit);
vardelta=diff(vardelta);
%% Select data to use
y=hosp';
%X=[];
X(1,:)=temp;
X(2,:)=contact;
X(3,:)=tested;
X(4,:)=vacc;
X(5,:)=varbrit;
X(6,:)=vardelta;
%X(1,:)=strin;

%T0 = 14;
% 
%y = y(:,1:T0);
%X = X(:,1:T0);
 
%% PARX(1,1)
%Benchmark model

p = 1; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
[r,T] = size(X);
g = p+q+r+1; % determine number of start values used
y = y(:,1:T);

% Set initial values
theta0 = [0.5;0.3;0.1;0.3;0.4;0.1;0.8;0.5;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81];
theta0 = theta0(1:g);

% Set output save location
filename = 'Outputs/PARX_22/1308firstdifference.xlsx';

% Use parx_rob function to get outputs
[OUT, pearson_residuals, pred, zero_prob, meany, rpit,upper_pred] = output_parx_rob(y,X,theta0,p,q,s,filename); 

%Y(1:T0),X(:,1:T0)

%% Predicted defaults per month

% Set dates for plotting
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
k = p+1;
g = y(k:465);
Date = Date(k:end);

%%
%Creating the subsample length, removing observations at the start of the
%period, responding to the lag length
%k is lag length + 1, to start from the first observation with sufficient
%lags. g is the subsample length, the endpoint 465 is quite a "manual"
%soluition, and could benefit from being calculated by code determining the
%length of the full sample. 

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
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/predictions.jpg')

%% Predicted defaults per month

k = p+1;
g = y(k:465);

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

%% Residual Autocorrelation of pearson residuals

figure; 
    autocorr(pearson_residuals)
    %title('Residual Autocorrelation Function', 'FontSize', 20);
    legend('Predictions','True', 'FontSize', 14);
    xlabel('Lag', 'FontSize', 14);     
    ylabel('Sample Autocorrelation', 'FontSize', 14);
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
    %saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/ACF_residuals.jpg')
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
    plot(Date, zero_prob)
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
    grid on
    
    %set(fig1,'Position',[100 100 1000 500])
    
    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/zero_prob.jpg')

%% PIT

% Number of bins in histogram
nbins=25;

figure;
    histogram(rpit,nbins)
    title('Randomized Probability Integral Transform', 'FontSize', 14)
    %xlabel('Her', 'FontSize', 14)
    %ylabel('Her også', 'FontSize', 14)
    
    ax = gca;
    ax.FontSize = 14; 

% Kolmogorov-Sminov test
[h,p] = kstest(norminv(rpit))

saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/pit_histogram.jpg')
