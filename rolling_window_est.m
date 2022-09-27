%% Rolling estimation and forecasting of COVID-19 hospitalizations using PARX(1,1) with X = tested(+) and PAR(1,1)
clear all
close all
clc
%% Import data
hosp     = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');     
temp     = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_18.txt'); 
contact  = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested   = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
vacc     = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_18.txt');
strin    = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit  = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAR MODEL
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select data for model
Y=hosp'; % Transpose to fit into out
X=[]; % Defined as empty matrix to run a PAR frame work

T0 = 465; %Why not 465 - because after doing first diffence only 464 datapoints are available 
Y = Y(:,1:T0); % y is shortend in order to allow us to run the model. 

%X(1,:)=max(0,diff(tested));

% To include exogenous variables into model X needs to be changed to
% X(:,1:t) in loop for estimations

% Determine length of data (number of days)
T = length(Y);
%Denne er vel overflødig?
%% Initial values of theta used in numerical optimization computed

% Set Theta_init
%theta_init = [0.2284;  0.1842; 0.1990; 0.1579];
theta_init = [0.2284;  0.1842; 0.1990];

% Amount of initial parameters need to be adjusted to amount of parameters included

% Determine size of theta_init
N_init = length(theta_init);

% Set parameter length omega + alpha + beta + X
parameters = 3;

% Create matrices to put values into
theta_est = zeros(parameters,T);

%% Iterative estimation and forecasting
% We re-estimate the model adding one additional observation at a time;
% these are then used to compute one-step ahead forecast of lambda_t

% Get rolling time window of 100 days
x  = 1:T; % Length of dataset (464 data points)
w  = 120; % rolling window size of 30 days (month)
T0 = 0; % days prior to evaluation of first estimation (e.g. start model at day 101)

series = zeros(floor(((length(x)-w)/w)*w),w); % Preallocation of series (100 days wide and rows that are multiple of 30 days window => length of 420)
for i  = 1:length(series)
    series(i,:) = x(i:i+w-1);
end
% This series array holds the index of days occuring in each series
% e.g. first row contains day 101 until day 200
% e.g. second row contains day 102 until day 201 and so forth...
% window of series keeps rolling until index day 464 occurs

%% 
%T0 = 120; % Optional: Set start observation for forecast

% Set options for the fminsearch function
tolerance = 1e-8;
options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun',tolerance,'TolX', tolerance); %,'PlotFcns',@optimplotfval

% Create empty matrices to fill data into
logL_tmp_T  = zeros(1, length(x)-w -T0);
theta_est_T = zeros(parameters, length(x)-w -T0);
theta_sqr_est_T = zeros(parameters, length(x)-w -T0);

% Uncomment to run for only 10 days
% T=T0+10;

%Specify laglength
p=1;
%%
% Loop over model estimations day by day, from T0 to end
for t = 1:size(series,1)

    % Search for values of theta. Save value for each iteration
    % The index of x X(1:t) is inserted when exogenous variables are considered and removed when only a PAR model is forecasted
    [theta_sqr_est_T(:,t),logL_tmp_T(t)] = fminsearch('logL_PARX',sqrt(theta_init),options,Y(series(t,1):series(t,w)),X,p,1);

    % Imposing numerical restrictions
    theta_est_T(:,t) = theta_sqr_est_T(:,t).^2;
    % Set theta init to found theta values
    theta_init = theta_est_T(:,t);
 
    % Forecast lambda value for next day
    lambda_f(t+1) = PARX_forecast(theta_est_T(:,t),Y(series(t,1):series(t,w)),X,1,1);

end 

%% Determine forecast error KLIC and MSFE
ForecastError = Y(w+1:end)-lambda_f(1:end-1);  % starts compairing at day windowlength +1 in data Y with first lambda forecast. Ends at second to last lambda forecast as this is where data Y ends.

logf = Y(w+1:end).*log(lambda_f(1:end-1)) - lambda_f(1:end-1);
%logf = Y(1:t+1).*log(lambda_f(1:t+1)) - lambda_f(1:t+1); , den gamle 

KLIC = zeros(1,t);
%før t+1
MSFE = zeros(1,t);

% This loop computes the accumulated MSFE + KLIC for each rolling window 
for i = 1:t
    %i = 1:t+1
    MSFE(i) = sum(ForecastError(1:i).^2)/(i); % "2" added as 401 is 0 whereby the observation is invalid. 
    KLIC(i) = sum(logf(2:i))/(i); % "2" added instead of 1, as the first two observations of lambda_f return inf, whereby it ruins the observations
%   KLIC(t) = -sum(logf(T0+1:t))/(t-T0); %Furthermore - is removed because
%   it is not in the original equation , but everything becomes negative if
%   it is not there? maybe absolute value
end

% Here the specific measure of the model in question is saved to provide data for comparringson plot further down 
writematrix(MSFE', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/MSFE_PAR.xlsx')
writematrix(KLIC', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/KLIC_PAR.xlsx')

%% Plot stability of parameters and forecast error

% Set dates for plotting interval
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

% Plot forecast errors
figure;
    yyaxis left
    plot(Date(1:t),MSFE)
    %ylabel('Parameter value')
    %plot(Date(1:t+1),MSFE)

    yyaxis right
    plot(Date(1:t),KLIC)
    %ylabel('Parameter value')
    %plot(Date(1:t+1),KLIC)
    
    %title('Out-of-sample fit: MSE and Score evaluation', 'FontSize', 14)
    legend('MSFE PAR','KLIC PAR', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARX MODEL
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

% Import data
hosp     = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');     
temp     = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_18.txt'); 
contact  = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested   = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
vacc     = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_18.txt');
strin    = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit  = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');

%% Select data for model
Y = hosp'; % Transpose to fit into out
X = max(0,diff(tested))';  % positive change in tested procentage

T0 = 464;
Y  = Y(:,1:T0); % y is shortend in order to allow us to run the model. 

% Determine length of data (number of days)
T = length(Y);

% Initial values of theta used in numerical optimization computed

% Set Theta_init
theta_init = [0.0645;  0.415; 0.577;  0.0025];
%theta_init = [0.2284;  0.1842; 0.1990;  0.1579];
%theta_init = [0.2284;  0.1842; 0.1990];
% Amount of initial parameters need to be adjusted to amount of parameters included

% Determine size of theta_init
N_init = length(theta_init);

% Set parameter length omega + alpha + beta + X
parameters = 4;

% Create matrices to put values into
theta_est = zeros(parameters,T);

%% Iterative estimation and forecasting
% We re-estimate the model adding one additional observation at a time;
% these are then used to compute one-step ahead forecast of lambda_t

% Get rolling time window of 100 days
x  = 1:T; % Length of dataset (464 data points)
w  = 30; % rolling window size of 30 days (month)
T0 = 0; % days prior to evaluation of first estimation (e.g. start model at day 101)

series = zeros(floor(((length(x)-w)/w)*w),w); % Preallocation of series (100 days wide and rows that are multiple of 30 days window => length of 420)
for i  = 1:length(series)
    series(i,:) = x(i:i+w-1);
end

% This series array holds the index of days occuring in each series
% e.g. first row contains day 101 until day 200
% e.g. second row contains day 102 until day 201 and so forth...
% window of series keeps rolling until index day 464 occurs

%%
%T0 = 120; % Optional: Set start observation for forecast

% Set options for the fminsearch function
tolerance = 1e-8;
options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun',tolerance,'TolX', tolerance); %,'PlotFcns',@optimplotfval

% Create empty matrices to fill data into
logL_tmp_T  = zeros(1, length(x)-w -T0);
theta_est_T = zeros(parameters, length(x)-w -T0);
theta_sqr_est_T = zeros(parameters, length(x)-w -T0);

% Uncomment to run for only 10 days
% T=T0+10;

%Specify laglength
p=1;

lambda_f_PARX=[]
% Loop over model estimations day by day, from T0 to end
for t = 1:size(series,1)

    % Search for values of theta. Save value for each iteration
    % The index of x X(1:t) is inserted when exogenous variables are considered and removed when only a PAR model is forecasted
    [theta_sqr_est_T(:,t),logL_tmp_T(t)] = fminsearch('logL_PARX',sqrt(theta_init),options,Y(series(t,1):series(t,w)),X(series(t,1):series(t,w)),p,1);

    % Imposing numerical restrictions
    theta_est_T(:,t) = theta_sqr_est_T(:,t).^2;

    % Set theta init to found theta values
    theta_init = theta_est_T(:,t);
 
    % Forecast lambda value for next day
    lambda_f_PARX(t+1) = PARX_forecast(theta_est_T(:,t),Y(series(t,1):series(t,w)),X(series(t,1):series(t,w)),1,1);
    
end 

%% Determine forecast error KLIC and MSFE
ForecastError_PARX = Y(w+1:end)-lambda_f_PARX(1:end-1);
logf_PARX= Y(w+1:end).*log(lambda_f_PARX(1:end-1)) - lambda_f_PARX(1:end-1);
%ForecastError_PARX = Y(w+1:length(lambda_f_PARX)+w)-lambda_f_PARX;  % starts compairing at day 31 in data Y with first lambda forecast

writematrix(ForecastError_PARX', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Outputs/PARX/forecasterror2.xlsx')

KLIC_PARX = zeros(1,t);
MSFE_PARX = zeros(1,t);

% This loop computes the accumulated MSFE + KLIC for each rolling window 
for i = 1:t
    MSFE_PARX(i) = sum(ForecastError_PARX(1:i).^2)/(i); % "2" added as 401 is 0 whereby the observation is invalid. 
    KLIC_PARX(i) = sum(logf_PARX(2:i))/(i); % "2" added instead of 1, as the first two observations of lambda_f return inf, whereby it ruins the observations
%   KLIC(t) = -sum(logf(T0+1:t))/(t-T0); %Furthermore - is removed because
%   it is not in the original equation
end

% Here the specific measure of the model in question is saved to provide data for comparringson plot further down 
writematrix(MSFE_PARX', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/MSFE_PARX.xlsx')
writematrix(KLIC_PARX', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/KLIC_PARX.xlsx')

% Plot stability of parameters and forecast error

%% Figure that plots KLIC and MSFE for a single model 

% Set dates for plotting interval
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

% Plot forecast errors
figure;
    yyaxis left
    %plot(Date(1:t+1),MSFE); hold on
    plot(Date(1:t),MSFE_PARX); hold off
    %ylabel('Parameter value')

    yyaxis right
    %plot(Date(1:t+1),KLIC); hold on
    plot(Date(1:t),KLIC_PARX); hold off
    %ylabel('Parameter value')
    
    %title('Out-of-sample fit: MSE and Score evaluation', 'FontSize', 14)
    legend('MSFE PARX','KLIC PARX', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
    %% Figure that compares KLIC and MSFE for different models

    %MSFE
MSFE_PAR = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/MSFE_PAR.xlsx'); 
MSFE_PARX = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/MSFE_PARX.xlsx');

T0 = w+1; % The point of first forecast
T = length(X); % The length of the included exogeneus variable is used as this series is shorter than the series of the dependent variable due to FD

% Set dates for plotting interval
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

figure; 
    plot(Date(T0:T),MSFE_PAR(1:end-1),Date(T0:T),MSFE_PARX(1:end))
%    ylabel('Parameter value')
%   title('Comparing MSFE PARX and MSFE PAR', 'FontSize', 14)
    legend('MSFE PAR','MSFE PARX', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/Comparing_MSFE.jpg')
%%
 % KLIC
KLIC_PAR = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/KLIC_PAR.xlsx'); 
KLIC_PARX = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/KLIC_PARX.xlsx');
%%
figure; 
    plot(Date(T0:T),KLIC_PAR(1:end-1),Date(T0:T),KLIC_PARX(1:end))
%    ylabel('Parameter value')
    
%   title('Comparing KLIC PARX and KLIC PAR', 'FontSize', 14)
    legend('KLIC PAR','KLIC PARX', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/Comparing_KLIC.jpg')

%% Figure that illustrates best model performance
scores = zeros(length(MSFE_PARX),2);

for i =1:length(MSFE_PARX)
%     if par_est1(i,1) > alpha  % Model 1
%         par_est1(i,2) = 0;
%     end
%         
%     if par_est2(i,1) > alpha % Model 2
%         par_est2(i,2) = 0;
%     end
    
% What model is best    
    if MSFE_PAR(i) < MSFE_PARX(i)
        scores(i,1) = 1;
    else
        scores(i,2) = 2;  
    end
end

days = linspace(1,length(MSFE_PARX),length(MSFE_PARX));

figure;
plot(days,scores(:,1),'linestyle','none','marker','+','color','b','LineWidth',0.7); hold on
plot(days,scores(:,2),'linestyle','none','marker','+','color','r','LineWidth',0.7); hold off
ylim([0.5,2.5]);
xlim([0,days(end)])
xlabel('Days')
yticks([1 2])
yticklabels({'PAR','PARX'})
ax = gca;
ax.FontSize = 14; 
grid minor
% yticks([0 50 100])
% yticklabels({'y = 0','y = 50','y = 100'})
%     plot(Date(z),hospp(z),'linestyle','none','marker','*','LineWidth',0.7)
%     ylim([3,10])

%% FEJL
figure;

plot(days,scores(:,2),'linestyle','none','marker','+','color','b','LineWidth',0.7); hold on
plot(days,scores(:,1),'linestyle','none','marker','+','color','r','LineWidth',0.7); hold off

ylim([0.5,2.5]);
xlim([0,days(end)])
xlabel('Days')
yticks([1 2])
yticklabels({'PAR','PARX'})
ax = gca;
ax.FontSize = 14; 
grid minor
% yticks([0 50 100])
% yticklabels({'y = 0','y = 50','y = 100'})
%     plot(Date(z),hospp(z),'linestyle','none','marker','*','LineWidth',0.7)
%     ylim([3,10])


%%
 
 % Plot parameter stability 
figure; 
    plot(Date(T0+1:T),theta_est_T(1,1:T-T0))
    %title('Parameter stability')
    legend('\omega', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

 saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/Parmeters_stability_omega_PARX_Rolling.jpg')
figure; 
    plot(Date(T0+1:T),theta_est_T(2,1:T-T0))
    %title('Parameter stability')
    legend('\alpha_1', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
 saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/Parmeters_stability_alpha_PARX_Rolling.jpg')
figure; 
    plot(Date(T0+1:T),theta_est_T(3,1:T-T0))
    %title('Parameter stability')
    legend('\beta', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

     saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/Parmeters_stability_beta_PARX_Rolling.jpg')

figure; 
    plot(Date(T0+1:T),theta_est_T(4,1:T-T0))
    %title('Parameter stability')
    legend('Tested(+)', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

     saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/Parmeters_stability_Tested+_PARX_Rolling.jpg')
     
%%

figure; 
     [AX,H1,H2] = plotyy(Date(T0:T),theta_est_T(4,1:end),Date(T0:T),MSFE_PARX(1:end))
     %title('Parameter stability')
     legend('Tested(+)','MSFE PARX')
     xlabel('Date', 'FontSize', 14)
     ylabel(AX(1),'Parameter value','FontSize', 14)
     ylabel(AX(2),'MSFE','Fontsize',14)
     set(AX,'Fontsize',14)
     %ax = gca;
     %ax.FontSize = 14;
     grid minor
     saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/Parmeters_stability_Tested+MSFE_PARX_Rolling.jpg')

     
%% Figure to illustrate best model performance MSFE

%Creating Date for our observation period. 
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

% Checking if exogenous varaible is significant (aplha = 0.5)
% alpha = 1.64;
alpha = 0;

MSFE_PAR = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/MSFE_PAR.xlsx'); 
MSFE_PARX = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/MSFE_PARX.xlsx');

%Cutting the last observation from our PAR Model due to reasons. 
MSFE_PAR = MSFE_PAR(1:end-1,:)
T0 = 100
T = 464

scores = zeros(length(MSFE_PARX),3);
%scores = zeros(length(MSFE_PARX),2);

for i =1:length(MSFE_PARX)
%     if par_est1(i,1) > alpha  % Model 1
%         par_est1(i,2) = 0;
%     end
%         
%     if par_est2(i,1) > alpha % Model 2
%         par_est2(i,2) = 0;
%     end
    
% What model is best    
    if MSFE_PAR(i) < MSFE_PARX(i)
        scores(i,1) = 1;
    elseif MSFE_PAR(i) == MSFE_PARX(i)
        scores(i,3) = 3;
    else
        scores(i,2) = 2;  
    end
end

figure;
plot(Date(T0+1:T),scores(:,1),'linestyle','none','marker','+','color','#4DBEEE','LineWidth',0.7); hold on
plot(Date(T0+1:T),scores(:,2),'linestyle','none','marker','+','color','#EDB120','LineWidth',0.7); hold on
plot(Date(T0+1:T),scores(:,3),'linestyle','none','marker','+','color','#77AC30','LineWidth',0.7); hold off


ylim([0.5,3.5]);
%xlim([0,days(end)])
xlabel('Date')
yticks([1 2 3])
yticklabels({'PAR','PARX','PAR=PARX'})
ax = gca;
ax.FontSize = 14; 
grid off
% yticks([0 50 100])
% yticklabels({'y = 0','y = 50','y = 100'})
%     plot(Date(z),hospp(z),'linestyle','none','marker','*','LineWidth',0.7)
%     ylim([3,10])

%
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/MSFE_Comparison.png')
%% Figure to illustrate best model performance KLIC!!!

%Creating Date for our observation period. 
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

% Checking if exogenous varaible is significant (aplha = 0.5)
% alpha = 1.64;
alpha = 0;

KLIC_PAR = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/KLIC_PAR.xlsx'); 
KLIC_PARX = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/KLIC_PARX.xlsx');

%Cutting the last observation from our PAR Model due to reasons. 
KLIC_PAR = KLIC_PAR(1:end-1,:)
T0 = 100
T = 464

scores = zeros(length(KLIC_PARX),3);
%scores = zeros(length(MSFE_PARX),2);

for i =1:length(KLIC_PARX)
%     if par_est1(i,1) > alpha  % Model 1
%         par_est1(i,2) = 0;
%     end
%         
%     if par_est2(i,1) > alpha % Model 2
%         par_est2(i,2) = 0;
%     end
    
% What model is best    
    if KLIC_PAR(i) > KLIC_PARX(i)
        scores(i,1) = 1;
    elseif KLIC_PAR(i) == KLIC_PARX(i)
        scores(i,3) = 3;
    else
        scores(i,2) = 2;  
    end
end

days = linspace(1,length(KLIC_PAR),length(KLIC_PAR));

figure;
plot(Date(T0+1:T),scores(:,1),'linestyle','none','marker','+','color','#4DBEEE','LineWidth',0.7); hold on
plot(Date(T0+1:T),scores(:,2),'linestyle','none','marker','+','color','#EDB120','LineWidth',0.7); hold on
plot(Date(T0+1:T),scores(:,3),'linestyle','none','marker','+','color','#77AC30','LineWidth',0.7); hold off


ylim([0.5,3.5]);
%xlim([0,days(end)])
xlabel('Date')
yticks([1 2 3])
yticklabels({'PAR','PARX','PAR=PARX'})
ax = gca;
ax.FontSize = 14; 
grid off
% yticks([0 50 100])
% yticklabels({'y = 0','y = 50','y = 100'})
%     plot(Date(z),hospp(z),'linestyle','none','marker','*','LineWidth',0.7)
%     ylim([3,10])
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/RollingWindow/KLIC_Comparison.png')