%% Rolling estimation and forecasting of COVID-19 hospitalizations using PARX(1,1) with X = tested(+) and PAR(1,1)
clear all
close all
clc
%% Import data
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');     
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_18.txt'); 
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_18.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');

%% Select data for model
Y=hosp'; % Transpose to fit into out
%X=[]; % Defined as empty matrix to run a PAR frame work

T0 = 464;
Y = Y(:,1:T0); % y is shortend in order to allow us to run the model. 

X(1,:)=max(0,diff(tested));

% To include exogenous variables into model X needs to be changed to
%X(:,1:t) in loop for estimations

% Determine length of data (number of days)
T = length(Y);

%% Initial values of theta used in numerical optimization computed

% Set Theta_init
theta_init = [0.02284;  0.01842; 0.01990; 0.01579];
%theta_init = [0.2284;  0.1842; 0.1990; 0.1579];
%theta_init = [0.2284;  0.1842; 0.1990];

% Determine size of theta_init
N_init = length(theta_init);

% Set parameter length omega + alpha + beta + X
parameters = 4;

% Create matrices to put values into
theta_est = zeros(parameters,T);

%% Initialization of algorithm - OBS
%TOLERANCE ER MEGET LAV FOR AT TESTE AT KODEN K�RER. TolFun samt TolX skal
%s�ttes op inden koden k�res 'rigtigt'. Default v�rdien er 1e-4. 
tolerance = 1e-8;
options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun',tolerance,'TolX', tolerance); %'Display', 'iter','PlotFcns',@optimplotfval
T0 = 365; % Set number of days to include

% Search for minimum values of theta 
[theta_sqr_est, logL_tmp] = fminsearch('logL_PARX', sqrt(theta_init), options, Y(1:T0),X(1:T0),1,1);

% Imposing numerical restrictions
theta_est = theta_sqr_est.^2;
% Set theta init to found theta values
theta_init = theta_est;

%Specify lag-length
p=1

% Forecast value T+1 using estimated thate values
lambda_f(T0+1) = PARX_forecast(theta_est,Y(1:T0),X(1:T0),p,1); 
 
%% Iterative estimation and forecasting
% We re-estimate the model adding one additional observation at a time;
% these are then used to compute one-step ahead forecast of lambda_t

T0 = 120; % Set start observation for forecast

% Set options for the fminsearch function
tolerance = 1e-8;
options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun',tolerance,'TolX', tolerance); % 'Display', 'iter','PlotFcns',@optimplotfval
% Create empty matrices to fill data into
logL_tmp_T = zeros(1, T-T0);
theta_est_T = zeros(parameters, T-T0);
theta_sqr_est_T = zeros(parameters, T-T0);

% Uncomment to run for only 10 days
% T=T0+10;

%Specify laglength
p=1;

% Loop over model estimations day by day, from T0 to end
for t = T0+1:T

    % Search for values of theta. Save value for each iteration
    % The index of x X(1:t) is inserted when exogenous variables are considered and removed when only a PAR model is forecasted
    [theta_sqr_est_T(:,t-T0),logL_tmp_T(t-T0)] = fminsearch('logL_PARX',sqrt(theta_init),options,Y(1:t),X(:,1:t),p,1);

    % Imposing numerical restrictions
    theta_est_T(:,t-T0) = theta_sqr_est_T(:,t-T0).^2;
    % Set theta init to found theta values
    theta_init = theta_est_T(:,t-T0);

    % Forecast lambda value for next day
    lambda_f(t+1) = PARX_forecast(theta_est_T(:,t-T0),Y(1:t),X(:,1:t),1,1);

end 

%% Save Matrices

writematrix(lambda_f', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/lambda_f_2para_PARX_2.xlsx')
writematrix(logL_tmp_T, '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/logL_tmp_T_PAR.xlsx')
writematrix(Y', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/Y_PAR.xlsx')
writematrix(theta_est_T', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/theta_est_T_PAR.xlsx')
writematrix(theta_sqr_est_T', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/theta_sqr_est_T_2para_PAR.xlsx')

%% Determine forecast error KLIC and MSFE
ForecastError(T0+1:T) = Y(T0+1:T)-lambda_f(T0+1:T);
logf(T0+1:T) = Y(T0+1:T).*log(lambda_f(T0+1:T)) - lambda_f(T0+1:T);
 
KLIC = zeros(1,T);
MSFE = zeros(1,T);
 
for t = T0+1:T
    MSFE(t) = sum(ForecastError(T0+1:t).^2)/(t-T0); % "2" added as 401 is 0 whereby the observation is invalid. 
    KLIC(t) = sum(logf(T0+2:t))/(t-T0); % "2" added instead of 1, as the first two observations of lambda_f return inf, whereby it ruins the observations
%   KLIC(t) = -sum(logf(T0+1:t))/(t-T0); %Furthermore - is removed because
%   it is not in the original equation
end

writematrix(MSFE', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/MSFE_PAR.xlsx')
writematrix(KLIC', '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/KLIC_PAR.xlsx')
%% Plot stability of parameters and forecast error
close all;

% Set dates for plotting interval
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

% Plot forecast errors
figure;
    yyaxis left
    plot(Date(T0+1:T),MSFE(T0+1:T))
    ylabel('Parameter value')

    yyaxis right
    plot(Date(T0+1:T),KLIC(T0+1:T))
    ylabel('Parameter value')
    
    title('Out-of-sample fit: MSE and Score evaluation', 'FontSize', 14)
    legend('MSFE','KLIC', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
 saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/Forecast_errors_PAR.jpg')
% Plot parameter stability 
figure; 
    plot(Date(T0+1:T),theta_est_T(1,1:T-T0))
    title('Parameter stability')
    legend('\omega', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

 saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/Parmeters_stability_omega_PAR.jpg')
figure; 
    plot(Date(T0+1:T),theta_est_T(2,1:T-T0))
    title('Parameter stability')
    legend('\alpha_1', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
 saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/Parmeters_stability_alpha_PAR.jpg')
figure; 
    plot(Date(T0+1:T),theta_est_T(3,1:T-T0))
    title('Parameter stability')
    legend('\beta', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

     saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/Parmeters_stability_beta_PAR.jpg')
% figure; plot(Date(T0+1:T),theta_est_T(4,1:T-T0))
% %figure; plotyy(Date(T0+1:T),theta_est_T(4,1:T-T0),Date(T0+1:T),theta_est_T(5,1:T-T0))
%     title('Parameter stability')
%     legend('Tested (+)','Vaccines')
%     xlabel('Date', 'FontSize', 14)
%     ylabel('Parameter value', 'FontSize', 14)
%     ax = gca;
%     ax.FontSize = 14; 
%     grid minor
% 
% figure; plotyy(Date(T0+1:T),theta_est(7,T0+1:T),Date(T0+1:T),theta_est(8,T0+1:T))
% title('Parameter stability')
% legend('Test','VAX')
%% Comparing MSFE model

T0 = 120;
T = length(MSFE_PAR);

% Set dates for plotting interval
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

MSFE_PAR = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/MSFE_PAR.xlsx'); 
MSFE_PARX = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PARX/MSFE_PARX.xlsx');

figure; 
    plot(Date(T0+1:T),MSFE_PAR(T0+1:T),Date(T0+1:T),KLIC_PARX(T0+1:T))
    ylabel('Parameter value')
    
%   title('Comparing MSFE PARX and MSFE PAR', 'FontSize', 14)
    legend('MSFE PAR','MSFE PARX', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Comparing_MSFE.jpg')
%% Comparing KLIC model
KLIC_PAR = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/KLIC_PAR.xlsx'); 
KLIC_PARX = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PARX/KLIC_PARX.xlsx');

figure; 
    plot(Date(T0+1:T),KLIC_PAR(T0+1:T),Date(T0+1:T),KLIC_PARX(T0+1:T))
    ylabel('Parameter value')
    
%   title('Comparing KLIC PARX and KLIC PAR', 'FontSize', 14)
    legend('KLIC PAR','KLIC PARX', 'FontSize', 14)
    xlabel('Date', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

    saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Comparing_KLIC.jpg')
%% OBS NY TEST

for i=1:length(lambda_f) %her har du tilf�jet -1 fordi?
    c(i)= 1/poisscdf(10,lambda_f(i)); %change for LESS than 10 hospitilazation. For MORE than 1 - xxx
end

%% Rolling forecast experiment
% We re-estimate the model adding one additional observation at a time;
% these are then used to compute one-step ahead forecast of lambda_t

T0 = 120; % Set start observation for forecast

% Set options for the fminsearch function
tolerance = 1e-8;
options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun',tolerance,'TolX', tolerance, 'Display', 'iter','PlotFcns',@optimplotfval);

% Create empty matrices to fill data into
logL_tmp_T = zeros(1, T-T0);
theta_est_T = zeros(parameters, T-T0);
theta_sqr_est_T = zeros(parameters, T-T0);

% Uncomment to run for only 10 days
T=T0+10;

%Specify laglength
p=1;


% loops through dataset index 1 to length of T
% the loop removes the first index dynamically
for t = T0+1:T
    
    %takes the 50th first elements from dataset
    for h = T0:50

        % Search for values of theta. Save value for each iteration
        % The index of x X(1:t) is inserted when exogenous variables are considered and removed when only a PAR model is forecasted
        [theta_sqr_est_T(:,t-T0),logL_tmp_T(t-T0)] = fminsearch('logL_PARX',sqrt(theta_init),options,Y(1:t),X,p,1);
        
        % Imposing numerical restrictions
        theta_est_T(:,t-T0) = theta_sqr_est_T(:,t-T0).^2;
    
        % Set theta init to found theta values
        theta_init = theta_est_T(:,t-T0);

        % Forecast lambda value for next day
        lambda_f(t+1) = PARX_forecast(theta_est_T(:,t-T0),Y(1:t),X,1,1);
    end
    %remover first element
    t(2:end)

end 

%%
%createRollingWindow

%% Matlab version

% Simulate time points
x = linspace(1,463,463);

% Prealloacte
series   = zeros(363,100);
par_est1 = zeros(363,2);  %loglike, exo
par_est2 = zeros(363,2);


%%
%har lige �ndret til 1 og 99 for at pr�ve
for i = 1:length(x)-99
    series(i,:) = x(i:i+99);
end
    
%model simulation
par_est1(:,:) = randn(363,2);
par_est2(:,:) = randn(363,2);

%%
% Checking if exogenous varaible is significant (aplha = 0.5)
%alpha = 1.64;
alpha = 0;

MSFE_PARX = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PARX/MSFE_PARX.xlsx');
MSFE_PAR = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PAR/MSFE_PAR.xlsx');
%%
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

days = linspace(1,length(MSFE_PAR),length(MSFE_PAR));
%%
figure;
plot(days,scores(:,1),'linestyle','none','marker','+','color','r','LineWidth',0.7); hold on
plot(days,scores(:,2),'linestyle','none','marker','+','color','b','LineWidth',0.7); hold off
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
