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
X=[]; % Defined as empty matrix to run a PAR frame work

T0 = 464;
Y = Y(:,1:T0); % y is shortend in order to allow us to run the model. 

% X(1,:)=max(0,diff(tested));

% To include exogenous variables into model X needs to be changed to
% X(:,1:t) in loop for estimations

% Determine length of data (number of days)
T = length(Y);

%% Initial values of theta used in numerical optimization computed

% Set Theta_init
theta_init = [0.2284;  0.1842; 0.1990; 0.1579];
%theta_init = [0.2284;  0.1842; 0.1990];

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
w  = 100; % rolling window size of 100 days
T0 = 100; % days prior to evaluation of first estimation (start model at day 101)

series = zeros(length(x)-w -T0,w); % Preallocation of series (100 days wide and rows until day 464 occurs => length of 264)
for i  = 1:length(x)-w -T0
    series(i,:) = x(i+w:i+w+T0-1);
end
% This series array holds the index of days occuring in each series
% e.g. first row contains day 101 until day 200
% e.g. second row contains day 102 until day 201 and so forth...
% window of series keeps rolling until index day 464 occurs

%%
%T0 = 120; % Set start observation for forecast

% Set options for the fminsearch function
tolerance = 1e-2;
options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun',tolerance,'TolX', tolerance, 'Display', 'iter','PlotFcns',@optimplotfval);

% Create empty matrices to fill data into
logL_tmp_T = zeros(1, length(x)-w -T0);
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