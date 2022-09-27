
%%
clear 
clc
%% Slettes
lambda = randi(30,50,1);
date = (1:50)';

x = 10;
y = poisscdf(x, lambda);
prob = 1 - y;

figure
line(date,prob)
xlabel('Observation')
ylabel('Cumulative Probability')

%% 
% Set dates interval
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
T0 = 120;
T = 464;
date = Date(T0+1:T+1); 

%import lambda values
lambda = xlsread('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Prefered_PARX_model_forecast/Resultater_PARX/lambda_f_PARX_VaR.xlsx');    
%%
%Estimating probability of P(Y > x|lambda) (The probability that tomorrows
%realized value, Y, is greater than the value x given our estimated lambda.

%Setting the parameter limit of interest, x, (prob that Y is smaller than x)
x = 10;
P = poisscdf(x, lambda(T0+1:T+1));
%Convert to greater than x rather then below to estimate the predicted risk
%of getting a value higher than x. 
prob = 1 - P;

%Plotting our results. 
figure
line(date,prob, 'LineWidth', 0.8)  %Linewidth is adjusted to be thicker than standard  
xlabel('Date', 'FontSize', 14)
ylabel('Cumulative Probability', 'FontSize', 14)

ax = gca;
ax.YColor = 'k';
ax.FontSize = 14; 

grid minor

saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/VaR_10.png')

%% Same plot with quantiles

alpha = 0,025;
Q = poisscdf(alpha,lambda(T0+1:T+1))
