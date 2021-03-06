%% PARX model estimation
clc 
clear all 
close all
%% generated data
m = 32.22; % mean
data = poissrnd(m,200,1);

% generated norm data
k = zeros(200,1);
k(1)= 0.5;
for i = 2:200
    k(i) = normrnd(k(i-1)*0.5,1);
end 

y = data';
x = k';
%x = [];

% model = 'PARX';
% %model = 'PARX'
% if model == 'PAR'
%     x = [];
% else
%     x = k';
% else

%% Coefficients (The chosen parameters)
% omega = 0.159;
% alpha = 0.435;
% beta  = 0.545;
% gamma = 0.0;

omega = 0.0645;
alpha = 0.415;
beta  = 0.577;
gamma = 0.0025;

t = 1:200; % time points

%% PARX(1,1)
%Benchmark model

p = 1; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
[r,T] = size(x);
g = p+q+r+1; % determine number of start values used

% Set initial values
%theta0 = [0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81];
theta0 = [0.5;0.3;0.1;0.26;0.64;0.05;0.71;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81];
theta0 = theta0(1:g);

% Set output save location
filename = 'simu.xlsx';

% Use parx_rob_simu function to get outputs
[OUT, pearson_residuals, pred, probability, meany, rpit,confidence_interval] = output_parx_rob_simu(y,x,theta0,p,q,s,filename,omega, alpha, beta, gamma); 

%Output_PARX_simu applies the chosen parameters to data and so does not
%conduct maximum likelihood.

%% Plot predicted values vs. actual observations
% Set dates for plotting
t1 = 1;
t2 = 199;
Date = (t1:t2)';
data = data(2:end);

figure; 

bar(Date, data)
hold on
plot(Date,pred,'k', 'LineWidth',1.5)
ylabel('Count','FontSize',  14)
xlabel('Time','FontSize',  14)
grid minor
% hold off

title('Specified PARX(1,1) vs Generated Data','FontSize', 14);
legend('Generated Data','Specified PARX(1,1)', 'FontSize', 14);

saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/simu_PARX_1711.jpg')
% ax.FontSize = 14;
%% PLot of data distributions

% Number of bins in histogram
nbins=30;

figure;
    histogram(data,nbins)
    title('', 'FontSize', 14)
    xlabel('Data draw', 'FontSize', 14)
    ylabel('Counts', 'FontSize', 14)
    
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
    
    %% PARX model estimation

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting parameters to data insted of forcing predeterment parametervalues

clc 
clear all 
close all
%% generated data
m = 32.22; % mean
data = poissrnd(m,200,1);

% generated norm data
k = zeros(200,1);
k(1)= 0.5;
for i = 2:200
    k(i) = normrnd(k(i-1)*0.5,1);
end 

y = data';
X = k';
%x = [];

% model = 'PARX';
% %model = 'PARX'
% if model == 'PAR'
%     x = [];
% else
%     x = k';
% else
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
filename = 'not_simu.xlsx';

% Use parx_rob_simu function to get outputs
[OUT, pearson_residuals, pred, probability, meany, rpit,confidence_interval, theta] = output_parx_rob(y,X,theta0,p,q,s,filename); 

%Output_PARX_simu applies the chosen parameters to data and so does not
%conduct maximum likelihood.
%% Plotting predictions

% Set dates for plotting
t1 = 1;
t2 = 199;
Date = (t1:t2)';
data = data(2:end);
% Plot predicted values vs. actual observations
figure; 

bar(Date, data)
hold on
plot(Date,pred,'k', 'LineWidth',1.5)
ylabel('Count','FontSize',  14)
xlabel('Time','FontSize',  14)
% hold off

title('Specified PAR(1,1) vs Generated Data','FontSize', 14);
legend('Generated Data','Specified PAR(1,1)', 'FontSize', 14);
% ax.FontSize = 14;
% saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/predictions.jpg')