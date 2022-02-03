% PAR model estimation
%Version tilpasset efter Anders mail den 12/1.  
clc 
clear all 
close all

%%
%Loop to generate several simulations at once
for n = 1:1
% Definition of pre-determined parameters and coefficients
% Coefficients (The chosen parameters)
omega = 0.0645;
alpha = 0.415;
beta  = 0.577;
gamma = 0.000;


% Set output save location
filename = 'simu.xlsx';

% PAR(1,1)
%Benchmark model

% generated data
m = omega/(1-alpha-beta);
x = [];

% Preallocate size of data_y
data_y = zeros(199,1);

% Use first observation in data as 1st datapoint in simulation
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');
data_y = [hosp(1);data_y];

p = 1; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
data_y = data_y';
[r,T] = size(data_y);
g = p+q+r+1; % determine number of start values used

% Set initial values
theta0 = [0.5;0.3;0.1;0.26;0.64;0.05;0.71;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81];
theta0 = theta0(1:g); 

% Preallocate size of lambda
lambdas = zeros(200,1);

% To initialize the simulation, lambda0 and y0 is defined and lambdas is preallocated to 200 observations:
lambdas(1) = m;
data_y = data_y(1);
%%
for i = 1:199
    data   = poissrnd(m,1,1);     % Draw from conditional distribution (m will be updated iteratively)
    data_y = [data_y,data];       % Concatinate data with newly drawn datapoint
    %data_x = x(1:length(data_y)); % Size the exogenous data to match in size
    
    % Use parx_rob_simu function to get outputs
    [OUT, pearson_residuals, pred, probability, meany, rpit,confidence_interval] = output_parx_rob_simu(data_y,x,theta0,p,q,s,filename,omega, alpha, beta, gamma); 
    
    m = pred(end);     % Update lambda for new conditional poisson distibution 
    lambdas(i+1) = m;  % Save lambda in array
    
end 
%Output_PARX_simu applies the chosen parameters to data and so does not
%conduct maximum likelihood.
writematrix(lambdas, '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Outputs/PARX/lambdas_simulation.xlsx')

%% Plot predicted values vs. actual observations
% Set dates for plotting
t1 = 1;
t2 = 199;
Date = (t1:t2)';
data = data_y(2:end);

figure; 

bar(Date, data)
hold on
plot(Date,pred,'k', 'LineWidth',2)
ylabel('Count','FontSize',  14)
xlabel('Time','FontSize',  14)
grid minor
% hold off
ax = set(gca,'FontSize', 14);

legend('Simulated Data','Simulated PAR(1,1)', 'FontSize', 14);

saveas(gcf,sprintf('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/simu_PARX_2001%f.jpg',n))
% ax.FontSize = 14;hal


%% PLot of data distributions

% Number of bins in histogram
nbins=30;

figure;
    histogram(data_y,nbins)
    title('', 'FontSize', 14)
    xlabel('Data draw', 'FontSize', 14)
    ylabel('Counts', 'FontSize', 14)
    %title('Simulated data')
    
    ax = gca;
    ax.FontSize = 14; 
    grid minor
saveas(gcf,sprintf('/Users/krmmm/Documents/Dokumenter_Mac/simu_PARX_distribution%f.jpg',n))

figure;
    histogram(hosp,nbins)
    title('', 'FontSize', 14)
    xlabel('Data draw', 'FontSize', 14)
    ylabel('Counts', 'FontSize', 14)
    %title('Hosp dataset')
    
    ax = gca;
    ax.FontSize = 14; 
    grid minor
    
    
saveas(gcf,sprintf('/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/real_PARX_distribution%f.jpg',n))
end
    
%%
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