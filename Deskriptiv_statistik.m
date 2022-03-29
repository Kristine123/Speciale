clc
clear
%%
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');

%%
%temp=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/Temperatur_366.txt');
%X=temp;  
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_18.txt'); 
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_18.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');

%%

%Deskriptiv statetstik 

%% Temperature
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
x = Date;
y = hosp;
figure;
yyaxis left
plot(x,y)
xlim(datetime([2020 2021],[03 08],1))

ylabel('Daily Hospitalizations', 'FontSize', 14)
z = temp;
yyaxis right
plot(x,z)

xlabel('Date', 'FontSize', 14)
ylabel('Temperature °C', 'FontSize', 14)
ax = gca;
ax.FontSize = 14; 
grid minor

%title('Hospitalizations & Temperature')
legend('Hospitalizations','Temperature (rhs)', 'FontSize',14)
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/Temperature.jpg')
%% Temperature
% temp1=abs(min(0,diff(temp)));
% temp2=max(0,diff(temp));
% 
% t1 = datetime(2020,03,26,8,0,0);
% t2 = datetime(2021,07,3,8,0,0);
% Date = (t1:t2-1);
% figure;
% yyaxis left
% plot(Date,temp1)
% xlim(datetime([2020 2021],[03 08],1))
% 
% ylabel('Daily Hospitalizations', 'FontSize', 14)
% yyaxis right
% plot(Date,temp2)
% 
% xlabel('Date', 'FontSize', 14)
% ylabel('Temperature °C', 'FontSize', 14)
% ax = gca;
% ax.FontSize = 14; 
% grid minor
% 
% %title('Hospitalizations & Temperature')
% legend('Hospitalizations','Temperature (rhs)', 'FontSize',14)
% saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/Temperature.jpg')
%% Restrictions 
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
figure;
yyaxis left
plot(Date,hosp)
xlim(datetime([2020 2021],[03 08],1))

ylabel('Daily Hospitalizations', 'FontSize', 14)
yyaxis right
plot(Date,strin)

xlabel('Date', 'FontSize', 14)
ylabel('Stringency Index', 'FontSize', 14)
ax = gca;
ax.FontSize = 14; 
grid minor

%title('Hospitalizations & Temperature')
legend('Hospitalizations','Restrictions (rhs)', 'FontSize',14)
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/Restrictions.jpg')
%% Contact Number 
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
figure;
yyaxis left
plot(Date,hosp)
xlim(datetime([2020 2021],[03 08],1))

ylabel('Daily Hospitalizations', 'FontSize', 14)
yyaxis right
plot(Date,contact)

xlabel('Date', 'FontSize', 14)
ylabel('Contact Number', 'FontSize', 14)
ax = gca;
ax.FontSize = 14; 
grid minor

%title('Hospitalizations & Temperature')
legend('Hospitalizations','Contact Number (rhs)', 'FontSize',14)
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/Contact.jpg')
%% Tested
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
figure;
yyaxis left
plot(Date,hosp)
xlim(datetime([2020 2021],[03 08],1))

ylabel('Daily Hospitalizations', 'FontSize', 14)
yyaxis right
plot(Date,tested)

xlabel('Date', 'FontSize', 14)
ylabel('Percent Tested', 'FontSize', 14)
ax = gca;
ax.FontSize = 14; 
grid minor

%title('Hospitalizations & Temperature')
legend('Hospitalizations','Tested (rhs)', 'FontSize',14)
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/Tested.jpg')
%% Vaccins
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
figure;
yyaxis left
plot(Date,hosp)
xlim(datetime([2020 2021],[03 08],1))

ylabel('Daily Hospitalizations', 'FontSize', 14)
yyaxis right
plot(Date,vacc)

xlabel('Date', 'FontSize', 14)
ylabel('Percent Vaccinated', 'FontSize', 14)
ax = gca;
ax.FontSize = 14; 
grid minor

%title('Hospitalizations & Temperature')
legend('Hospitalizations','Vaccinated (rhs)', 'FontSize',14)
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/Vaccinated.jpg')

%% Variants
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
figure;
yyaxis left
plot(Date,hosp)
xlim(datetime([2020 2021],[03 08],1))

ylabel('Daily Hospitalizations', 'FontSize', 14)
yyaxis right

ax = gca;
ax.YColor = 'k';
ax.FontSize = 14; 
plot(Date ,varbrit ,Date ,vardelta,'g-','LineWidth',1)

xlabel('Date', 'FontSize', 14)
ylabel('Percent of Positive Cases', 'FontSize', 14)

grid minor

%title('Hospitalizations & Temperature')
legend('Hospitalizations','British - B.1.1.7 (rhs)','Delta - B.1.617.2 (rhs)','Location','northwest','FontSize',14)
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/Variations.jpg')
%% For the descriptive statistics ACF

figure;

autocorr(hosp, 'NumLags', 60)

%% For the descriptive Statistics PACF

figure; parcorr(hosp,'NumLags',60)

%% For the descriptive Statistics , plot of Hospitalizations
figure
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
x = Date;
y = hosp;
plot(x,y)
xlim(datetime([2020 2021],[03 08],1));
title('Hospitalizations');

%% Hospitalizations
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
figure;
plot(Date,hosp)
xlim(datetime([2020 2021],[03 08],1))
ylabel('Daily Hospitalizations', 'FontSize', 14)

ax = gca;
ax.YColor = 'k';
ax.FontSize = 14; 

xlabel('Date', 'FontSize', 14)

grid minor

%title('Hospitalizations & Temperature')
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer_2408/Hospitalizations.jpg')

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
 
% Plot parameter stability 
figure; 
    plot(Date(T0+1:T),theta_est_T(1,1:T-T0))
    title('Parameter stability')
    legend('\omega', 'FontSize', 20)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

figure; 
    plot(Date(T0+1:T),theta_est_T(2,1:T-T0))
    title('Parameter stability')
    legend('\alpha_1', 'FontSize', 20)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

figure; 
    plot(Date(T0+1:T),theta_est_T(3,1:T-T0))
    title('Parameter stability')
    legend('\beta', 'FontSize', 20)
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

figure; plot(Date(T0+1:T),theta_est_T(4,1:T-T0))
%figure; plotyy(Date(T0+1:T),theta_est_T(4,1:T-T0),Date(T0+1:T),theta_est_T(5,1:T-T0))
    title('Parameter stability')
    legend('Tested (-)','Vaccines')
    xlabel('Date', 'FontSize', 14)
    ylabel('Parameter value', 'FontSize', 14)
    ax = gca;
    ax.FontSize = 14; 
    grid minor

% 
% figure; plotyy(Date(T0+1:T),theta_est(7,T0+1:T),Date(T0+1:T),theta_est(8,T0+1:T))
% title('Parameter stability')
% legend('Test','VAX')