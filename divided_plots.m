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
%% Create data and 2-by-1 tiled chart layout
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

x1 = Date(1:63);
x2 = Date(64:278);
x3 = Date(279:391);
x4 = Date(392:end);

y1 = hosp(1:63);
y2 = hosp(64:278);
y3 = hosp(279:391);
y4 = hosp(392:end);
tiledlayout(1,4)

% First plot
ax1 = nexttile;
plot(ax1,x1,y1)
title(ax1,'26th mar to 28th may 2020', 'FontSize', 14)
ylabel(ax1,'New Hospitalization')
xlim(datetime([2020 2020],[03 05],[26 28]))
ylim([0 250]);

ax = gca;
ax.FontSize = 14; 
grid minor
xtickangle(45)

% Second plot
ax2 = nexttile;
plot(ax2,x2,y2)
title(ax2,'29th may to 29th dec 2020', 'FontSize', 14)
xlim(datetime([2020 2020],[05 12],[29 29]))
ylim([0 250]);

ax = gca;
ax.FontSize = 14; 
grid minor
xtickangle(45)

% Third plot
ax3 = nexttile;
plot(ax3,x3,y3)
title(ax3,'30th dec 2020 to 21th apr 2021', 'FontSize', 14)
xlim(datetime([2020 2021],[12 04],[30 21]))
ylim([0 250]);
xtickangle(45)

ax = gca;
ax.FontSize = 14; 
grid minor

% Last plot
ax4 = nexttile;
plot(ax4,x4,y4)
title(ax4,'22th apr to 13th july 2021', 'FontSize', 14)
xlim(datetime([2021 2021],[04 07],[22 13]))
ylim([0 250]);

ax = gca;
ax.FontSize = 14; 
grid minor
xtickangle(45)