%% Creation of stationary exogeneus variables
%%
clc 
clear all
close all
%% Load data
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');     
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_16.txt'); 
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_16.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_16.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');
%%
% X=[];
X(1,:)=abs(min(0,diff(temp)));
X(2,:)=max(0,diff(temp));

X(3,:)=abs(min(0,diff(contact)));
X(4,:)=max(0,diff(contact));

X(5,:)=abs(min(0,diff(tested)));
X(6,:)=max(0,diff(tested));

X(7,:)=abs(diff(vacc));

X(8,:)=abs(min(0,diff(varbrit)));
X(9,:)=max(0,diff(varbrit));

X(10,:)=abs(min(0,diff(vardelta)));
X(11,:)=max(0,diff(vardelta));

%% Initial test
%ADF test
kk = abs(diff(vacc));
% k = log(max(kk,1e-99));
%x=diff(x);
[h,pValue] = adftest(kk)

%% Initial plot
figure
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
%x = Date;
y = kk;
plot(y)
%xlim(datetime([2020 2021],[03 08],1));
title('Initial plot');

%% Transform
%First difference

y= diff(y);

%Could also try log!
%Or 2nd diff

% %Decomposing into positive and negative values:
% %X = y;
% idx = X<0; % create logical index
% neg = X(idx);
% pos = X(~idx);
% %Takeing the absolute value of the negative series
% abs_neg = abs(neg);

%Du mangler at fylde nuller ind på de manglende pladser...
%ELLER - sørge for at dimensionerne i X beholdes pos = max(0, X);
%% 2nd try
pos = max(0,y);
neg = min(0,y);
% Take the absolute value of negative series:
abs_neg = abs(neg);

%% Plot to see effect

figure
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1+1:t2);
plot(Date,pos)
xlim(datetime([2020 2021],[03 08],1));
title('after transformation positive plot');


figure
plot(Date,abs_neg)
xlim(datetime([2020 2021],[03 08],1));
title('after transformation negative plot');

%% Test
%ADF test
[h,pValue] = adftest(abs_neg)

[h,pValue] = adftest(pos)
