clc
clear


%%
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');

%%      
%temp=importdata('/Users/krmmm/Documents/Dokumenter - MacBook Pro/Speciale/Model_data/Temperatur_366.txt');
%X=temp;  
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_16.txt'); 
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_16.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_16.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
%%

y=hosp';
%X(1,:)=temp;
X(1,:)=temp;
X(2,:)=contact;
X(3,:)=tested;
X(4,:)=vacc;

%OBS y og X er her transponeret ifht den originale kode for at passe ind i
%output_parx_rob funktionen


%%
%PARX(1,1)
%Benchmark model

p = 2; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[0.2;0.1;0.2;0.2;0.54;0.01;0.81;0.81;0.81;0.81]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');

filename = 'Outputs/PARX_22/forsogingenrestrektion.xlsx';

%function[OUT,pears,pred]=output_parx_rob(y,X,theta0,p,q,s,filename)  
[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  


%% Predicted defaults per month
t1 = datetime(2020,03,28,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);

%Creating the subsample length, removing observations at the start of the
%period, responding to the lag length
%k is lag length + 1, to start from the first observation with sufficient
%lags. g is the subsample length, the endpoint 465 is quite a "manual"
%soluition, and could benefit from being calculated by code determining the
%length of the full sample. 
k = p+1;
g=y(k:465);
figure; plot(Date,pred,Date,g)
xlim(datetime([2020 2021],[03 08],1));
title('Predictions vs Fitted');
legend('Predictions','True');
%saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/predictions.jpg')

%% Residual Autocorrelation of pearson residuals
%Do you need to understand the confidence intervals? probaly do to also
%apply them to predictions

figure; autocorr(pearson_residuals)
%saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/ACF_residuals.jpg')
%% Ljungbox test - is it okay that you are using pearson residuals? yes why not?
%e
[h,pValue] = lbqtest(pearson_residuals,'lags',[7,14,21])
%e^2 - ARCH effects
[h,pValue] = lbqtest(pearson_residuals.^2,'lags',[7,14,21])

%% Zero counts and probabilities 
%Here you have used matlabs cdfpoisson specification and put x=0, yup plot
%needs some love

%figure; plot(Date,zero_prob)
figure
t1 = datetime(2020,03,28,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
x = Date;
y = zero_prob;
yyaxis left
plot(x,y)
xlim(datetime([2020 2021],[03 08],1))

z =(hosp==0);

yyaxis right
plot(x(z),hosp(z),'linestyle','none','marker','*')
ylim([0,10])
%axis off
set(gca,'ytick',[])
title('Zeros and corresponding estimated probabilities')
legend('Estimated probability of zero count','Empirical zero counts')
%saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/zero_prob.jpg')


%% PIT
%Du har konstrueret rpit inde i selve rob filen
%kolmogorov-smirnov
[h,p] = kstest(rpit)

%et histogram kunne være fedt
figure
t1 = datetime(2020,03,28,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
x = Date;
y = rpit;
nbins=25;
histogram(y,nbins)

%Bincounts
%Counts = h.Values

%Bin width
%h = histogram(C,'BarWidth',0.5)
%Giver måske bedre mening at specificere width end antal søjler
saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/Midlertidige_Figurer/pit_histogram.jpg')


