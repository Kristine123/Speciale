%% Forsøg med OLS
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');     
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_18.txt'); % Ikke transformerede. Trans er i temp_16
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_18.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');

%%
y = hosp;
X = temp;
fitlm(X,y)

%%
X = vacc;
fitlm(X,y)
%%
X = tested;
fitlm(X,y)
%%
X = tested;
fitlm(X,y)
%%
X=abs(min(0,diff(temp)));
T0 = 464;
% T0 is set to lenght of time series -1 when including the exogenes variables. This is necessary as taking first difference of exogenes variables results
% in a shorter time series

y = hosp';
y = y(:,1:T0); 

%%
X=abs(min(0,diff(temp)));
T0 = 464;
% T0 is set to lenght of time series -1 when including the exogenes variables. This is necessary as taking first difference of exogenes variables results
% in a shorter time series

y = hosp';
y = y(:,1:T0); 

fitlm(X,y)
%%
y = hosp;
X = strin;
fitlm(X,y)
