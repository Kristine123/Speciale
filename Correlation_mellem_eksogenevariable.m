%%
clc
clear
%% Data 
hosp=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');     
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_18.txt'); % Ikke transformerede. Trans er i temp_16
contact=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_18.txt');
vacc=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_18.txt');
strin=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
varbrit=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/BritiskVariant.txt');
vardelta=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/DeltaVariant.txt');

%% Matrix to be used when calculation correlation values
X(:,1)=hosp;
X(:,2)=temp;
X(:,3)=contact;
X(:,4)=tested;
X(:,5)=vacc;
X(:,6)=varbrit;
X(:,7)=vardelta;
X(:,8)=strin;

%% Matrix to be used when calculating VIF values
% X(:,1)=temp;
% X(:,2)=contact;
% X(:,3)=tested;
% X(:,4)=vacc;
% X(:,5)=varbrit;
% X(:,6)=vardelta;

%% Transformed dataseries first difference
X(:,1)=abs(min(0,diff(temp)));
X(:,2)=max(0,diff(temp));

X(:,3)=abs(min(0,diff(contact)));
X(:,4)=max(0,diff(contact));

X(:,5)=abs(min(0,diff(tested))); 
X(:,6)=max(0,diff(tested));

X(:,7)=abs(diff(vacc)); 
X(:,8)=abs(min(0,diff(varbrit)));

X(:,9)=max(0,diff(varbrit));
X(:,10)=abs(min(0,diff(vardelta)));

X(:,11)=max(0,diff(vardelta));
X(:,12)=hops;
%% Covarians matrix
C = cov(X)

%% Correlation matrix
[R,P] = corrcoef(X)
%xlswrite('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/correlation.xlsx',R)
%% Correlation plot
figure
corrplot(X,'testR','on', 'varNames',{'Hospitaslity','Temperature','Contact','Tested','Vaccinated','British','Delta','Restrictions'})
%% VIF values
VIF = diag(inv(R))'