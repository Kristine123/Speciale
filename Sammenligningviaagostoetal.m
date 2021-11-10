%Forsøg på at estimere deres data med begge tilsendte koder og se hvad der passer med paperet

%%
%Load data
%def count, rv, credit_spread
def=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/PARX_codes/PARX_code_October2014/Code&Script/defcountUS.txt');
rv=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/PARX_codes/PARX_code_October2014/Code&Script/rv.txt');
%spread=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/PARX_codes/PARX_code_October2014/Code&Script/credit_spread.txt');
li=importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/PARX_codes/PARX_code_October2014/Code&Script/Leading_index.txt');

%%
%Format til model gammel kode
y=def;
X=rv;
X(:,2)=li;
%%
%%
%Gamle kode model 

%Model with two covariates:

theta0=[0.1;0.2;0.2;0.2;0.2;0.2];
lb=-Inf*ones(6,1);
opt=optimset('Algorithm','interior-point');

p = 2; % autoreg led
q = 1; % intensity lagged
s = 1; % Excel sheet 
%k = p+q+1;

filename = 'Outputs/agostoetal30.xlsx';

out1 = output(y,X,theta0,lb,opt,p,q,s,filename); 
%output([10 12 15],   X,  theta0,lb,opt,p,q,s,filename);
%output(data,cov,theta0,lb,opt,p,q,s,filename)


%%
%Format til ny kode:
y=def';
X(1,:)=rv;
X(2,:)=li;

%%
%Nye kode model
%PARX(2,1)

p = 2; % autoreg led
q = 1; % intensity lagged
s = 5; % Excel sheet 
%g = p+q+5;

%theta0 = 0.1*ones(g,1);
theta0=[0.15095;0.18053;0.1098;0.4098;0.011;0.5]; 
%den er ikke stabil overfor ændringer i intitial arbitraty values
% opt=optimset('Algorithm','interior-point');
%Denne er kørt med to forskellige indstillinger i options omkring
%algorimte

filename = 'Outputs/PARX/agosto_30.xlsx';

%[OUT]=output_parx_rob(y,X,theta0,lb,opt,p,q,s,filename)
[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  
%[OUT, pearson_residuals, pred, zero_prob, meany, rpit] = output_parx_rob(y,X,theta0,p,q,s,filename);  

%% Predicted defaults per month
%t1 = datetime(2020,03,26,8,0,0);
%t2 = datetime(2021,07,3,8,0,0);
%Date = (t1:t2);
%Date = linspace(366,1,1)

Date= linspace(1,359,359)';
new=def(2:360);

figure; plot(Date,pred,Date,new)
%xlim(datetime([2020 2021],[03 08],1));
%title('Predictions vs Fitted');
legend('Predictions','True');


%% Residual Autocorrelation of pearson residuals
%Do you need to understand the confidence intervals? probaly do to also
%apply them to predictions

figure; autocorr(pearson_residuals)

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
%t1 = datetime(2020,03,26,8,0,0);
%t2 = datetime(2021,07,3,8,0,0);
%Date = (t1:t2);
Date= linspace(1,358,358)';
yyaxis left
plot(Date,zero_prob)
%xlim(datetime([2020 2021],[03 08],1))
z =(def==0);

yyaxis right
plot(x(z),def(z),'linestyle','none','marker','*')
ylim([0,10])
%axis off
set(gca,'ytick',[])
title('Zeros and corresponding estimated probabilities')
legend('Estimated probability of zero count','Empirical zero counts')
%saveas(gcf,'/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Figurer/zero_prob.jpg')


%% PIT
%Du har konstrueret rpit inde i selve rob filen
%kolmogorov-smirnov
[h,p] = kstest(rpit)

%et histogram kunne være fedt
figure
%t1 = datetime(2020,03,26,8,0,0);
%t2 = datetime(2021,07,3,8,0,0);
%Date = (t1:t2);
Date= linspace(1,360,360)';
u = Date;
nbins=25;
histogram(rpit,nbins)

%Bincounts
%Counts = h.Values

%Bin width
%h = histogram(C,'BarWidth',0.5)
%Giver måske bedre mening at specificere width end antal søjler
