%% Rolling estimation and forecasting of US defaults using PARX(2,1) model
%% with X = [RV,LI]
clear all


%% Load data
hosp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/hosp_16.txt');      
temp = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/temp_16.txt'); 
contact = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/contact_16.txt');
tested = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/Tested_pcroganti_16.txt');
vacc = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/vacc_not_16.txt');
strin = importdata('/Users/krmmm/Documents/Dokumenter_Mac/Speciale/Model_data/stringency_16.txt');
%%

Y=hosp';
X=[];
% X(1,:)=temp;
% X(2,:)=contact;
% X(3,:)=tested;
% X(4,:)=vacc;

%OBS y og X er her transponeret ifht den originale kode for at passe ind i
%output_parx_rob funktionen

%%
%Leading_minus = -Leading_index.*(Leading_index < 0);
 
%Y = def';
%X = [rv';Leading_minus'];
 
[dimX,T] = size(Y);
%Date = linspace(1982,2011,T);
t1 = datetime(2020,03,26,8,0,0);
t2 = datetime(2021,07,3,8,0,0);
Date = (t1:t2);
%% Initial values of theta used in numerical optimization computed
%theta_init(:,1) = [0.4421;0.1119;0.4296;0.2903;0.0000;0.2890];
%theta_init(:,2) = [0.7994;0.2152;0.4341;zeros(3,1)];
%theta_init(:,3) = [0.0000; 0.1302; 0.0693; 0.7337; 56.9762; 2.3939];
%theta_init(:,4) = [zeros(3,1); 0.8168; 99.2496; 0.6913];
%theta_init(:,5) = [0.2284; 0.1842; 0.1990; 0.4827; 23.3633; 0.8029];
%k = 6;
%for i = 1:5
%   for j = i+1:5
%        theta_init(:,k) = (theta_init(:,i) + theta_init(:,j))/2;
%        k = k+1;
%    end
%end
 
%[dimtheta,N_init] = size(theta_init);

 
%theta_est = zeros(8,T);
%logL      = zeros(1,T);
%logL_tmp  = zeros(1,N_init);
%theta_sqrt_est = zeros(8,N_init);
 
%Hertil giver det mening hvad der bliver lavet, pÂ en besvÊrlig mÂde laver
%den matrice med forskellige mulige start-vÊrdier, k¯rer modeller med dem
%alle sammen og gemmer den med lavest loglikelihood. Det der nedenfor er
%lidt mÊrkeligt for der overskrives? de f¯rste 5 rÊkker med identiske
%rÊkker, bortset fra den ¯verste, der er defineret mÊrkeligt.
% theta_init(:,1) = [0.4914    0.1331    0.4606    0.2092    0.0000    0.3822    0.543     0.321]';
% theta_init(:,2) = [0.7994;0.2152;0.4341;zeros(5,1)];
% theta_init(:,3) = [0.0000; 0.1302; 0.0693; 0.7337; 56.9762; 2.3939; 2.352; 0.462];
% theta_init(:,4) = [zeros(3,1); 0.8168; 99.2496; 0.6913; 0.420; 3.234];
% theta_init(:,5) = [0.2284; 0.1842; 0.1990; 0.4827; 23.3633; 0.8029; 0.363; 2.534];

% theta_init(:,1) = [0.4914;  0.1331; 0.4606];
% theta_init(:,2) = [0.7994;  0.2152; 0.4341];
% theta_init(:,3) = [0.0000;  0.1302; 0.0693];
% theta_init(:,4) = [0;       0;      0];
% theta_init(:,5) = [0.2284;  0.1842; 0.1990];

% Set Theta_init
theta_init = [0.2284;  0.1842; 0.1990];

% Determine size of theta_init
[dimtheta,N_init] = size(theta_init);

% Create matrices to put values into
theta_est = zeros(3,T);
logL      = zeros(1,T);
% logL_tmp  = zeros(1, N_init);
% theta_sqrt_est = zeros(3, N_init);


%% Initialization of algorithm
%TOLERANCE ER MEGET LAV FOR AT TESTE AT KODEN KØRER. TolFun samt TolX skal
%sættes op inden koden køres 'rigtigt'. Default værdien er 1e-4. 
options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-1,'TolX', 1e-1, 'Display', 'iter','PlotFcns',@optimplotfval);
T0 = 365; % Set number of days to include

[theta_sqr_est, logL_tmp] = fminsearch('logL_PARX', sqrt(theta_init), options, Y(1:T0),X,1,1);

% for i = 1:N_init
%     [theta_sqr_est(:,i),logL_tmp(i)] = fminsearch('logL_PARX',sqrt(theta_init(:,i)),options,Y(1:T0),X,1,1);
% end
% 

%NEDENSTÅENDE VERSION ER GAMMEL OG BRUGES TIL NÅR DER ER EKSOGENE VARIABLE
%INKLUDERET. FORSKELLEN LIGGER I X(:,1:T0). Denne giver når vi ikke har
%eksogene variable med fejlkoden: Index in position 2 exceeds array bounds.
%ENDVIDERE ER STARTVÆRDI NOGET LORT OG LIGEGYLDIGT.

% for i = 1:N_init
%     [theta_sqr_est(:,i),logL_tmp(i)] = fminsearch('logL_PARX',sqrt(theta_init(:,i)),options,Y(1:T0),X(:,1:T0),1,1);
% end


 %logL_PARX(theta,Y,X,p,q)
% [logLmin,index] = min(logL_tmp);
% [theta_sqr_est(:,index),logL(T0)] = fminsearch('logL_PARX',theta_sqr_est(:,index),options,Y(1:T0),X,1,1);
%[theta_sqr_est(:,index),logL(T0)] = fminsearch('logL_PARX',theta_sqr_est(:,index),options,Y(1:T0),X(:,1:T0),1,1);
theta_est(:,T0) = theta_sqr_est.^2; 
theta_init(:,N_init+1) = theta_est(:,T0);
theta_init = theta_init';
N_init = N_init + 1;
 


% [logLmin,index] = min(logL_tmp);
% [theta_sqr_est(:,index),logL(T0)] = fminsearch('logL_PARX',theta_sqr_est(:,index),[],Y(1:T0),X(:,1:T0),2,1);
% theta_est(:,T0) = theta_sqr_est(:,index).^2; 
% theta_init(:,N_init+1) = theta_est(:,T0);
% N_init = N_init + 1;
% 
% %forecast of lambda_t for t = T0
% lambda_f(T0+1) = PARX_forecast(theta_est(:,T0),Y(1:T0),X(:,1:T0),2,1); 



%Den laver f¯rst en "normal estimation" derefter propper den de vÊrdier
%ind i en forecast model og forecaster en periode ud i fremtiden... Det g¯r
%den sÂ en helt masse gange, sÂ vi kan evaluere performance
 
%forecast of lambda_t for t = T0

%lambda_f(T0+1) = PARX_forecast(theta_est(:,T0),Y(1:T0),X(:,1:T0),1,1); 
lambda_f(T0+1) = PARX_forecast(theta_est(:,T0),Y(1:T0),X,1,1); 

%lambda_forecast = PARX_forecast(theta,Y,X,p,q)
 
%% Iterative estimation and forecasting
% We re-estimate the four models adding one additional observation at a
% time; these are then used to compute one-step ahead forecast of lambda_t
T0 = 300
options = optimset('Algorithm','interior-point','MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-8,'TolX', 1e-8, 'Display', 'iter','PlotFcns',@optimplotfval);
disp('time period')
%Den skriver bare time period
 
for t = T0+1:T
    disp(t)
    
    for i = 1:N_init
        %X(:,1:t) ÆNDRET TIL X GRUNDET INGEN EKSOGENE
        [theta_sqr_est(:,i),logL_tmp(i)] = fminsearch('logL_PARX',sqrt(theta_init(:,i)),options,Y(1:t),X,1,1);
    end
    
    [logLmin,index] = min(logL_tmp);
    %X(:,1:t) ÆNDRET TIL X GRUNDET INGEN EKSOGENE
    [theta_sqr_est(:,index),logL(t)] = fminsearch('logL_PARX',theta_sqr_est(:,index),options,Y(1:t),X,1,1);
    theta_est(:,t) = theta_sqr_est(:,index).^2;
    theta_init(:,N_init) = theta_est(:,t);

 %X(:,1:t) ÆNDRET TIL X GRUNDET INGEN EKSOGENE
    lambda_f(t+1) = PARX_forecast(theta_est(:,t),Y(1:t),X,1,1);  %forecast
 
    %if (compu == 'UCL   ')
       % save 'U:\Fag\Programmel\MatLab\PARX\Results\theta_est_iter.mat' theta_est
       % save 'U:\Fag\Programmel\MatLab\PARX\Results\lambda_f.mat' lambda_f
   % elseif (compu == 'laptop')
        %save 'C:\Users\dennis\Documents\Fag\Programmel\MatLab\PARX\Results\theta_est_iter.mat' theta_est
        %save 'C:\Users\dennis\Documents\Fag\Programmel\MatLab\PARX\Results\lambda_f.mat' lambda_f
   % end
   
   %De to sidste her er ligegyldige, blot til at vise periode.
   disp(i)
   disp(t)
end 
 
save '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Outputs/PARX/theta_est_iter.mat' theta_est 
save '/Users/krmmm/Documents/Dokumenter_Mac/MATLAB/PARX_1/Outputs/PARX/lambda_f.mat' lambda_f
 
ForecastError(T0+1:T) = Y(T0+1:T)-lambda_f(T0+1:T);
logf(T0+1:T) = Y(T0+1:T).*log(lambda_f(T0+1:T)) - lambda_f(T0+1:T);
 
KLIC = zeros(1,T);
MSFE = zeros(1,T);
 
for t = T0+1:T
    MSFE(t) = sum(ForecastError(T0+1:t).^2)/(t-T0);
    KLIC(t) = -sum(logf(T0+1:t))/(t-T0);
end
 
figure; plotyy(Date(T0+1:T),MSFE(T0+1:T),Date(T0+1:T),KLIC(T0+1:T))
title('Out-of-sample fit: MSE and Score evaluation')
legend('MSFE','KLIC')
 
figure; plot(Date(T0+1:T),theta_est(1,T0+1:T))
title('Parameter stability')
legend('\omega')

figure; plot(Date(T0+1:T),theta_est(2,T0+1:T))
title('Parameter stability')
legend('\alpha_1','\alpha_2')

figure; plot(Date(T0+1:T),theta_est(3,T0+1:T))
title('Parameter stability')
legend('\beta')

% figure; plotyy(Date(T0+1:T),theta_est(5,T0+1:T),Date(T0+1:T),theta_est(6,T0+1:T))
% title('Parameter stability')
% legend('TEMP','Kontakt tal')
% 
% figure; plotyy(Date(T0+1:T),theta_est(7,T0+1:T),Date(T0+1:T),theta_est(8,T0+1:T))
% title('Parameter stability')
% legend('Test','VAX')