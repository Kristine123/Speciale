%%Lånt fra github
%https://github.com/cgroll/mff/blob/master/matlab_for_finance_files/matlab_script_3.m

%% Simulation basis

n = 1000;     % sample size
nBins = 20;   % number of points per bin

% simulate from uniform distribution
simU = rand(n, 1);

% check appropriateness
hist(simU, nBins)

% calculate expected number of bins
expNum = n/nBins;

% include horizontal line at expected number of points per bin
line([0 1], expNum*[1 1], 'Color', 'r')
title(['sample size: ' num2str(n) ' / expected points per bin: '...
    num2str(expNum)])

%% PIT example

% init params
n2 = 10;    % sample size
mu = 4;     % params distribution 
sigma = 1;

% show first n2 points of simU on y-axis
close
scatter(zeros(n2, 1), simU(1:n2), 'r.')
hold on;

% plot cumulative distribution function
grid = 0:0.01:8;
vals = normcdf(grid, mu, sigma);
plot(grid, vals)

% scatter resulting values after inverse PIT
X = norminv(simU, mu, sigma);
scatter(X(1:n2), zeros(n2, 1), 'r.')

% plot lines between points to indicate relation
for ii=1:n2
    % horizontal line, beginning in uniform sample point
    line([0 X(ii)], [simU(ii) simU(ii)], 'Color', 'r',...
       'LineStyle', ':')
    
    % vertical line, beginning in new point, up to function
    line([X(ii) X(ii)], [0 simU(ii)], 'Color', 'r', 'LineStyle', ':')
end
title('Inverse probability integral transformation')