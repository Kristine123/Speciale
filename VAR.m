%% VAR

for i=1:length(lambda_f) %her har du tilføjet -1 fordi?
    c(i)= 1/poisscdf(10,lambda_f(i)); %change for LESS than 10 hospitilazation. For MORE than 1 - xxx
end

