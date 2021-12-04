function GMPT_pipeline
addpath('../GMPT_Common_function')

N = 30; % community size
diff_C = 0.05; % Connectivity of ecological network
delta = 0.2; % variance of aij
diag = -1; % the value of aii
VarianceType = 2; % 1: Travis; 2: Variance
exchange_rate = [0.1 0.3 0.5]; % phenotype = 2 * length(exchange_rate) + 2
Num_sample = 20;
simulation = 5;

Many_times_different_sample_size(N,diff_C,delta,diag,VarianceType,exchange_rate,Num_sample,simulation)

end