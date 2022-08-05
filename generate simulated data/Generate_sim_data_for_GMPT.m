function Generate_sim_data_for_GMPT
N = 30;
diff_C = 0.08; 
delta = 0.2;
diag = -1;
noise = 0.05;
shuffle_SF_Amatrix = 0;
VarianceType = 2; % 1: Travis; 2: Variance
exchange_rate = [0.1 0.3 0.5]; % phenotype = 2 * length(exchange_rate) + 2
Num_sample = 20;
simulation = 5;
dynamical_parameter = 'connectivity';

for i = 1 : length(diff_C)
    Many_times_different_sample_size_v2(dynamical_parameter,N,diff_C(i),delta,diag,VarianceType,exchange_rate,Num_sample,simulation,noise,shuffle_SF_Amatrix);
end

end