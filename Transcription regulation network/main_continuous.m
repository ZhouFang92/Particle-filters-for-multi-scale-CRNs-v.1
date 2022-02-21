% this function computes the result in continuous-time observation case.

load matlab.mat
M=10000;

%% Particle filters based on full models 
% (NOTICE: It takes more than 20 hours to output the result. Please run it 
%  on a cluster server. NEVER run it on your PC!)

% Continuous-time observation
tic
 [T_filter_continuous,filter_continuous_FM, filterSD_continuous_FM,parameters_CF]= particle_filter_full_model_continuous_time(TYC, YC, M, 0.5);
time_cf=toc;

%% Particle filters based on limit models

%Continuous-time observations
tic
 [T_filter_continuous,filter_continuous_RM, filterSD_continuous_RM, parameters_CR]= particle_filter_reduced_model_continuous_time(TYC, YC, M, 0.5);
time_cr=toc;

%% save

save continuous