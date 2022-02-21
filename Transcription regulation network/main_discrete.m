% this function computes the result in discrete-time observation case.

load matlab.mat
M=10000;

%% Particle filters based on full models 
% (NOTICE: It takes more than 20 hours to output the result. Please run it 
%  on a cluster server. NEVER run it on your PC!)

% Discrete-time observation
tic
 [T_filter,filter_full_model, filterSD_full_model, parameters_DF]=particle_filter_full_model_discrete_time(TY, Y, M);
time_df=toc;


%% Particle filters based on limit models

%Discrete-time observations
tic
  [T_filter,filter_reduced_model, filterSD_reduced_model, parameters_DR]=particle_filter_reduced_model_discrete_time(TY, Y, M);
time_dr=toc;


%% save
save discrete