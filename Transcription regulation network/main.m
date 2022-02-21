%Discription: in this example, we consider a transcription regulation
%network that has five species: the protein monomer (S1), the protein dimer, 
%the mRNA (S3), the inactivated DNA (S4), and the activated DNA (S5).
% reaction scheme: 
% S3-> S1+S3; S1->0; S5->S5+S3; S3->0
% S2+S4->S5; S5->S2+S4; 2S1->S2; S2->2S1
% Observation: h(x)=min(x(1)+2*x(2),1000);

%% initialization
%clear;
M=10000; %number of particles
sampling_time_observation=2; % resampling frequnecy
X0=[1;10;2;0;1];   %Initial conditions
K(1)=0.45;         %Model parameters
K(2)=0.75; 
K(3)=0.8;
K(4)=0.390;
K(5)=2;
K(6)=0.379;
K(7)=7.5;
K(8)=0.5;
% To randomly generate the initial conditions and model parameters
%[K,X0]=system_parameters;
K=K';
FT=90;      %Total time length
delta=0.5;  %The period of the resampling time

%% Generate a reference process and observation (It may take several minutes)
 tic
 [K,XF,TXF]=modified_next_reaction_method_full_model(K,X0,FT);
 av=moving_average_of_gene(TXF,XF,1); 
 toc
 [Y,TY]=observation_discrete_time(TXF,XF,sampling_time_observation);
 [YC,TYC]=observation_continuous_time(TXF,XF);

%% Particle filters based on full models 
% (NOTICE: It takes more than 10 hours to output the result. Please run it 
%  on a cluster server. NEVER run it on your PC!)

% Discrete-time observation
tic
 [T_filter,filter_full_model, filterSD_full_model, parameters_DF]=particle_filter_full_model_discrete_time(TY, Y, M);
time_df=toc;

% Continuous-time observation
tic
 [T_filter_continuous,filter_continuous_FM, filterSD_continuous_FM,parameters_CF]= particle_filter_full_model_continuous_time(TYC, YC, M, 0.5);
time_cf=toc;

%% Particle filters based on limit models

%Discrete-time observations
tic
[T_filter,filter_reduced_model, filterSD_reduced_model, parameters_DR]=particle_filter_reduced_model_discrete_time(TY, Y, M);
time_dr=toc;

%Continuous-time observations
tic
 [T_filter_continuous,filter_continuous_RM, filterSD_continuous_RM, parameters_CR]= particle_filter_reduced_model_continuous_time(TYC, YC, M, 0.5);
time_cr=toc;


%% Plot the result

%Discrete-time observations
plot_result(TY,Y,K,TXF,XF,av,T_filter,filter_full_model, filterSD_full_model,filter_reduced_model, filterSD_reduced_model);

%Continuous-time observations
plot_result_continuous(TYC,YC,K,TXF,XF,av,T_filter_continuous,filter_continuous_FM, filterSD_continuous_FM,filter_continuous_RM, filterSD_continuous_RM);

%% plot particle distribution

plot_particle_distribution(K,parameters_CF,parameters_CR,parameters_DF,parameters_DR)


%% Plot the CPU time

CPU_time(time_cf,time_cr, time_df, time_dr);

%% Distance

distance_discrete=relative_L1_distance(filter_full_model, filterSD_full_model,filter_reduced_model, filterSD_reduced_model);
distance_continuous=relative_L1_distance(filter_continuous_FM, filterSD_continuous_FM,filter_continuous_RM, filterSD_continuous_RM);
