%%Discription: in this example, we consider a simple gene expression system
%%that has four species inactivated gene (S1), activated gene (S2), and the
%%mRNA (S3), and the protein (S4).There are four reactions involved:
%%S1-->S2, S2-->S1, S2-->S2+S3, S3-->S3+S4, S3-->0, S4-->0. The gene
%%product can be measured by a fluorescent microscope, where
%%h(x)=10*min(x(4),100).
%%Now, we consider the observation noise to have the Laplace distribution
%%instead of the Guassian distribution

%% Generate a reference process and observation
M=100000; %number of particles
sampling_time_observation=2; % The sampling period of discrete-time observations
X0=[0;1;2;2];   %True initial conditions
K(1)=0.014;     %The model parameters of the system
K(2)=0.0084;
K(3)=0.715;
K(4)=0.390;
K(5)=0.199;
K(6)=0.379;
% To randomly generate the initial conditions and model parameters
%[K,X0]=system_parameters;
K=K';
FT=90;      %Time length

[K,XF,TXF]=modified_next_reaction_method_full_model(K,X0,FT); % trajectory
[Y,TY]=observation_discrete_time(TXF,XF,sampling_time_observation); % Observation

%% Particle filter based on full model
% Discrete-time observation
tic
[T_filter,filter_full_model, filterSD_full_model, parameters_DF]=particle_filter_full_model_discrete_time(TY, Y, M);
time_df=toc   

%% Perticle filter limit model
% Discrete-time observation
tic
 [T_filter,filter_reduced_model, filterSD_reduced_model,parameters_DR]=particle_filter_reduced_model_discrete_time(TY, Y, M);
time_dr=toc

%% Plot the result

% Discrete-time 
plot_result(TY,Y,K,TXF,XF,T_filter,filter_full_model, filterSD_full_model,filter_reduced_model, filterSD_reduced_model);


%% Relative L2 distance

distance_discrete=relative_L1_distance(filter_full_model,filterSD_full_model,filter_reduced_model,filterSD_reduced_model);

%% save

% save double_exponential_noise