% The particle filter using full model, where the inputs are the initial
% condition of the system X0 obervations (TY, Y) and the number of particles M.
% The outputs are the filters (T_filter,filter) and the estimation 
% of the standard deviation.

function [T_filter,filter, filterSD, parameters]=particle_filter_reduced_model_discrete_time(TY, Y, M)

TY=[0,TY];
Y=[0,Y];
T_filter=TY;

% initialize particles
parfor j=1:M
    [k(:,j),v(:,j)]=system_parameters;
    w(j)=1/M;
end

filter=zeros(10,size(TY,2));
filterSD=zeros(10,size(TY,2));

    % calculate filter at the initial time
    i=1;
    filterTemp=filter(:,i);
    filterSDTemp=filterSD(:,i);
    sumW=sum(w);
    parfor j=1:M 
       filterTemp=filterTemp+w(j)*[k(:,j);v(:,j)]/sumW;
       filterSDTemp=filterSDTemp+w(j)*[k(:,j);v(:,j)].^2/sumW;
    end
    filter(:,i)=filterTemp;
    filterSD(:,i)=filterSDTemp-filterTemp.^2;

parameters.diversity=zeros(size(TY,2),1); 
parameters.diversity(i)=numel(unique(k(1,:)));
    

for i=2:size(TY,2)
   
    %sampling set
      parfor j=1:M
        [K,XF1,TXF1]=modified_next_reaction_method_limit_model(k(:,j),v(:,j),TY(i)-TY(i-1));
        v(:,j)=XF1(:,size(XF1,2));
        k(:,j)=K;
      end
    
    %update weights
    parfor j=1:M
        w(j)=weights_calculation_discrete_time(w(j),v(:,j),Y(i));
    end
    
    %compute filter
    filterTemp=filter(:,i);
    filterSDTemp=filterSD(:,i);
    sumW=sum(w);
    parfor j=1:M 
       filterTemp=filterTemp+w(j)*[k(:,j);v(:,j)]/sumW;
       filterSDTemp=filterSDTemp+w(j)*[k(:,j);v(:,j)].^2/sumW;
    end
    filter(:,i)=filterTemp;
    filterSD(:,i)=filterSDTemp-filterTemp.^2;
    
    % compute the diversity
    parameters.diversity(i)=numel(unique(k(1,:)));
    if i==size(TY,2)
       parameters.particles= kv(1:6,:);
    end
    
    %resampling
    [kv,w]=resampling_step([k;v],w);
    k=kv(1:6,:);
    v=kv(7:10,:);
    
    %TY(i)
    
end


