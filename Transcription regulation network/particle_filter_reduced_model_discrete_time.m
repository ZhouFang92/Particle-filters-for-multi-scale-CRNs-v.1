% The particle filter using full model, where the inputs are the initial
% condition of the system X0 obervations (TY, Y) and the number of particles M.
% The outputs are the filters (T_filter,filter) and the estimation 
% of the standard deviation.

function [T_filter,filter, filterSD]=particle_filter_reduced_model_discrete_time(TY, Y, M)

TY=[0,TY];
Y=[0,Y];
T_filter=TY;

filter=zeros(13,size(TY,2));
filterSD=zeros(13,size(TY,2));


% initialize particles

for j=1:M
    [k(:,j),v(:,j)]=system_parameters;
end

    % calculate filter at the initial time
    i=1;
    for j=1:M 
       filter(:,i)=filter(:,i)+[k(:,j);v(:,j)]/M;
       filterSD(:,i)=filterSD(:,i)+[k(:,j);v(:,j)].^2/M;
    end
    filterSD(:,i)=filterSD(:,i)-(filter(:,i)).^2;

    Pi=[1/5 2/5 0 0 0; 2/5 4/5 0 0 0; 0 0 1 0 0; 0 0 0 1/2 1/2; 0 0 0 1/2 1/2];

for j=1:M    
    v(:,j)=Pi*v(:,j);
    w(j)=1/M;
end

for i=2:size(TY,2)
   
    %sampling set
      for j=1:M
        [K,XF1,TXF1]=modified_next_reaction_method_limit_model(k(:,j),v(:,j),TY(i)-TY(i-1));
        v(:,j)=XF1(:,size(XF1,2));
        k(:,j)=K;
      end
    
    %update weights
    for j=1:M
        w(j)=weights_calculation_discrete_time(w(j),v(:,j),Y(i));
    end
    
    %compute filter
    for j=1:M 
       filter(:,i)=filter(:,i)+w(j)*[k(:,j);v(:,j)]/sum(w);
       filterSD(:,i)=filterSD(:,i)+w(j)*[k(:,j);v(:,j)].^2/sum(w);
    end
    filterSD(:,i)=filterSD(:,i)-(filter(:,i)).^2;
    filterSD(13,i)=filter(13,i)-(filter(13,i)).^2;
    
    %resampling
    [kv,w]=resampling_step([k;v],w);
    k=kv(1:8,:);
    v=kv(9:13,:);
    
end