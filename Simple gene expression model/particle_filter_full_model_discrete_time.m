% The particle filter using full model, where the inputs are the initial
% condition of the system X0 obervations (TY, Y) and the number of particles M.
% The outputs are the filters (T_filter,filter) and the estimation 
% of the standard deviation.

function [T_filter,filter, filterSD]=particle_filter_full_model_discrete_time(TY, Y, M)

TY=[0,TY];
Y=[0,Y];
T_filter=TY;
filter=zeros(10,size(TY,2));
filterSD=zeros(10,size(TY,2));



% initialize particles

for j=1:M
    [k(:,j),v(:,j)]=system_parameters;
    w(j)=1/M;
end


diversity(1)=numel(unique(k(1,:)));

    % calculate filter at the initial time
i=1;
    for j=1:M 
       filter(:,i)=filter(:,i)+[k(:,j);v(:,j)]/M;
       filterSD(:,i)=filterSD(:,i)+[k(:,j);v(:,j)].^2/M;
    end
    filterSD(:,i)=filterSD(:,i)-(filter(:,i)).^2;


for i=2:size(TY,2)
   
    %sampling set
       parfor j=1:M
         [K,XF1,TXF1]=modified_next_reaction_method_full_model(k(:,j),v(:,j),TY(i)-TY(i-1));
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
    
    %resampling
    [kv,w]=resampling_step([k;v],w);
    k=kv(1:6,:);
    v=kv(7:10,:);
    
    %TY(i)
    
    %diversity(i)=numel(unique(k(1,:)));
    
end


