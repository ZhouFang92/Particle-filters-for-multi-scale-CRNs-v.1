function [T_filter,filter, filterSD,parameters]= particle_filter_full_model_continuous_time(TY, Y, M, delta)

resampling_period= ceil (delta/(TY(2)-TY(1))); %resampling time
dTY=TY(2)-TY(1);  % time grid
T_filter=TY;
filter=zeros(13,size(TY,2));
filterSD=zeros(13,size(TY,2));


% initialize particles

parfor j=1:M
    [k(:,j),v(:,j)]=system_parameters;
    w(j)=1/M;
end

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
    
    
for iteration=1 : floor(size(TY,2)/resampling_period)
   
    %sampling
    parfor j=1:M
         [K,XF1,TXF1]=modified_next_reaction_method_full_model(k(:,j),v(:,j),resampling_period*dTY);
         v(:,j)=XF1(:,size(XF1,2));
         k(:,j)=K;
         v_temp(:,:,j)=extract_continuous_path (XF1,TXF1,TY(2)-TY(1));
    end
    
    
    %calculate weights
    parfor j=1:M
      w_temp_temp=calculate_weights_continuous_time(v_temp(:,:,j),Y,TY,(i-1)*dTY);
      w_temp_temp=w_temp_temp';
      w_temp(:,j)=w_temp_temp;
      w(j)=w_temp_temp(size(w_temp_temp,1));
    end
    
    
    
    %calculate filter
    %calculate filter
    l=i;
    for i=l+1:l+resampling_period
           filterTemp=filter(:,i);
           filterSDTemp=filterSD(:,i);
           filterSDTemp_P=0;
           sumW=sum(w_temp(i-l,:));
           parfor j=1:M 
              filterTemp=filterTemp+w_temp(i-l,j)*[k(:,j);v_temp(:,i-l,j)]/sumW;
              filterSDTemp=filterSDTemp+w_temp(i-l,j)*[k(:,j);v_temp(:,i-l,j)].^2/sumW;
              filterSDTemp_P=filterSDTemp_P+w_temp(i-l,j)*(v_temp(1,i-l,j)+2*v_temp(2,i-l,j))^2/sumW;
           end
           filter(:,i)=filterTemp;
           filterSD(:,i)=filterSDTemp-filterTemp.^2;
           filterSD(9,i)=filterSDTemp_P-(filterTemp(9)+2*filterTemp(10))^2;
    end
    

    % compute the diversity
    parameters.diversity(i)=numel(unique(k(1,:)));
    if i==size(TY,2)
       parameters.particles= kv(1:6,:);
    end   
    
    %resample
    [kv,w1]=resampling_step([k;v],w);
    w=w1;
    k=kv(1:8,:);
    v=kv(9:13,:);
    
    TY(i)
   
      
    
end
