function [T_filter,filter, filterSD]= particle_filter_full_model_continuous_time(TY, Y, M, delta)

%SD is actually the variance.

resampling_period= ceil (delta/(TY(2)-TY(1))); %resampling time
dTY=TY(2)-TY(1);  % time grid
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
    
    
for iteration=1 : floor(size(TY,2)/resampling_period)
   
    %sampling
    for j=1:M
         [K,XF1,TXF1]=modified_next_reaction_method_full_model(k(:,j),v(:,j),resampling_period*dTY);
         v(:,j)=XF1(:,size(XF1,2));
         k(:,j)=K;
         v_temp(:,:,j)=extract_continuous_path (XF1,TXF1,TY(2)-TY(1));
    end
    
    
    %calculate weights
    for j=1:M
      w_temp(:,j)=calculate_weights_continuous_time(v_temp(:,:,j),Y,TY,(i-1)*dTY);
      w(j)=w_temp(size(w_temp,1),j);
    end
    
    
    
    %calculate filter
    l=i;
    for i=l+1:l+resampling_period
       for j=1:M 
          filter(:,i)=filter(:,i)+w_temp(i-l,j)*[k(:,j);v_temp(:,i-l,j)]/sum(w_temp(i-l,:));
          filterSD(:,i)=filterSD(:,i)+w_temp(i-l,j)*[k(:,j);v_temp(:,i-l,j)].^2/sum(w_temp(i-l,:));
       end
       filterSD(:,i)=filterSD(:,i)-(filter(:,i)).^2;
    end
    

    
    %resample
    [kv,w1]=resampling_step([k;v],w);
    w=w1;
    k=kv(1:6,:);
    v=kv(7:10,:);
    
   % TY(i)
   
   %diversity(i)=numel(unique(k(1,:))); 
    
end
