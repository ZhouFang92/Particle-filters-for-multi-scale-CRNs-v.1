% This function use modified next reaction method to simulate the full
% model. It returns state vector X and time vector TV, given final time FT,
% initial conditions X0.

% The reaction network: DNA-->DNA*  DNA*-->DNA  DNA*-->DNA*+mRNA
% mRNA-->RNA+FP RNA-->0  FP-->0


function [K,X,TV]=modified_next_reaction_method_full_model(K,X0,FT)

%initialization (k is the scaled reaction constant)
N=100;
%k1=0.0140;%0.700;
%k2=0.0084;%0.430;
%k3=0.715;
%k4=0.390;
%k5=0.199;
%k6=0.379;
k1=K(1);
k2=K(2);
k3=K(3);
k4=K(4);
k5=K(5);
k6=K(6);
X(:,1)=X0;
TV(1)=0;
for j=1:6
  T(j)=0;
  r(j)=rand;
  a(j)=propensity(X(:,1),j,k1,k2,k3,k4,k5,k6);
  if j== 4 | j==6
      a(j)=a(j)*N;
  end   
  P(j)=log(1/r(j));
  DeltaT(j)=(P(j)-T(j))/a(j);
end
[deltaT,mu] = min(DeltaT);
num_step=1;


while TV(num_step)+deltaT < FT % update while the final time is not reached.
     TV(num_step+1)=TV(num_step)+deltaT; %update time
     DX=reaction_direction(mu); 
     if mu== 4 | mu==6
         DX=DX/N;
     end
     X(:,num_step+1)=X(:,num_step)+DX; %update the state
     if X(4,num_step+1)<0              % avoid negative concentration
        X(4,num_step+1)=0;
     end    
     num_step=num_step+1;
     for j=1:6
         T(j)=T(j)+deltaT*a(j);       %update the used time for each reaction
     end
     r(mu)=rand;
     P(mu)=P(mu)+log(1/r(mu));         %update the next reaction time for mu
     for j=1:6                   %recalculate the propensity
           a(j)=propensity(X(:,num_step),j,k1,k2,k3,k4,k5,k6);
           if j== 4 | j==6
              a(j)=a(j)*N;
           end   
           DeltaT(j)=(P(j)-T(j))/a(j);
     end
     [deltaT,mu] = min(DeltaT);  %update the next time to fire a reaction
end


TV(num_step+1)=FT; %update the final time
X(:,num_step+1)=X(:,num_step);


