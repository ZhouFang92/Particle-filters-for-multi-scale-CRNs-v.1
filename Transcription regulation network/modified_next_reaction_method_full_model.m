% This function use modified next reaction method to simulate the full
% model. It returns state vector X and time vector TV, given final time FT,
% initial conditions X0.

%S3-> S1+S3; S1->0; S5->S5+S3; S3->0
%S2+S4->S5; S5->S2+S4; 2S1->S2; S2->2S1

function [K,X,TV]=modified_next_reaction_method_full_model(K,X0,FT)

%initialization (k is the scaled reaction constant)
N=100;
%k1=K(1);
%k2=K(2);
%k3=K(3);
%k4=K(4);
%k5=K(5);
%k6=K(6);
%k7=K(7);
%k8=K(8);
X(:,1)=X0;
TV(1)=0;
for j=1:8
  T(j)=0;
  r(j)=rand;
  a(j)=propensity(X(:,1),j,K,N);
  P(j)=log(1/r(j));
  DeltaT(j)=(P(j)-T(j))/a(j);
end
[deltaT,mu] = min(DeltaT);
num_step=1;


while TV(num_step)+deltaT < FT % update while the final time is not reached.
     TV(num_step+1)=TV(num_step)+deltaT; %update time
     DX=reaction_direction(mu); 
     DX=[1/N;1/N;1;1;1].*DX;
     X(:,num_step+1)=X(:,num_step)+DX; %update the state
     if X(1,num_step+1)<0              % avoid negative concentration
        X(1,num_step+1)=0;
     end    
     if X(2,num_step+1)<0              % avoid negative concentration
        X(2,num_step+1)=0;
     end   
     num_step=num_step+1;
     for j=1:8
         T(j)=T(j)+deltaT*a(j);       %update the used time for each reaction
     end
     r(mu)=rand;
     P(mu)=P(mu)+log(1/r(mu));         %update the next reaction time for mu
     for j=1:8                   %recalculate the propensity
           a(j)=propensity(X(:,num_step),j,K,N);
           DeltaT(j)=(P(j)-T(j))/a(j);
     end
     [deltaT,mu] = min(DeltaT);  %update the next time to fire a reaction
     %TV(num_step)
end


TV(num_step+1)=FT; %update the final time
X(:,num_step+1)=X(:,num_step);


