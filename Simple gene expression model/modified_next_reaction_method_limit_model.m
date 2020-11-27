% This function use modified next reaction method to simulate the limit
% model. Specifically, I apply the modified next reaction method to 
% simulate the dynamics of gene and RNA, and simulate the dynamic of 
% fluorescent protein by ODE solver. It returns state vector X and time 
% vector TV, given final time FT,initial conditions X0.

% The reaction network: DNA-->DNA*  DNA*-->DNA  DNA*-->DNA*+mRNA
% mRNA-->RNA+FP RNA-->0  FP-->0


function [K,X,TV]=modified_next_reaction_method_limit_model(K,X0,FT)

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
  if j== 4 | j==6
      DeltaT(j)=FT;
  else
       T(j)=0;
       r(j)=rand;
       a(j)=propensity(X(:,1),j,k1,k2,k3,k4,k5,k6);
       P(j)=log(1/r(j));
       DeltaT(j)=(P(j)-T(j))/a(j);
  end   
end
[deltaT,mu] = min(DeltaT);
num_step=1;
opts = odeset('RelTol',1e-1);



while TV(num_step)+deltaT < FT % update while the final time is not reached.
    %[microT,x]=ode45('drift_term',[TV(num_step),TV(num_step)+deltaT],X(:,num_step)); 
    %microT=microT';
    %x=x';
    [microT,x]=ode_Euler(k4,k6,X(:,num_step),TV(num_step),TV(num_step)+deltaT); %update the protein molecule
    TV=[TV,microT];    %update time
    X=[X,x];
    num_step=size(TV,2);
    DX=reaction_direction(mu); 
    X(:,num_step)=X(:,num_step)+DX; % update DNA and protein
    if X(4,num_step)<0
        X(4,num_step)=0;
    end  
    for j=1:6
         if j== 4 | j==6
         else
           T(j)=T(j)+deltaT*a(j);       %update the used time for each reaction
         end
     end
     r(mu)=rand;
     P(mu)=P(mu)+log(1/r(mu));         %update the next reaction time for mu
     for j=1:6                   %recalculate the propensity and the next reaction time
           if j== 4 | j==6
              DeltaT(j)=FT;
           else
              a(j)=propensity(X(:,num_step),j,k1,k2,k3,k4,k5,k6);
              DeltaT(j)=(P(j)-T(j))/a(j);
           end   
     end
     [deltaT,mu] = min(DeltaT);  %update the next time to fire a reaction
end

%update the final time
%[microT,x]=ode45('drift_term',[TV(num_step),FT],X(:,num_step));
%microT=microT';
%x=x';
[microT,x]=ode_Euler(k4,k6,X(:,num_step),TV(num_step),FT);
TV=[TV,microT];    %update time
X=[X,x];


