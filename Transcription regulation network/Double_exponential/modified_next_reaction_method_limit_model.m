% This function use modified next reaction method to simulate the limit
% model. Specifically, I apply the modified next reaction method to 
% simulate the dynamics of gene and RNA, and simulate the dynamic of 
% fluorescent protein by ODE solver. It returns state vector X and time 
% vector TV, given final time FT,initial conditions X0.

% The reaction network: DNA-->DNA*  DNA*-->DNA  DNA*-->DNA*+mRNA
% mRNA-->RNA+FP RNA-->0  FP-->0


function [K,X,TV]=modified_next_reaction_method_limit_model(K,X0,FT)

delta=0.05;
N=100;
Pi=[1/5 2/5 0 0 0; 2/5 4/5 0 0 0; 0 0 1 0 0; 0 0 0 1/2 1/2; 0 0 0 1/2 1/2];
X(:,1)=Pi*X0;
TV(1)=0;
DeltaT=ones(8,1)/0;
for j=3:4
  T(j)=0;
  r(j)=rand;
  a(j)=propensity_reduced_system(X(:,1),j,K);
  P(j)=log(1/r(j));
  DeltaT(j)=(P(j)-T(j))/a(j);
end
[deltaT,mu] = min(DeltaT);
num_step=1;



while TV(num_step)+delta <= FT % update while the final time is not reached.
    
    TV(num_step+1)=TV(num_step)+delta;
    num_step=num_step+1;
    i=num_step;
    
    %update continuous process
    dx1=K(1)/5*X(3,i-1)-K(2)/5*(5*X(1,i-1)-2*psi_function(X(:,i-1),K));
    X(1,i)=X(1,i-1)+dx1*delta;
    X(2,i)=X(2,i-1)+2*dx1*delta;
    X(3,i)=X(3,i-1);
    X(5,i)=K(5)*psi_function(X(:,i),K)*(X(4,i-1)+X(5,i-1))/(K(5)*psi_function(X(:,i),K)+K(6));%X(5,i-1);
    X(4,i)=X(4,i-1)+X(5,i-1)-X(5,i);
    %propensity_3=2*K(3)*K(5)*psi_function(X(:,i-1),K)*X(4,i-1)/(K(5)*psi_function(X(:,i-1),K)+K(6));
    T(3)=T(3)+propensity_reduced_system(X(:,i),3,K)*delta;
    T(4)=T(4)+propensity_reduced_system(X(:,i),4,K)*delta;
    
    %update discrete process
    
    if T(3)>P(3)
       X(3,i)=X(3,i)+1;
       r(3)=rand;
       P(3)=P(3)+log(1/r(3));
    end
    
    if T(4)>P(4)
       X(3,i)=X(3,i)-1;
       r(4)=rand;
       P(4)=P(4)+log(1/r(4));
    end
    
end

