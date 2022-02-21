function [K,X0]=system_parameters

% This function return system parameters sampled from the given
% distribution

K(1)=0.01+0.01*rand;
K(2)=0.007+0.003*rand;
K(3)=0.7+0.2*rand;
K(4)=0.3+0.1*rand;
K(5)=0.1+0.2*rand;
K(6)=0.3+0.1*rand;

if rand<1/3
    X0(1)= 0;
else
    X0(1)= 1;
end
X0(2)= 1-X0(1);
X0(3)= poissrnd(2);
X0(4)= poissrnd(2);
