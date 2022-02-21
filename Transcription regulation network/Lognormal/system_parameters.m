function [K,X0]=system_parameters

% This function return system parameters sampled from the given
% distribution

K(1)=0.4+0.1*rand;
K(2)=0.6+0.3*rand;
K(3)=0.7+0.2*rand;
K(4)=0.3+0.2*rand;
K(5)=1+2*rand;
K(6)=0.3+0.2*rand;
K(7)=6+3*rand;
K(8)=0.4+0.2*rand;

if rand<1/10
    X0(4)= 1;
else
    X0(4)= 0;
end
X0(5)= 1-X0(4);
X0(3)= poissrnd(2);
X0(2)= poissrnd(10);
X0(1)= poissrnd(1);

%X0=[1;10;2;0;1];
