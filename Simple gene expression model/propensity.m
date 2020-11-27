% return the value of propensity given the state X, the number of
% reaction j, and reaction constants kj.

function a=propensity(X,j,k1,k2,k3,k4,k5,k6) 

if j == 1 
    a=k1*X(1);
end

if j == 2
    a=k2*X(2);
end

if j==3
    a=k3*X(2);
end

if j==4
    a=k4*X(3);
end

if j==5
    a=k5*X(3);
end

if j==6
    a=k6*X(4);
end

