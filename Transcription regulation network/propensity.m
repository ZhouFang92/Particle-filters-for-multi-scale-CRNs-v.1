% return the value of propensity given the state X, the number of
% reaction j, and reaction constants kj.

function a=propensity(X,j,K,N) 

if j == 1 
    a=K(1)*X(3)*N;
end

if j == 2
    a=K(2)*X(1)*N;
end

if j==3
    a=K(3)*X(5);
    return
end

if j==4
    a=K(4)*X(3);
end

if j==5
    a=K(5)*X(2)*X(4)*N;
end

if j==6
    a=K(6)*X(5)*N;
end

if j==7
    a=abs(K(7)*N*X(1)*(N*X(1)-1));
end

if j==8
    a=K(8)*X(2)*N*N;
end
