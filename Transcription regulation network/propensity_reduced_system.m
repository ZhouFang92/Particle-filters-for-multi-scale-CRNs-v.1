function a=propensity_reduced_system(X,j,K)

if j==3
    a=K(3)*K(5)*psi_function(X,K)*(X(4)+X(5))/(K(5)*psi_function(X,K)+K(6));
end

if j==4
   a=K(4)*X(3);
end