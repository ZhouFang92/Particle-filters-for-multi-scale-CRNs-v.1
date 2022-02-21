% update weights where w is the previous weights, x is 
% the simulated state and y is observation. 
% the update is modified according to the Laplace distribution

function w=weights_t_distribution(w,x,y) 

w=w*tpdf(abs(y-h_function(x)),4);