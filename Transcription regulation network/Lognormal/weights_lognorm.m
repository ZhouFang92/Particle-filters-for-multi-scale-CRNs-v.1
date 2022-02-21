% update weights where w is the previous weights, x is 
% the simulated state and y is observation. 
% the update is modified according to the lognorm distribution

function w=weights_lognorm(w,x,y) 

w=w*lognpdf(y-h_function(x),0,1/2);