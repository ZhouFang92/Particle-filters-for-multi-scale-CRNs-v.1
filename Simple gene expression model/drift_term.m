function a=drift_term(t,x) %x1 corresponds to X_3 (RNA) and x2 corresponds to protein.

k4=0.390;

k6=0.379;
a=[0;0; 0;k4*x(3)-k6*x(4)];