function [t,x]=ode_Euler(k3,k4,x0,ti,tf) %initial condition 

delta=0.05;
size=(tf-ti)/delta;

t(1)=ti;
x(:,1)=x0;
x01=x0(1);
x02=x0(2);
x03=x0(3);
i=1;

if size>=2 
for i=2:size
   t(i)=t(i-1)+delta;
   x(1,i)=x01;
   x(2,i)=x02;
   x(3,i)=x03;
   x(4,i)=x(4,i-1)+ delta*(k3*x03-k4*x(4,i-1));
end
else 
    i=1;
end

   t(i+1)=tf;
   x(1,i+1)=x01;
   x(2,i+1)=x02;
   x(3,i+1)=x03;
   x(4,i+1)=x(4,i)+ (tf-t(i))*(k3*x03-k4*x(4,i));
