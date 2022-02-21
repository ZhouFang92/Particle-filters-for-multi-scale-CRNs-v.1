function [t,x]=ode_Euler(K,x0,ti,tf) %initial condition 

delta=0.05;
size=(tf-ti)/delta;

t(1)=ti;
x(:,1)=x0;
x01=x0(1);
x02=x0(2);
x03=x0(3);
x04=x0(4);
x05=x0(5);
i=1;

if size>=2 
for i=2:size
   dx1=K(1)/5*x(3,i-1)-K(2)/5*(5*x(1,i-1)-2*psi_function(x(:,i-1),K));
   t(i)=t(i-1)+delta;
   x(1,i)=x(1,i-1)+dx1*delta;
   x(2,i)=x(2,i-1)+2*dx1*delta;
   x(3,i)=x03;
   x(4,i)=x04;
   x(5,i)=x05;
end
else 
    i=1;
end

   t(i+1)=tf;
   dx1=K(1)/3*x(3,i)-K(2)/3*(5*x(1,i)-2*psi_function(x(:,i),K));
   x(1,i+1)=x(1,i)+dx1*(tf-t(i));
   x(2,i+1)=x(2,i)+2*dx1*(tf-t(i));
   x(3,i+1)=x03;
   x(4,i+1)=x04;
   x(5,i+1)=x05;
