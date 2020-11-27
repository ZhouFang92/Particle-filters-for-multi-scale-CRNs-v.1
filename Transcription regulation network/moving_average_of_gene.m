function av=moving_average_of_gene(T,X,t) 
% T is the time vector, X is the state vector, and t is the time interval
% to calculate the average.

tmax=T(size(T,2));

j=1;
for i=1:size(T,2)
    if T(i)<=t
       av(i)=mean(X(5,1:i));
    else
       av(i)=av(i-1)*(i-j)+X(5,i);
       while T(j)<T(i)-t
           j=j+1;
           av(i)=av(i)-X(5,j);
       end
       av(i)=av(i)/(i-j+1);
       %av(i)=mean(X(5,j:i));
    end
    %T(i)
end

