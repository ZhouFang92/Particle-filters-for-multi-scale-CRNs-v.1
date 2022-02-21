function w_temp=calculate_weights_continuous_time(X_temp,Y,TY,initial_time)

base=find(TY <= initial_time,1,'last');  %the index before the iteration
dTY=TY(2)-TY(1);

w_temp(1)= exp( -(Y(base+1)-Y(base))^2+h_function(X_temp(:,1))*(Y(base+1)-Y(base))-(h_function(X_temp(:,1)))^2/2*dTY );

for i=2:size(X_temp,2)
    w_temp(i)= w_temp(i-1)*exp( -(Y(base+1)-Y(base))^2+h_function(X_temp(:,i))*(Y(base+i)-Y(base+i-1))-(h_function(X_temp(:,i)))^2/2*dTY );
end
