function X_temp= extract_continuous_path (X,TX,delta)

for i=1: floor(TX(size(TX,2))/delta)
    tmax=find(TX <= delta*i,1,'last'); 
    X_temp(:,i)= X(:,tmax);
end