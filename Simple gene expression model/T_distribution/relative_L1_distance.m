function distance=relative_L1_distance(filter1,filter1SD,filter2,filter2SD)

% This function compute the L2 distance between two filters

for i=1:10
    delta=filter1(i,:)-filter2(i,:);
    distance(i,1)=sum(abs(delta));%sqrt(sum(((filter1(i,:)-filter2(i,:))./filter1(i,:)).^2)); 
    distance(i,1)=distance(i,1)/(sum(filter1(i,:)));
end

for i=1:10
    delta=sqrt(filter1SD(i,:))-sqrt(filter2SD(i,:));
    distance(i,2)=sum(abs(delta));%sqrt(sum(((filter1(i,:)-filter2(i,:))./filter1(i,:)).^2)); 
    distance(i,2)=distance(i,2)/sum(sqrt(filter1SD(i,:)));
end
