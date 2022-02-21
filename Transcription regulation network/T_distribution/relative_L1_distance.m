function distance=relative_L1_distance(filter1,filter1SD,filter2,filter2SD)

% This function compute the L2 distance between two filters

% protein

mean1=filter1(9,:)+2*filter1(10,:);
mean2=filter2(9,:)+2*filter2(10,:);
delta=mean1-mean2;
distance(1,1)=sum(abs(delta))/(sum(mean1));

variance1=filter1SD(9,:);
variance1(1)=filter1SD(9,1)+4*filter1SD(10,1);
variance2=filter2SD(9,:)+4*filter2SD(10,:)+4*sqrt(filter2SD(9,:).*filter2SD(10,:));
variance2(1)=filter2SD(9,1)+4*filter2SD(10,1);
delta=sqrt(variance1)-sqrt(variance2);
distance(1,2)=sum(abs(delta))/sum(sqrt(variance1));

% mRNA

mean1=filter1(11,:);
mean2=filter2(11,:);
delta=mean1-mean2;
distance(2,1)=sum(abs(delta))/(sum(mean1));

variance1=filter1SD(11,:);
variance2=filter2SD(11,:);
delta=sqrt(variance1)-sqrt(variance2);
distance(2,2)=sum(abs(delta))/sum(sqrt(variance1));

% Activated Gene

mean1=filter1(13,:);
mean2=filter2(13,:);
delta=mean1-mean2;
distance(3,1)=sum(abs(delta))/(sum(mean1));

variance1=filter1SD(13,:);
variance2=filter2SD(13,:);
delta=sqrt(variance1)-sqrt(variance2);
distance(3,2)=sum(abs(delta))/sum(sqrt(variance1));

