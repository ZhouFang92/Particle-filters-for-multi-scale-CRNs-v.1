function [Y,TY]=observation_continuous_time(TVF,XF)

delta=0.25;

TY=[0:delta:TVF(size(TVF,2))];

Y(1)=0;

for i=2:size(TY,2)
    tmax=find(TVF <= delta*i,1,'last');       % Find the last update before the end of interval
    Y(i)=Y(i-1)+h_function(XF(:,tmax))*delta+normrnd(0,delta);
end