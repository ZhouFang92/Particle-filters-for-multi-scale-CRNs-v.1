function plot_result(TY,Y,K,TXF,XF,av,T_filter,filter_true,filter_true2,filter_approximate,filter_approximate2)

N=100;
K=K';
TY=TY*100/60;
TXF=TXF*100/60;
T_filter=T_filter*100/60;

%Resampling the data
  %resampling_size=150;
  XF_final=XF(:,size(TXF,2));
  TXF_final=TXF(:,size(TXF,2));
  av_final=av(:,size(av,2));
  XF=XF(:,1:fix(size(TXF,2)/200):end);
  XF(:,size(XF,2)+1)=XF_final;
  av=av(:,1:fix(size(TXF,2)/200):end);
  av(:,size(av,2)+1)=av_final;
%  filter_true=filter_true(:,1:fix(size(T_filter,2)/resampling_size):end);
%  filter_true2=filter_true2(:,1:fix(size(T_filter,2)/resampling_size):end);
%  filter_approximate=filter_approximate(:,1:fix(size(T_filter,2)/resampling_size):end);
%  filter_approximate2=filter_approximate2(:,1:fix(size(T_filter,2)/resampling_size):end);
  TXF=TXF(1:fix(size(TXF,2)/200):end);
%  T_filter=T_filter(1:fix(size(T_filter,2)/resampling_size):end);
  TXF(:,size(TXF,2)+1)=TXF_final;
  

T(1)=TXF(1);
X(:,1)=XF(:,1);
for i=2:size(TXF,2)
   T(2*i-3)=TXF(i);
   X(:,2*i-3)=XF(:,i-1);
   av1(:,2*i-3)=av(:,i);
   T(2*i-2)=TXF(i);
   X(:,2*i-2)=XF(:,i);
   av1(:,2*i-2)=av(:,i);
end
TXF=T;
XF=X;
av=av1;

  

%for i=1:size(T_filter,2)
%   t=find( TXF <= T_filter(i), 1, 'last' );
%   XF1(:,i)= XF(:,t);
%   TXF1(i) = T_filter(i);
%end
%XF=XF1;
%TXF=TXF1;

subplot(4,2,1);

plot(TY,Y,'black--*')
xlabel('Time (minutes)','FontSize',16)
ylabel('Observations','FontSize',16)
title('(A)','FontSize',16)
xticks([0:30:150])

subplot(4,2,2)
X=[1:8];
K_error=[abs(filter_true(1:8,size(filter_true,2))-K)./K, abs(filter_approximate(1:8,size(filter_true,2))-K)./K];
SD_error=[abs(filter_true2(1:8,size(filter_true,2)))./K, abs(filter_approximate2(1:8,size(filter_true,2)))./K];
b=bar(X,K_error,0.9);
axis([0.5 8.5 0 0.4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'$k''_1$','$k''_2$','$k''_3$','$k''_4$','$k''_5$','$k''_6$','$k''_7$','$k''_8$'},'FontSize',13);
set(b(1),'FaceColor','r');
set(b(2),'FaceColor','b');
legend(b,{'Particle filters (original models)','Particle filters (reduced models)'},'Location','northeast','FontSize',12);
ylabel({'RE of parameters'; 'at the final time'},'FontSize',16)
xlabel('Reaction constants','FontSize',16)
title('(B)','FontSize',16)




%low_error_bar=max(filter_true(7,:)-sqrt(filter_true2(7,:)),0);
%upper_error_bar=min(filter_true(7,:)+sqrt(filter_true2(7,:)),1);
%error_bar_matrix=[(low_error_bar)', (upper_error_bar-low_error_bar)'];
%h=area(T_filter,error_bar_matrix,-1,'LineStyle',':');
%h(1).FaceColor = [1 1 1];
%h(2).FaceColor = [1 0 0];
%h(2).FaceAlpha = 0.1;
%hold on
%true_value=plot(TXF,XF(1,:),'black','LineWidth',1);
%hold on;
%bayesian_estimation=plot(T_filter, filter_true(7,:),'r--*','LineWidth',1);
%hold on
%approximate_filter_plot=plot(T_filter, filter_approximate(7,:),'b--*','LineWidth',1);
%hold on
%axis([0 T_filter(size(T_filter,2)) 0 1])
%set(gca,'FontSize',12);
%yticks([0:0.2:1])
%box off;     
%xlabel('Time (seconds)','FontSize',16)
%ylabel('Inactivated DNA','FontSize',16)
%title('(A)','position',[25 -0.25],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filter (original model)','Paticle filter (reduced model)'},'Location','northwest','FontSize',12);
%legend('Exact value','True filter','approximate filter')


subplot(4,2,3);
true_value=plot(TXF,N*XF(1,:)+2*N*XF(2,:),'black','LineWidth',2);
hold on;
bayesian_estimation=plot(T_filter, N*filter_true(9,:)+2*N*filter_true(10,:),'r--*','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, N*filter_approximate(9,:)+2*N*filter_approximate(10,:),'b--*','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 19*N 75*N])
yticks([20*N:10*N:70*N])
xticks([0:30:150])
set(gca,'FontSize',12); 
xlabel('Time (minutes)','FontSize',16)
ylabel('Total mass of proteins','FontSize',16)
title('(C)','FontSize',16)
legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Particle filters (reduced models)'},'Location','northwest','FontSize',12);
ax=gca;
ax.YAxis.Exponent = 2;

subplot(4,2,4);
bayesian_estimation=plot(T_filter, sqrt(N^2*filter_true2(9,:)+2^2*N^2*filter_true2(10,:)),'r--*','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, sqrt(N^2*filter_approximate2(9,:)+2*N^2*sqrt(filter_approximate2(9,:).*filter_approximate2(10,:))+2^2*N^2*filter_approximate2(10,:)),'b--*','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0*N 10*N])
yticks([0*N:2*N:10*N])
xticks([0:30:150])
set(gca,'FontSize',12); 
xlabel('Time (minutes)','FontSize',16)
ylabel({'SD of the total' ; 'mass of proteins'},'FontSize',16)
%title('(B)','position',[25 -0.25],'FontSize',16);
legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Particle filters (reduced models)'},'Location','northwest','FontSize',12);
ax=gca;
ax.YAxis.Exponent = 2;
title('(D)','FontSize',16)


subplot(4,2,5);
true_value=plot(TXF,XF(3,:),'black','LineWidth',2);
hold on;
bayesian_estimation=plot(T_filter, filter_true(11,:),'r--*','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, filter_approximate(11,:),'b--*','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 fix(max(XF(3,:)))+1])
yticks([0:1:fix(max(XF(3,:))+3)])
xticks([0:30:150])
set(gca,'FontSize',12);
xlabel('Time (minutes)','FontSize',16)
ylabel('mRNA','FontSize',16)
%title('(C)','position',[25 -6.5/4],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Paticle filters (reduced models)'},'Location','northeast','FontSize',12);
title('(E)','FontSize',16)

subplot(4,2,6);
bayesian_estimation=plot(T_filter, sqrt(filter_true2(11,:)),'r--*','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, sqrt(filter_approximate2(11,:)),'b--*','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 fix(max(XF(3,:)))+1])
yticks([0:1:fix(max(XF(3,:))+3)])
xticks([0:30:150])
set(gca,'FontSize',12);
xlabel('Time (minutes)','FontSize',16)
ylabel('SD of mRNA','FontSize',16)
title('(F)','FontSize',16)

subplot(4,2,7);
true_value=plot(TXF,av,'black','LineWidth',2);
hold on;
bayesian_estimation=plot(T_filter, filter_true(13,:),'r--*','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, filter_approximate(13,:),'b--*','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0.5 1.09])
yticks([0:0.1:1])
xticks([0:30:150])
set(gca,'FontSize',12);
xlabel('Time (minutes)','FontSize',16)
ylabel('Bounded DNA','FontSize',16)
%title('(C)','position',[25 -6.5/4],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Paticle filters (reduced models)'},'Location','northeast','FontSize',12);
title('(G)','FontSize',16)

subplot(4,2,8);
bayesian_estimation=plot(T_filter, sqrt(filter_true2(13,:)),'r--*','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, sqrt(filter_approximate2(13,:)),'b--*','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 0.5])
yticks([0:0.1:1])
xticks([0:30:150])
set(gca,'FontSize',12);
xlabel('Time (minutes)','FontSize',16)
ylabel('SD of bounded DNA','FontSize',16)
title('(H)','FontSize',16)