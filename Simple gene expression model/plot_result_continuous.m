function plot_result_continuous(TY,Y,K,TXF,XF,T_filter,filter_true,filter_true2,filter_approximate,filter_approximate2)

N=100;

T(1)=TXF(1);
X(:,1)=XF(:,1);
for i=2:size(TXF,2)
   T(2*i-3)=TXF(i);
   X(:,2*i-3)=XF(:,i-1);
   T(2*i-2)=TXF(i);
   X(:,2*i-2)=XF(:,i);
end
TXF=T;
XF=X;

if size(K,1)==1
   K=K';
end
  

%for i=1:size(T_filter,2)
%   t=find( TXF <= T_filter(i), 1, 'last' );
%   XF1(:,i)= XF(:,t);
%   TXF1(i) = T_filter(i);
%end
%XF=XF1;
%TXF=TXF1;

subplot(4,2,1);
plot(TY,Y,'black','LineWidth',2)
xlabel('Time (minutes)','FontSize',16)
ylabel('Observations','FontSize',16)
title('(A)','FontSize',16)

subplot(4,2,2);
X=[1:6];
K_error=[abs(filter_true(1:6,size(filter_true,2))-K)./K, abs(filter_approximate(1:6,size(filter_true,2))-K)./K];
SD_error=[abs(filter_true2(1:6,size(filter_true,2)))./K, abs(filter_approximate2(1:6,size(filter_true,2)))./K];
b=bar(X,K_error,0.9);
axis([0.5 6.5 0 0.3]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'$k''_1$','$k''_2$','$k''_3$','$k''_4$','$k''_5$','$k''_6$'},'FontSize',13);
set(b(1),'FaceColor','r');
set(b(2),'FaceColor','b');
legend(b,{'Particle filters (original models)','Particle filters (reduced models)'},'Location','northeast','FontSize',12);
ylabel({'RE of parameters'; 'at the final time'},'FontSize',16)
xlabel('Reaction constants','FontSize',16)
title('(B)','FontSize',16)

%hold on;
%x=zeros(1,12);
%y=x;
%sd=x;
%for i=1:6
%  x(2*i-1)= i-0.15;
%  x(2*i)= i+0.15;
%end
%for i=1:6
%  y(2*i-1)= K_error(i,1);
%  y(2*i)= K_error(i,2);
%end
%for i=1:6
%  sd(2*i-1)=SD_error(i,1);
%  sd(2*i)=SD_error(i,2);
%end
%errorbar(x,y,sd);
%v=cell2mat(get(cell2mat(get(b,'child')),'vertices')');
%errorbar(X,K_error,K_error-SD_error,K_error+SD_error);

%subplot(3,2,3);

%for i=1:size(filter_true,2)
%error_K(i)=norm(filter_true(1:6,i)-filter_approximate(1:6,i));
%mean_K(i)= norm(filter_true(1:6,i));
%error_K(:,i)=filter_true(1:6,i)-K;
%to_plot(:,i)=abs(error_K(:,i)./K);
%end

%for i=1:6
%    h(i)=plot(T_filter,to_plot(i,:),'--*','LineWidth',1);
%    hold on;
%end
%legend(h,{'$k''_1$','$k''_2$','$k''_3$','$k''_4$','$k''_5$','$k''_6$'},'Interpreter','Latex','Location','northwest','FontSize',12);

%axis([0 T_filter(size(T_filter,2)) 0 max(max(to_plot))]);
%xlabel('Time (seconds)','FontSize',16)
%ylabel({'Relative errors in'; 'infered parameters'},'FontSize',16)
%title('(A)','position',[25 -0.25],'FontSize',16);

%subplot(3,2,3);

%for i=1:size(filter_true,2)
%error_K(i)=norm(filter_true(1:6,i)-filter_approximate(1:6,i));
%mean_K(i)= norm(filter_true(1:6,i));
%error_K(:,i)=filter_approximate(1:6,i)-K;
%to_plot(:,i)=abs(error_K(:,i)./K);
%end

%for i=1:6
%    h(i)=plot(T_filter,to_plot(i,:),'--*','LineWidth',1);
%    hold on;
%end
%legend(h,{'$k''_1$','$k''_2$','$k''_3$','$k''_4$','$k''_5$','$k''_6$'},'Interpreter','Latex','Location','northwest','FontSize',12);

%axis([0 T_filter(size(T_filter,2)) 0 max(max(to_plot))]);
%xlabel('Time (seconds)','FontSize',16)
%ylabel({'Relative errors in'; 'infered parameters'},'FontSize',16)
%title('(A)','position',[25 -0.25],'FontSize',16);


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
true_value=plot(TXF,XF(2,:),'black','LineWidth',2);
hold on;
bayesian_estimation=plot(T_filter, filter_true(8,:),'r','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, filter_approximate(8,:),'b','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 1+0.1])
yticks([0:0.2:1])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',16)
ylabel('Activated DNA','FontSize',16)
%title('(B)','position',[25 -0.25],'FontSize',16);
legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Particle filters (reduced models)'},'Location','south','FontSize',12);
title('(C)','FontSize',16)

subplot(4,2,4);
bayesian_estimation=plot(T_filter, sqrt(filter_true2(8,:)),'r','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, sqrt(filter_approximate2(8,:)),'b','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 1+0.1])
yticks([0:0.2:1])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',16)
ylabel({'SD of activated DNA'},'FontSize',16)
%title('(B)','position',[25 -0.25],'FontSize',16);
legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Particle filters (reduced models)'},'Location','southwest','FontSize',12);
title('(D)','FontSize',16)


subplot(4,2,5);
true_value=plot(TXF,XF(3,:),'black','LineWidth',2);
hold on;
bayesian_estimation=plot(T_filter, filter_true(9,:),'r','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, filter_approximate(9,:),'b','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 fix(max(XF(4,:)))+3])
yticks([0:1:fix(max(XF(3,:))+3)])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',16)
ylabel('mRNA','FontSize',16)
%title('(C)','position',[25 -6.5/4],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Paticle filters (reduced models)'},'Location','northeast','FontSize',12);
print -dtiff 'x.tif'
title('(E)','FontSize',16)

subplot(4,2,6);
bayesian_estimation=plot(T_filter, sqrt(filter_true2(9,:)),'r','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, sqrt(filter_approximate2(9,:)),'b','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 fix(max(XF(4,:)))+3])
yticks([0:1:fix(max(XF(3,:))+3)])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',16)
ylabel({'SD of mRNA'},'FontSize',16)
%title('(B)','position',[25 -0.25],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Particle filters (reduced models)'},'Location','southwest','FontSize',12);
title('(F)','FontSize',16)

subplot(4,2,7);
true_value=plot(TXF,N*XF(4,:),'black','LineWidth',2);
hold on
bayesian_estimation=plot(T_filter, N*filter_true(10,:),'r','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, N*filter_approximate(10,:),'b','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 (fix(max(XF(4,:)))+1.5)*N])
yticks([0:100:(fix(max(XF(4,:)))+1.5)*N])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',16)
ylabel('Fluorescent proteins','FontSize',16)
ax=gca;
ax.YAxis.Exponent = 2;
%title('(D)','position',[25 -500/4],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filter (original model)','Paticle filter (reduced model)'},'Location','northeast','FontSize',12);
title('(G)','FontSize',16)

subplot(4,2,8);
bayesian_estimation=plot(T_filter, N*sqrt(filter_true2(10,:)),'r','LineWidth',1);
hold on
approximate_filter_plot=plot(T_filter, N*sqrt(filter_approximate2(10,:)),'b','LineWidth',1);
hold on
axis([0 T_filter(size(T_filter,2)) 0 (fix(max(XF(4,:)))+1.5)*N])
yticks([0:100:(fix(max(XF(4,:)))+1.5)*N])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',16)
ylabel({'SD of FP'},'FontSize',16)
%title('(B)','position',[25 -0.25],'FontSize',16);
ax=gca;
ax.YAxis.Exponent = 2;
title('(H)','FontSize',16)
