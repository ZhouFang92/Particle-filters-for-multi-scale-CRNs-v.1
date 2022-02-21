function plot_result(TY,Y,K,TXF,XF,T_filter,filter_true,filter_true2,filter_approximate,filter_approximate2)

% initialization

f = figure;
f.Units='pixels';
f.OuterPosition=[10 10 800 750];

font_size=12;
title_size=15;
marker_size=15;
lengend_font_size=9;

% Basic settings

N=100;

if size(K,1)==1
   K=K';
end
  



% observations
subplot(4,2,1);
plot(TY,Y,'k.--','MarkerSize',marker_size)
xlabel('Time (minutes)','FontSize',font_size)
ylabel('Observations','FontSize',font_size)
title('(A)','FontSize',title_size)
axis([0 90 0 60])
xticks([0:10:90])
yticks([0:20:60])

subplot(4,2,2);
X=[1:6];
Relative_K=[abs(filter_true(1:6,size(filter_true,2)))./K, abs(filter_approximate(1:6,size(filter_true,2)))./K];
SD_error=[abs(sqrt(filter_true2(1:6,size(filter_true,2))))./K, abs(sqrt(filter_approximate2(1:6,size(filter_true,2))))./K];
b=bar(X,Relative_K,0.9);
hold on
plot([0.5,6.5],[1,1],'k--','LineWidth',0.5)
hold on
errorbar([X'-0.15,X'+0.15],Relative_K,SD_error,'k.','LineWidth',1,'MarkerSize',5)
axis([0.5 6.5 0 1.9]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'$k''_1$','$k''_2$','$k''_3$','$k''_4$','$k''_5$','$k''_6$'},'FontSize',font_size);
set(b(1),'FaceColor','r');
set(b(2),'FaceColor','b');
legend(b,{'Particle filters (original models)','Particle filters (reduced models)'},'Location','northwest','FontSize',lengend_font_size);
ylabel({'Relative value'; '(estimate/truth)'},'FontSize',font_size)
xlabel('Reaction constants','FontSize',font_size)
title('(B)','FontSize',title_size)


subplot(4,2,3);
true_value=stairs(TXF,XF(2,:),'k','LineWidth',2);
hold on;
bayesian_estimation=plot(T_filter, filter_true(8,:),'r.--','LineWidth',1,'MarkerSize',marker_size);
hold on
approximate_filter_plot=plot(T_filter, filter_approximate(8,:),'b.--','LineWidth',1,'MarkerSize',marker_size);
hold on
axis([0 T_filter(size(T_filter,2)) 0 1+0.1])
yticks([0:0.2:1])
xticks([0:10:90])
xlabel('Time (minutes)','FontSize',font_size)
ylabel('Activated DNA','FontSize',font_size)
lgd=legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Particle filters (reduced models)'},'Location','south','FontSize',lengend_font_size);
title('(C)','FontSize',title_size)

subplot(4,2,4);
bayesian_estimation=plot(T_filter, sqrt(filter_true2(8,:)),'r--.','LineWidth',1,'MarkerSize',marker_size);
hold on
approximate_filter_plot=plot(T_filter, sqrt(filter_approximate2(8,:)),'b--.','LineWidth',1,'MarkerSize',marker_size);
hold on
axis([0 T_filter(size(T_filter,2)) 0 1+0.1])
yticks([0:0.2:1])
xticks([0:10:90])
%box off; 
xlabel('Time (minutes)','FontSize',font_size)
ylabel({'SD of activated DNA'},'FontSize',font_size)
%title('(B)','position',[25 -0.25],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Particle filters (reduced models)'},'Location','southwest','FontSize',12);
title('(D)','FontSize',title_size)


subplot(4,2,5);
true_value=stairs(TXF,XF(3,:),'black','LineWidth',2);
hold on;
bayesian_estimation=plot(T_filter, filter_true(9,:),'r--.','LineWidth',1,'MarkerSize',marker_size);
hold on
approximate_filter_plot=plot(T_filter, filter_approximate(9,:),'b--.','LineWidth',1,'MarkerSize',marker_size);
hold on
axis([0 T_filter(size(T_filter,2)) 0 fix(max(XF(4,:)))+3])
yticks([0:1:fix(max(XF(3,:))+3)])
xticks([0:10:90])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',font_size)
ylabel('mRNA','FontSize',font_size)
%title('(C)','position',[25 -6.5/4],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Paticle filters (reduced models)'},'Location','northeast','FontSize',12);
title('(E)','FontSize',title_size)

subplot(4,2,6);
bayesian_estimation=plot(T_filter, sqrt(filter_true2(9,:)),'r--.','LineWidth',1,'MarkerSize',marker_size);
hold on
approximate_filter_plot=plot(T_filter, sqrt(filter_approximate2(9,:)),'b--.','LineWidth',1,'MarkerSize',marker_size);
hold on
axis([0 T_filter(size(T_filter,2)) 0 fix(max(XF(4,:)))+3])
yticks([0:1:fix(max(XF(3,:))+3)])
xticks([0:10:90])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',font_size)
ylabel({'SD of mRNA'},'FontSize',font_size)
%title('(B)','position',[25 -0.25],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filters (original models)','Particle filters (reduced models)'},'Location','southwest','FontSize',12);
title('(F)','FontSize',title_size)

subplot(4,2,7);
true_value=stairs(TXF,N*XF(4,:),'black','LineWidth',2);
hold on
bayesian_estimation=plot(T_filter, N*filter_true(10,:),'r--.','LineWidth',1,'MarkerSize',marker_size);
hold on
approximate_filter_plot=plot(T_filter, N*filter_approximate(10,:),'b--.','LineWidth',1,'MarkerSize',marker_size);
hold on
axis([0 T_filter(size(T_filter,2)) 0 (fix(max(XF(4,:)))+1.5)*N])
yticks([0:100:(fix(max(XF(4,:)))+1.5)*N])
xticks([0:10:90])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',font_size)
ylabel('Fluorescent proteins','FontSize',font_size)
ax=gca;
ax.YAxis.Exponent = 2;
%title('(D)','position',[25 -500/4],'FontSize',16);
%legend([true_value,bayesian_estimation,approximate_filter_plot],{'True','Particle filter (original model)','Paticle filter (reduced model)'},'Location','northeast','FontSize',12);
title('(G)','FontSize',title_size)

subplot(4,2,8);
bayesian_estimation=plot(T_filter, N*sqrt(filter_true2(10,:)),'r--.','LineWidth',1,'MarkerSize',marker_size);
hold on
approximate_filter_plot=plot(T_filter, N*sqrt(filter_approximate2(10,:)),'b--.','LineWidth',1,'MarkerSize',marker_size);
hold on
axis([0 T_filter(size(T_filter,2)) 0 (fix(max(XF(4,:)))+1.5)*N])
yticks([0:100:(fix(max(XF(4,:)))+1.5)*N])
xticks([0:10:90])
set(gca,'FontSize',12);
%box off; 
xlabel('Time (minutes)','FontSize',font_size)
ylabel({'SD of FP'},'FontSize',font_size)
%title('(B)','position',[25 -0.25],'FontSize',16);
ax=gca;
ax.YAxis.Exponent = 2;
title('(H)','FontSize',title_size)

