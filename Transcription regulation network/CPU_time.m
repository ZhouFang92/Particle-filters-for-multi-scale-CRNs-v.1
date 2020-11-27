function CPU_time(time_cf,time_cr, time_df, time_dr)

time_cf=log10(time_cf/3600)+3;
time_cr=log10(time_cr/3600)+3;
time_df=log10(time_df/3600)+3;
time_dr=log10(time_dr/3600)+3;
%X = categorical({'Continuous-time observations','Discrete-time observations'});
%X = reordercats(X,{'Continuous-time observations','Discrete-time observations'});
X=[1,2];
b=bar(X,[time_cf,time_cr; time_df, time_dr],0.9);
ylabel('Time in hours','FontSize',16);
%set(gca,'FontSize',12);
set(gca,'xticklabel',{'Continuous-time observations','Discrete-time observations'},'FontSize',12);
axis([0.5 2.5 0 5.2])
yticks([0:1:5])
set(gca,'yticklabel',{'0.001','0.01','0.1','1','10','100'},'FontSize',12);


set(b(1),'FaceColor','r');
set(b(2),'FaceColor','b');

legend(b,{'Particle filters (original models)','Particle filters (reduced models)'},'Location','northeast','FontSize',12);
title('CPU Time','FontSize',16);

text(X(1)-0.21,time_cf+0.1,[num2str(roundn(10^(time_cf-3),-2)), ' h']);
text(X(1)+0.1,time_cr+0.1,[num2str(roundn(10^(time_cr-3)*3600,-2)), ' s']);
text(X(2)-0.21,time_df+0.1,[num2str(roundn(10^(time_df-3),-2)), ' h']);
text(X(2)+0.1,time_dr+0.1,[num2str(roundn(10^(time_dr-3)*3600,-2)), ' s']);
