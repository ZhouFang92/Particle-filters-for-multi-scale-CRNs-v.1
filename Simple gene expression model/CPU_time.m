function CPU_time(time_cf,time_cr, time_df, time_dr)

%X = categorical({'Continuous-time observations','Discrete-time observations'});
%X = reordercats(X,{'Continuous-time observations','Discrete-time observations'});
X=[1,2];
b=bar(X,[time_cf,time_cr; time_df, time_dr]/60,0.9);
ylabel('Time in minutes','FontSize',16);
%set(gca,'FontSize',12);
set(gca,'xticklabel',{'Continuous-time observations','Discrete-time observations'},'FontSize',12);

set(b(1),'FaceColor','r');
set(b(2),'FaceColor','b');

legend(b,{'Particle filters (original models)','Particle filters (reduced models)'},'Location','northeast','FontSize',12);
title('CPU Time','FontSize',16);

text(X(1)-0.18,(time_cf+20)/60,num2str(time_cf/60,3));
text(X(1)+0.12,(time_cr+20)/60,num2str(time_cr/60,2));
text(X(2)-0.18,(time_df+20)/60,num2str(time_df/60,3));
text(X(2)+0.12,(time_dr+20)/60,num2str(time_dr/60,2));
