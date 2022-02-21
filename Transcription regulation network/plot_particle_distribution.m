function plot_particle_distribution(K,parameters_CF,parameters_CR,parameters_DF,parameters_DR)

% This function draws the scatter plot of the particles at final time. K3
% and K4

f = figure;
f.Units='pixels';
f.OuterPosition=[10 10 800 400];

font_size=12;
title_size=15;
marker_size=15;
lengend_font_size=9;
M=size(parameters_CF.particles(3,:),2);


% continuous-time observations: particle filter original model

subplot(2,4,1);
x=unique(parameters_CF.particles(3:4,:)','rows'); % find different particles
probability=zeros(size(x,1),1);
for j=1:size(x,1)
   probability(j)= size(find(parameters_CF.particles(3:4,:)'==x(j,:)),1)/M;
end
scatter1=scatter(x(:,1),x(:,2),'filled','red');
scatter1.AlphaData = probability;
scatter1.MarkerFaceAlpha = 'flat';
hold on
plot(K(3),K(4),'k.','MarkerSize',marker_size)
hold on
box on
axis([0.7 0.9 0.3 0.5])
xlabel('$k''_3$','FontSize',font_size,'Interpreter','latex')
ylabel('$k''_4$','FontSize',font_size,'Interpreter','latex')
title('(A)','FontSize',title_size)

subplot(2,4,5);
stairs([0:0.5:90],log10(parameters_CF.diversity(1:2:end)),'r-','LineWidth',1);
axis([0 90 log10(10) log10(10000)])
xticks([0:30:90])
xticklabels([0:30:90]*100/60)
yticks(log10([10 100 1000 10000 100000]))
yticklabels({'10' '10^2', '10^3','10^4','10^5'})
xlabel('Time (minutes)','FontSize',font_size)
title('(E)','FontSize',title_size)
%ylabel({'Diversity of particles'; 'for estimating parameters'},'FontSize',font_size)

% continuous-time observations: particle filter reduced model
subplot(2,4,2);
x=unique(parameters_CR.particles(3:4,:)','rows'); % find different particles
probability=zeros(size(x,1),1);
for j=1:size(x,1)
   probability(j)= size(find(parameters_CR.particles(3:4,:)'==x(j,:)),1)/M;
end
scatter1=scatter(x(:,1),x(:,2),'filled','blue');
scatter1.AlphaData = probability;
scatter1.MarkerFaceAlpha = 'flat';
hold on
plot(K(3),K(4),'k.','MarkerSize',marker_size)
hold on
box on
axis([0.7 0.9 0.3 0.5])
xlabel('$k''_3$','FontSize',font_size,'Interpreter','latex')
ylabel('$k''_4$','FontSize',font_size,'Interpreter','latex')
title('(B)','FontSize',title_size)

subplot(2,4,6);
stairs([0:0.5:90],log10(parameters_CR.diversity(1:2:end)),'b-','LineWidth',1);
axis([0 90 log10(10) log10(10000)])
xticks([0:30:90])
xticklabels([0:30:90]*100/60)
yticks(log10([10 100 1000 10000 100000]))
yticklabels({'10' '10^2', '10^3','10^4','10^5'})
xlabel('Time (minutes)','FontSize',font_size)
title('(F)','FontSize',title_size)


% discrete-time observations: particle filter full model
subplot(2,4,3);
x=unique(parameters_DF.particles(3:4,:)','rows'); % find different particles
probability=zeros(size(x,1),1);
for j=1:size(x,1)
   probability(j)= size(find(parameters_DF.particles(3:4,:)'==x(j,:)),1)/M;
end
scatter1=scatter(x(:,1),x(:,2),'filled','red');
scatter1.AlphaData = probability;
scatter1.MarkerFaceAlpha = 'flat';
hold on
plot(K(3),K(4),'k.','MarkerSize',marker_size)
hold on
box on
axis([0.7 0.9 0.3 0.5])
xlabel('$k''_3$','FontSize',font_size,'Interpreter','latex')
ylabel('$k''_4$','FontSize',font_size,'Interpreter','latex')
title('(C)','FontSize',title_size)

subplot(2,4,7);
stairs([0:2:90],log10(parameters_DF.diversity),'r-','LineWidth',1);
axis([0 90 log10(10) log10(10000)])
xticks([0:30:90])
xticklabels([0:30:90]*100/60)
yticks(log10([10 100 1000 10000 100000]))
yticklabels({'10' '10^2', '10^3','10^4','10^5'})
xlabel('Time (minutes)','FontSize',font_size)
title('(G)','FontSize',title_size)


% discrete-time observations: particle filter reduced model
subplot(2,4,4);
x=unique(parameters_DR.particles(3:4,:)','rows'); % find different particles
probability=zeros(size(x,1),1);
for j=1:size(x,1)
   probability(j)= size(find(parameters_DR.particles(3:4,:)'==x(j,:)),1)/M;
end
scatter1=scatter(x(:,1),x(:,2),'filled','blue');
scatter1.AlphaData = probability;
scatter1.MarkerFaceAlpha = 'flat';
hold on
plot(K(3),K(4),'k.','MarkerSize',marker_size)
hold on
box on
axis([0.7 0.9 0.3 0.5])
xlabel('$k''_3$','FontSize',font_size,'Interpreter','latex')
ylabel('$k''_4$','FontSize',font_size,'Interpreter','latex')
title('(D)','FontSize',title_size)

subplot(2,4,8);
stairs([0:2:90],log10(parameters_DR.diversity),'b-','LineWidth',1);
axis([0 90 log10(10) log10(10000)])
xticks([0:30:90])
xticklabels([0:30:90]*100/60)
yticks(log10([10 100 1000 10000 100000]))
yticklabels({'10' '10^2', '10^3','10^4','10^5'})
xlabel('Time (minutes)','FontSize',font_size)
title('(H)','FontSize',title_size)

