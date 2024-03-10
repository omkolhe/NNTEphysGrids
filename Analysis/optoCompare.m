%% Comaparing Behaviour 
[p,h] = ranksum(IntanBehaviourBaseline.reactionTime,IntanBehaviourOpto.reactionTime)
plotBox2(IntanBehaviourBaseline.reactionTime,IntanBehaviourOpto.reactionTime);
% xL=xlim;
% yL=ylim;
% text(0.995*xL(2),0.995*yL(2),['p-val = ' num2str(p)],'HorizontalAlignment','right','VerticalAlignment','top')
ylabel('Reaction Time (s)'); title('M2 -> VM Opto');subtitle(['p-val = ' num2str(p)]);
xtix = {'Baseline','Opto'}; xtixloc = [1 2]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
set(gca,'TickDir','out','fontsize',14');

%% Comparing PA 
PABaseline = getPA(IntanBehaviourBaseline,0,1,0,parameters,0);
PAOpto = getPA(IntanBehaviourOpto,0,1,0,parameters,0);

figure();
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(squeeze(mean(PABaseline.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[0.7 0.7 0.7],'LineWidth',1.5); hold on;
plot(IntanBehaviourOpto.cueHitTrace(1).time,smooth(squeeze(mean(PAOpto.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xlim([-0.5 1.5]);box off;legend('Baseline','Opto');legend('boxoff');set(gca,'TickDir','out','fontsize',14');
title("Phase Alignment : M2 -> VM Opto");

%% Comparing PDG 
PGD.avgPGDBaseline = mean(vertcat(WavesBaseline.wavesHit.PGD),1);
PGD.avgPGDOpto = mean(vertcat(WavesOpto.wavesHit.PGD),1);

figure(); hold on;
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(PGD.avgPGDBaseline,50,'sgolay',20),'Color',[0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviourOpto.cueHitTrace(1).time,smooth(PGD.avgPGDOpto,50,'sgolay',20),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off; legend('Baseline','Opto');legend('boxoff');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');

%% Comparing Wave Properties - Speed 


%% Comparing Wave Properties - Direction
% Wave direction stats
dirCombBaseline = horzcat(WavesBaseline.wavesHit(1:end).waveDir);
avgDirBaseline = mean(dirCombBaseline);
dirCombOpto = horzcat(WavesOpto.wavesHit(1:end).waveDir);
avgDirOpto = mean(dirCombOpto);

[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);

% Fit VonMisesMixture Model 
nComponents = 2;    
fittedVmBaseline = fitmvmdist(dirCombBaseline',nComponents,'MaxIter', 250);
fittedVmOpto = fitmvmdist(dirCombOpto',nComponents,'MaxIter', 250);
angles = linspace(-pi, pi, 1000)';
likelihoodsBaseline = fittedVmBaseline.pdf(angles);
likelihoodsOpto = fittedVmOpto.pdf(angles);

figure('Name','Polar Histogram for wave direction in Baseline and Opto');
subplot(1,2,1);
h=polarhistogram(dirCombBaseline,18,'Normalization','pdf','FaceColor',[153/255 153/255 255/255],'FaceAlpha',1,'EdgeAlpha',0.3);hold on;
polarplot(angles,likelihoodsBaseline,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);rlim([0 0.4]);
polarplot([fittedVmBaseline.mu';fittedVmBaseline.mu'],[0 0;fittedVmBaseline.pdf(fittedVmBaseline.mu)'],'-r','LineWidth',1.5);
title('Wave Direction : Baseline');box off;
subplot(1,2,2);
polarhistogram(dirCombOpto,18,'Normalization','pdf','FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.3,'EdgeAlpha',0.3);hold on;
polarplot(angles,likelihoodsOpto,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);rlim([0 0.4]);
polarplot([fittedVmOpto.mu';fittedVmOpto.mu'],[0 0;fittedVmOpto.pdf(fittedVmOpto.mu)'],'-r','LineWidth',1.5);
title('Wave Direction : Opto');box off;



