%% Getting PGD stats across sessions 
pathname = uigetdir(pwd);
savepath = pathname;
mkdir([pathname '\Figures'])
pathname = fullfile(pathname);
directory = dir(pathname);
firstflag1 = 1;
firstflag2 = 1;
for idx  = 3:length(directory)-1
    file = directory(idx).name;
    path = directory(idx).folder;
    pathtomat1 = [path,'\',file,'\PGD.mat'];
    pathtomat2 = [path,'\',file,'\PA.mat'];
    if isfile(pathtomat1)
        load(pathtomat1);
        if firstflag1 == 1
            PGDComb = PGD;
            firstflag1 = 0;
        else
            PGDComb =  horzcat(PGDComb,PGD);
        end
    end
    if isfile(pathtomat2)
        load(pathtomat2);
        if firstflag2 == 1
            PAComb = PA;
            firstflag2 = 0;
        else
            PAComb =  horzcat(PAComb,PA);
        end
    end
end
clear PGD PA firstflag1 firstflag2 

%% Ploting avg PGD for different condition vs Cooling 
% For hits 
PGD.PGDHitsBaseline = PGDComb(1).Baseline.cueHits;
for i=2:size(PGDComb,2)
    PGD.PGDHitsBaseline = [PGD.PGDHitsBaseline;PGDComb(i).Baseline.cueHits];
end

PGD.PGDHitsCooling = PGDComb(1).Cooling.cueHits;
for i=2:size(PGDComb,2)
    PGD.PGDHitsCooling = [PGD.PGDHitsCooling;PGDComb(i).Cooling.cueHits];
end

PGD.PGDHitsRecovery = PGDComb(1).Recovery.cueHits;
for i=2:size(PGDComb,2)
    PGD.PGDHitsRecovery = [PGD.PGDHitsRecovery;PGDComb(i).Recovery.cueHits];
end
time = [-1.5:1/parameters.Fs:1.5]';
h=figure();hold on;
plot(time,smooth(mean(PGD.PGDHitsBaseline,1),100,'sgolay',4),'Color', [166/255 14/255 90/255],'LineWidth',1.5);
plot(time,smooth(mean(PGD.PGDHitsCooling,1),100,'sgolay',4),'Color', [53/255 189/255 206/255],'LineWidth',1.5);
plot(time,smooth(mean(PGD.PGDHitsRecovery,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(time,smooth(mean(PGD.PGDHitsBaseline,1)-(std(PGD.PGDHitsBaseline,0,1)/sqrt(size(PGD.PGDHitsBaseline,1))),100,'sgolay',4),'Color', [166/255 14/255 90/255 0.4],'LineWidth',0.5);
plot(time,smooth(mean(PGD.PGDHitsBaseline,1)+(std(PGD.PGDHitsBaseline,0,1)/sqrt(size(PGD.PGDHitsBaseline,1))),100,'sgolay',4),'Color', [166/255 14/255 90/255 0.4],'LineWidth',0.5);
plot(time,smooth(mean(PGD.PGDHitsCooling,1)-(std(PGD.PGDHitsCooling,0,1)/sqrt(size(PGD.PGDHitsCooling,1))),100,'sgolay',4),'Color', [53/255 189/255 206/255 0.4],'LineWidth',0.5);
plot(time,smooth(mean(PGD.PGDHitsCooling,1)+(std(PGD.PGDHitsCooling,0,1)/sqrt(size(PGD.PGDHitsCooling,1))),100,'sgolay',4),'Color', [53/255 189/255 206/255 0.4],'LineWidth',0.5);
plot(time,smooth(mean(PGD.PGDHitsRecovery,1)-(std(PGD.PGDHitsRecovery,0,1)/sqrt(size(PGD.PGDHitsRecovery,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(time,smooth(mean(PGD.PGDHitsRecovery,1)+(std(PGD.PGDHitsRecovery,0,1)/sqrt(size(PGD.PGDHitsRecovery,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
% xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD) for Hits');box off;  legend('Baseline','Cooling','Recovery');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');
saveas(h,[savepath '\Figures\PGDHitsvsCooling.png']);
saveas(h,[savepath '\Figures\PGDHitsvsCooling.fig']);


h=figure();hold on;
plot(time,smooth(zscore(mean(PGD.PGDHitsBaseline,1)),100,'sgolay',4),'Color', [166/255 14/255 90/255],'LineWidth',1.5);
plot(time,smooth(zscore(mean(PGD.PGDHitsCooling,1)),100,'sgolay',4),'Color', [53/255 189/255 206/255],'LineWidth',1.5);
plot(time,smooth(zscore(mean(PGD.PGDHitsRecovery,1)),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
% xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD) for Hits');box off;  legend('Baseline','Cooling','Recovery');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');
saveas(h,[savepath '\Figures\PGDHitsvsCoolingz.png']);
saveas(h,[savepath '\Figures\PGDHitsvsCoolingz.fig']);

%% Phase allignment 

PA.PAHitsBaseline = reshape(PAComb(1).PABaseline.Hit,parameters.rows*parameters.cols,[]);
for i=2:size(PAComb,2)
    PA.PAHitsBaseline = [PA.PAHitsBaseline;reshape(PAComb(i).PABaseline.Hit,parameters.rows*parameters.cols,[])];
end

PA.PAHitsCooling = reshape(PAComb(1).PACooling.Hit,parameters.rows*parameters.cols,[]);
for i=2:size(PAComb,2)
    PA.PAHitsCooling = [PA.PAHitsCooling;reshape(PAComb(i).PACooling.Hit,parameters.rows*parameters.cols,[])];
end

PA.PAHitsRecovery = reshape(PAComb(1).PARecovery.Hit,parameters.rows*parameters.cols,[]);
for i=2:size(PAComb,2)
    PA.PAHitsRecovery = [PA.PAHitsRecovery;reshape(PAComb(i).PARecovery.Hit,parameters.rows*parameters.cols,[])];
end

h=figure();
plot(time,smooth(squeeze(mean(PA.PAHitsBaseline,1,'omitnan')),50,'sgolay',20),'Color',[166/255 14/255 90/255],'LineWidth',1.5); hold on;
plot(time,smooth(squeeze(mean(PA.PAHitsCooling,1,'omitnan')),50,'sgolay',20),'Color',[53/255 189/255 206/255],'LineWidth',1.5); hold on;
plot(time,smooth(squeeze(mean(PA.PAHitsRecovery,1,'omitnan')),50,'sgolay',20),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); hold on;
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xlim([-0.5 1.5]);box off;legend('Baseline','Cool','Recovery');legend('boxoff');set(gca,'TickDir','out','fontsize',14');
title("Phase Alignment : M2 -> M1 Cool");set(gca,'TickDir','out','fontsize',14');
saveas(h,[savepath '\Figures\PAHitsvsCooling.png']);
saveas(h,[savepath '\Figures\PAHitsvsCooling.fig']);