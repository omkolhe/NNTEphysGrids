%% Work In Progress code - Scratch work for new fucntions

b = 16;
% Getting Amplitude of waves
wavePhases = angle(Waves.wavesHit(2).p(:,:,Waves.wavesHit(2).waveTime(b,1):Waves.wavesHit(2).waveTime(b,2)));
phi = mapAngle360(rad2deg(reshape(wavePhases,32,size(wavePhases,3))));

waveAmp = abs(Waves.wavesHit(2).p(:,:,Waves.wavesHit(2).waveTime(b,1):Waves.wavesHit(2).waveTime(b,2)));
amp = reshape(waveAmp,32,size(waveAmp,3));

figure;
subplot(2,1,1);
plot(phi');
subplot(2,1,2);
plot(amp');

xijtHit = cellfun(@(s) s(1,1,1),xgpHit);
xijtMiss = cellfun(@(s) s(1,1,1),xgpMiss);

infomap = reshape(MIPhase,[],size(MIPhase,3));
infopeak = mean(infomap,2);
[~,indexsort] = sort(infopeak,1,'descend');

sortedInfomap = infomap(indexsort,:);

figure();
title("Mututal Information across all electrodes - Phase")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,sortedInfomap); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Information (bits)';
xline(0,'-w','Cue','LabelVerticalAlignment','top');

figure();
for i=1:size(IntanBehaviour.missTrace,2)
    plot(IntanBehaviour.missTrace(i).time,IntanBehaviour.missTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.missTrace(1).time,mean(horzcat(IntanBehaviour.missTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Misses');box off;

%% Plotting histogram of reaction time
figure();
histogram(IntanBehaviour.reactionTime);

%% Sorting PA 
a = reshape(PAHit,[],size(PAHitz,3));
[~,aidx] = sort(sum(a,2),1,'descend');
b = a(aidx,:);

figure();
title("Phase Alignment across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,b); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');


%% Performing wave stats for a time window
st = 1500;
sp = 2000;
dirComb = horzcat(Waves.wavesHit(1:end).waveDir);
evalPointsHit = vertcat(Waves.wavesHit(1:end).evaluationPoints)';
dirComb = dirComb(evalPointsHit >=st & evalPointsHit <= sp);
avgDir = mean(dirComb);
dirCombMiss = horzcat(Waves.wavesMiss(1:end).waveDir);
evalPointsMiss = vertcat(Waves.wavesMiss(1:end).evaluationPoints)';
dirCombMiss = dirCombMiss(evalPointsMiss >=st & evalPointsMiss <= sp);
avgDirMiss = mean(dirCombMiss);

[p, t] = ranksum(mapAngle360(rad2deg(dirComb)), mapAngle360(rad2deg(dirCombMiss)));
% Print the results.
disp('Wave Direction')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);

figure('Name','Polar Histogram for wave direction in Hits and Misses');
subplot(1,4,1);
polarhistogram(dirComb,30,'FaceColor','blue');
title('Wave Direction :   Hits');box off;
subplot(1,4,4);
polarhistogram(dirCombMiss,30);
title('Wave Direction :   Misses');box off;

figure('Name','Wave Direction in   Hits and Misses');
group = [ones(size(dirComb'));2.*ones(size(dirCombMiss'))];
boxplot([mapAngle360(rad2deg(dirComb))';mapAngle360(rad2deg(dirCombMiss))'],group,'BoxStyle','filled','PlotStyle','compact');
set(gca,'XTickLabel',{'Hits','Misses'});box off;


figure('Name','Histogram for wave direction in   Hits and Misses');
subplot(2,1,1);
histogram(mapAngle360(rad2deg(dirComb)),60);
title('Wave Direction :  Hits');box off;
subplot(2,1,2);
histogram(mapAngle360(rad2deg(dirCombMiss)),60);
title('Wave Direction :  Misses');box off;


%% Phase allignement but for phase gradient direction
xgp = arrayfun(@(s) s.pd, Waves.wavesHit, 'UniformOutput', false);
PGAHit = getPhaseAlignment(xgp,parameters);
xgp = arrayfun(@(s) s.pd, Waves.wavesMiss, 'UniformOutput', false);
PGAMiss = getPhaseAlignment(xgp,parameters);

figure();
subplot(2,1,1);
title("Phase Alignment across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PGAHit,[],size(PGAHit,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Phase Alignment across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PGAMiss,[],size(PGAMiss,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');


figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PGAHit,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PGAMiss,[1 2])),'-k','LineWidth',1);
ylabel("Phase Gradient Alignment"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Phase Gradient Alignment for Hits');box off;% legend('Hits','Miss');

%% Mutual Information
pdHit = arrayfun(@(s) s.insts(:,:,:), Waves.wavesHit, 'UniformOutput', false);
pdMiss = arrayfun(@(s) s.insts(:,:,:), Waves.wavesMiss, 'UniformOutput', false);
[MIpd] = getMutualInformation(pdHit,pdMiss,parameters);

figure();
title("Mututal Information across all electrodes - Phase")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(MIpd,[],size(MIpd,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Information (bits)';
xline(0,'-w','Cue','LabelVerticalAlignment','top');

%%

for i = 1:4
    st = (i-1)*500;
    sp = (i)*500;
    dirComb1 = horzcat(Waves.wavesHit(1:end).waveDuration);
    evalPointsHit = vertcat(Waves.wavesHit(1:end).evaluationPoints)';
    WaveComb(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesMiss(1:end).waveDuration);
    evalPointsMiss = vertcat(Waves.wavesMiss(1:end).evaluationPoints)';
    WaveComb(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);
end



figure(); hold on;
plot(1:4,arrayfun(@(s) mean(s.dir,'all'), dirComb));
plot(1:4,arrayfun(@(s) mean(s.dir,'all'), dirCombMiss));

a = arrayfun(@(s) ranksum(s.Hit,s.Miss) , WaveComb);

%% 
for z = 1:10
    disp(z);
    if z == 5
        z = 2;
    end
end
