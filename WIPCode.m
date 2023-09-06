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

figure();
histogram(IntanBehaviour.reactionTime);
