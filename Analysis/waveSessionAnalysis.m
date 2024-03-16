%% Plotting rasters 
%  Waves during hits vs miss
plotWaveRaster(Waves.wavesHit,Waves.wavesMiss,IntanBehaviour.cueHitTrace,IntanBehaviour.cueMissTrace,parameters);
sgtitle('Wave Rasters for Hits vs Misses');
subplot(4,1,1);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');

%  Waves during Hits vs FA
plotWaveRaster(Waves.wavesHitReward,Waves.wavesFA,IntanBehaviour.hitTrace,IntanBehaviour.missTrace,parameters);
sgtitle('Wave Rasters for Hits vs FA');

%  Waves during Hits vs FA
plotWaveRaster(Waves.wavesMIHit,Waves.wavesMIFA,IntanBehaviour.MIHitTrace,IntanBehaviour.MIFATrace,parameters);
sgtitle('Wave Rasters for Hits vs FA');


%% Wave Direction
% Waves Hits vs Miss
W = vertcat(selectWaves(Waves.wavesHit,1,3000).waveStart);figure();
prop = arrayfun(@(s) s.waveDir, selectWaves(Waves.wavesHit,1,3000), 'UniformOutput', false);
ax1 = subplot(2,1,1);
rasterPlotPropColor(W,prop,[],1);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');xlim([1 size(W,2)]);
title('Wave Hits')

W = vertcat(selectWaves(Waves.wavesMiss,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(Waves.wavesMiss,1,3000), 'UniformOutput', false);
ax2 = subplot(2,1,2);
rasterPlotPropColor(W,prop,[],1);
title('Wave Miss')
linkaxes([ax1,ax2],'x');

% Plotting Wave Rasters colorcoded by Direction
srt = 1500;
stp = 1800;
figure();
subplot(1,2,1);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHit(1:end),srt,stp),[]);
title('Hits');
subplot(1,2,2);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMiss(1:end),srt,stp),[]);
title('Miss');
sgtitle('Wave Direction')

[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);
