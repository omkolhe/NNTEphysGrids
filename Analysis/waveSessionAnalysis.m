%% Plotting rasters 
%  Waves during hits 
plotWaveRaster(WavesBaseline.wavesHit,WavesOpto.wavesHit,IntanBehaviourBaseline.cueHitTrace,IntanBehaviourOpto.cueHitTrace,parameters);
sgtitle('Wave Rasters for Hits: Baseline vs Opto');
subplot(4,1,1);
RTTraceTime = (IntanBehaviourBaseline.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r');