% Grids
[SpectrogramGrid.hitAvgSpectrogram, SpectrogramGrid.hitSpectrogramCWT,SpectrogramGrid.AvgHitTrace ,SpectrogramGrid.fwt] = getAvgSpectogram(IntanBehaviour.cueHitTrace,parameters,[5 80],0);
[SpectrogramGrid.missAvgSpectrogram, SpectrogramGrid.missSpectrogramCWT,SpectrogramGrid.AvgMissTrace,SpectrogramGrid.fwt] = getAvgSpectogram(IntanBehaviour.cueMissTrace,parameters,[5 80],0);
[SpectrogramGrid.hitRewardAvgSpectrogram, SpectrogramGrid.hitRewardSpectrogramCWT,SpectrogramGrid.AvgHitRewardTrace ,SpectrogramGrid.fwt] = getAvgSpectogram(IntanBehaviour.hitTrace,parameters,[5 80],0);
[SpectrogramGrid.FAAvgSpectrogram, SpectrogramGrid.FASpectrogramCWT,SpectrogramGrid.AvgFATrace,SpectrogramGrid.fwt] = getAvgSpectogram(IntanBehaviour.missTrace,parameters,[5 80],0);
[SpectrogramGrid.hitMIAvgSpectrogram, SpectrogramGrid.hitMISpectrogramCWT,SpectrogramGrid.AvgHitMITrace ,SpectrogramGrid.fwt] = getAvgSpectogram(IntanBehaviour.MIHitTrace,parameters,[5 80],0);
[SpectrogramGrid.FAMIAvgSpectrogram, SpectrogramGrid.FAMISpectrogramCWT,SpectrogramGrid.AvgFAMITrace,SpectrogramGrid.fwt] = getAvgSpectogram(IntanBehaviour.MIFATrace,parameters,[5 80],0);

% Shanks
[SpectrogramShank.hitAvgSpectrogram, SpectrogramShank.hitSpectrogramCWT,SpectrogramShank.AvgHitTrace ,SpectrogramShank.fwt] = getAvgSpectogram(IntanBehaviour.cueHitTrace,parameters,[5 80],1);
[SpectrogramShank.missAvgSpectrogram, SpectrogramShank.missSpectrogramCWT,SpectrogramShank.AvgMissTrace,SpectrogramShank.fwt] = getAvgSpectogram(IntanBehaviour.cueMissTrace,parameters,[5 80],1);
[SpectrogramShank.hitRewardAvgSpectrogram, SpectrogramShank.hitRewardSpectrogramCWT,SpectrogramShank.AvgHitRewardTrace ,SpectrogramShank.fwt] = getAvgSpectogram(IntanBehaviour.hitTrace,parameters,[5 80],1);
[SpectrogramShank.FAAvgSpectrogram, SpectrogramShank.FASpectrogramCWT,SpectrogramShank.AvgFATrace,SpectrogramShank.fwt] = getAvgSpectogram(IntanBehaviour.missTrace,parameters,[5 80],1);
[SpectrogramShank.hitMIAvgSpectrogram, SpectrogramShank.hitMISpectrogramCWT,SpectrogramShank.AvgHitMITrace ,SpectrogramShank.fwt] = getAvgSpectogram(IntanBehaviour.MIHitTrace,parameters,[5 80],1);
[SpectrogramShank.FAMIAvgSpectrogram, SpectrogramShank.FAMISpectrogramCWT,SpectrogramShank.AvgFAMITrace,SpectrogramShank.fwt] = getAvgSpectogram(IntanBehaviour.MIFATrace,parameters,[5 80],1);

% Spectrogram = SpectrogramGrid;
% Spectrogram = SpectrogramShank;

% Global average spectogram
figure('Name','Trial Averaged Wavelet Spectrogram for Hits & Misses');
subplot(1,2,1);
plotSpectrogram(10*log10((squeeze(Spectrogram.hitAvgSpectrogram))),IntanBehaviour.cueHitTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Hits','Time (s)','Frequency (Hz)')
caxis([-2 13]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueHitTrace(1).time,Spectrogram.AvgHitTrace,'-w','LineWidth',2.5);
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]);
subplot(1,2,2);
plotSpectrogram(10*log10((squeeze(Spectrogram.missAvgSpectrogram))),IntanBehaviour.cueMissTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Misses','Time (s)','Frequency (Hz)')
caxis([-2 13]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueMissTrace(1).time,Spectrogram.AvgMissTrace,'-w','LineWidth',2.5);
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]); box off;

% Global average spectogram
figure('Name','Trial Averaged Wavelet Spectrogram for Hits & FA');
subplot(1,2,1);
plotSpectrogram(10*log10((squeeze(Spectrogram.hitRewardAvgSpectrogram))),IntanBehaviour.hitTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Hits','Time (s)','Frequency (Hz)')
caxis([-2 15]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.hitTrace(1).time,Spectrogram.AvgHitRewardTrace,'-w','LineWidth',2.5);
xline(0,'--r','Reward','LabelVerticalAlignment','top');
% xline(-1*mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]);
subplot(1,2,2);
plotSpectrogram(10*log10((squeeze(Spectrogram.FAAvgSpectrogram))),IntanBehaviour.missTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Misses','Time (s)','Frequency (Hz)')
caxis([-1 15]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.missTrace(1).time,Spectrogram.AvgFATrace,'-w','LineWidth',2.5);
% xline(0.5,'--r','Threshold Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]); box off;


% Global average spectogram
figure('Name','Trial Averaged Wavelet Spectrogram for Hits & FA with MI');
subplot(1,2,1);
plotSpectrogram(10*log10((squeeze(Spectrogram.hitMIAvgSpectrogram))),IntanBehaviour.MIHitTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Hits','Time (s)','Frequency (Hz)')
caxis([-2 15]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.MIHitTrace(1).time,Spectrogram.AvgHitMITrace,'-w','LineWidth',2.5);
xline(0,'--r','MI','LabelVerticalAlignment','top');
% xline(-1*mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]);
subplot(1,2,2);
plotSpectrogram(10*log10((squeeze(Spectrogram.FAMIAvgSpectrogram))),IntanBehaviour.MIFATrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for FA','Time (s)','Frequency (Hz)')
caxis([-2 15]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.MIFATrace(1).time,Spectrogram.AvgFAMITrace,'-w','LineWidth',2.5);
xline(0,'--r','MI','LabelVerticalAlignment','top');
% xline(0.5,'--r','Threshold Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]); box off;