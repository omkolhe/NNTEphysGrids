function [wavesStat] = getWaveStats(Waves,parameters,plot)

rows = parameters.rows;
cols = parameters.cols;

% Number of waves
WavesPerTrial = mean( horzcat(Waves.wavesHit(1:end).nWaves),'all');
totalWaves =  sum( horzcat(Waves.wavesHit(1:end).nWaves));

WavesPerTrialMiss = mean( horzcat(Waves.wavesMiss(1:end).nWaves),'all');
totalWavesMiss =  sum( horzcat(Waves.wavesMiss(1:end).nWaves));

% Wave Speed stats
speedComb = horzcat(Waves.wavesHit(1:end).speed);
avgSpeed = mean(speedComb);

speedCombMiss = horzcat(Waves.wavesMiss(1:end).speed);
avgSpeedMiss = mean(speedCombMiss);

% Perform the t-test.
[p, t] = ranksum(speedComb, speedCombMiss);
% Print the results.
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wave speeds in motion Initiation and termination');
    subplot(2,1,1);
    histfit(speedComb,100,'kernel');
    xline(avgSpeed,'-r',{'Mean speed = ' num2str(avgSpeed) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Initiation');
    subplot(2,1,2);
    histfit(speedCombMiss,100,'kernel');
    xline(avgSpeedMiss,'-r',{'Mean speed = ' num2str(avgSpeedMiss) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Termination');

    figure('Name','Wave speeds in motion Initiation and termination');
    group = [ones(size(speedComb')); 2.*ones(size(speedCombMiss'))];
    boxplot([speedComb';speedCombMiss'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Termination'});
    ylabel('Wave speed in cm/s');
end


% Wavelength stats
lComb = horzcat(Waves.wavesHit(1:end).wavelength);
avgl = mean(lComb);

lCombMiss = horzcat(Waves.wavesMiss(1:end).wavelength);
avglMiss = mean(lCombMiss);

% Perform the t-test.
[p, t] = ranksum(lComb, lCombMiss);
% Print the results.
disp('Wavelength')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wavelength in motion Initiation and termination');
    subplot(2,1,1);
    histfit(lComb,100,'kernel');
    xline(avgl,'-r',{'Mean wavelength = ' num2str(avgl) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength: Motion Initiation');
    subplot(2,1,2);
    histfit(lCombMiss,100,'kernel');
    xline(avglMiss,'-r',{'Mean wavelength = ' num2str(avglMiss) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength : Motion Termination');

    figure('Name','Wavelength in motion Initiation and termination');
    group = [ones(size(lComb')); 2.*ones(size(lCombMiss'))];
    boxplot([lComb';lCombMiss'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Termination'});
    xlabel('Wavelength in cm');
end


% Wave direction stats
dirComb = horzcat(Waves.wavesHit(1:end).waveDir);
avgDir = mean(dirComb);
dirCombMiss = horzcat(Waves.wavesMiss(1:end).waveDir);
avgDirMiss = mean(dirCombMiss);

[p, t] = ranksum(dirComb, dirCombMiss);
% Print the results.
disp('Wave Direction')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Polar Histogram for wave direction in Motion Initiation and Termination');
    subplot(2,1,1);
    polarhistogram(dirComb,30);
    title('Wave Direction : Motion Initiation');
    subplot(2,1,2);
    polarhistogram(dirCombMiss,30);
    title('Wave Direction : Motion Termination');

    figure('Name','Wave Direction in motion Initiation and termination');
    group = [ones(size(dirComb'));2.*ones(size(dirCombMiss'))];
    boxplot([rad2deg(dirComb)';rad2deg(dirCombMiss)'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Termination'});
end

% Wave source points stats
sourceComb = horzcat(Waves.wavesHit(1:end).source);
sourceDen = zeros(rows,cols);
sourceCombMiss = horzcat(Waves.wavesMiss(1:end).source);
sourceDenMiss = zeros(rows,cols);

for j=1:size(sourceComb,2)
    sourceDen(sourceComb(2,j),sourceComb(1,j)) = sourceDen(sourceComb(2,j),sourceComb(1,j)) + 1;
end
maxSourcePoint = max(sourceComb);

for j=1:size(sourceCombMiss,2)
    sourceDenMiss(sourceCombMiss(2,j),sourceCombMiss(1,j)) = sourceDenMiss(sourceCombMiss(2,j),sourceCombMiss(1,j)) + 1;
end
maxSourcePointMiss = max(sourceCombMiss);

if plot == 1
    figure('Name','Spatial map of source points in Motion Initiation and Termination'); 
    subplot(2,1,1);
    imagesc(sourceDen);set(gca,'YDir','normal');
    title('Spatial map of sources points : Motion Inititaion'); colorbar;
    subplot(2,1,2);
    imagesc(sourceDenMiss);set(gca,'YDir','normal');
    title('Spatial map of sources points : Motion Termination'); colorbar;
end 

wavesStat.evaluationPoints =  horzcat(Waves.wavesHit(1:end).evaluationPoints);
wavesStat.evaluationPointsMiss =  horzcat(Waves.wavesMiss(1:end).evaluationPoints);
wavesStat.speed = speedComb;
wavesStat.avgSpeed = avgSpeed;
wavesStat.speedMiss = speedCombMiss;
wavesStat.avgSpeedMiss = avgSpeedMiss;

wavesStat.velDir = dirComb;
wavesStat.avgDir = avgDir;
wavesStat.velDirMiss = dirCombMiss;
wavesStat.avgDirMiss = avgDirMiss;

wavesStat.sourcePoints = sourceComb;
wavesStat.maxSourcePoint = maxSourcePoint;
wavesStat.sourcePointsMiss = sourceCombMiss;
wavesStat.maxSourcePointMiss = maxSourcePointMiss;

wavesStat.WavesPerTrial = WavesPerTrial;
wavesStat.totalWaves = totalWaves;
wavesStat.WavesPerTrialMiss = WavesPerTrialMiss;
wavesStat.totalWavesMiss = totalWavesMiss;

end

