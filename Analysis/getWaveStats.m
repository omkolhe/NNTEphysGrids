function [wavesStat] = getWaveStats(Waves,parameters,plot)

rows = parameters.rows;
cols = parameters.cols;

% Number of waves
WavesPerTrial = mean( horzcat(Waves(1:end).nWaves),'all');
totalWaves =  sum( horzcat(Waves(1:end).nWaves));

% Wave Speed stats
speedComb = horzcat(Waves(1:end).speed);
avgSpeed = mean(speedComb);
if plot == 1
    figure();histogram(speedComb,100);
    xline(avgSpeed,'-r',{'Mean speed = ' num2str(avgSpeed) ' cm/s'});
    xlabel('Wave speed in cm/s');
    ylabel('Frequency');
    title('Histogram of wave speeds');
end

% Wave direction stats
dirComb = horzcat(Waves(1:end).waveDir);
avgDir = mean(dirComb);
if plot == 1
    figure();polarhistogram(dirComb,30);
    title('Polar Histogram of wave direction');
end

% Wave source points stats
sourceComb = horzcat(Waves(1:end).source);
sourceDen = zeros(rows,cols);
for j=1:size(sourceComb,2)
    sourceDen(sourceComb(2,j),sourceComb(1,j)) = sourceDen(sourceComb(2,j),sourceComb(1,j)) + 1;
end
maxSourcePoint = max(sourceComb);
if plot == 1
    figure(); imagesc(sourceDen);set(gca,'YDir','normal');
    title('Spatial map of sources points accross all trials'); colorbar;
end 

wavesStat.evaluationPoints =  horzcat(Waves(1:end).evaluationPoints);
wavesStat.speed = speedComb;
wavesStat.avgSpeed = avgSpeed;
wavesStat.velDir = dirComb;
wavesStat.avgDir = avgDir;
wavesStat.sourcePoints = sourceComb;
wavesStat.maxSourcePoint = maxSourcePoint;
wavesStat.WavesPerTrial = WavesPerTrial;
wavesStat.totalWaves =totalWaves;

end

