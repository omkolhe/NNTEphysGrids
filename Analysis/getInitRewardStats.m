function [wavesStat] = getInitRewardStats(Waves,parameters,plot)

rows = parameters.rows;
cols = parameters.cols;

cueTime = (parameters.Fs*parameters.windowBeforePull) + 1;
initWindow = 0.8*parameters.Fs;
rewardWindow = 0.8*parameters.Fs;

allWavesEvalPoints = horzcat(Waves.wavesHit(1:end).evaluationPoints);
allWavesSpeed = horzcat(Waves.wavesHit(1:end).speed);
allWavesDir = horzcat(Waves.wavesHit(1:end).waveDir);
allWavesWaveLength = horzcat(Waves.wavesHit(1:end).wavelength);
allWavesSource = horzcat(Waves.wavesHit(1:end).source);

initWavesIndx = find(allWavesEvalPoints<cueTime & allWavesEvalPoints>cueTime-initWindow);
initWavesSpeed = allWavesSpeed(initWavesIndx);
initWavesDir = allWavesDir(initWavesIndx);
initWavesWavelength = allWavesWaveLength(initWavesIndx);
initWavesSource = allWavesSource(:,initWavesIndx);


rewardWavesIndx = find(allWavesEvalPoints>cueTime & allWavesEvalPoints<cueTime+rewardWindow);
rewardWavesSpeed = allWavesSpeed(rewardWavesIndx);
rewardWavesDir = allWavesDir(rewardWavesIndx);
rewardWavesWavelength = allWavesWaveLength(rewardWavesIndx);
rewardWavesSource = allWavesSource(:,rewardWavesIndx);

% Number of waves 
nWavesInit = size(initWavesIndx);
nWavesReward = size(rewardWavesIndx);

% Wave speed stats
avgSpeedInit = mean(initWavesSpeed);
avgSpeedReward = mean(rewardWavesSpeed);

[p,t] = ranksum(initWavesSpeed,rewardWavesSpeed);
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);

if plot == 1
    figure('Name','Histogram of wave speeds in motion Initiation and Reward');
    subplot(2,1,1);
    histfit(initWavesSpeed,100,'kernel');
    xline(avgSpeedInit,'-r',{'Mean speed = ' num2str(avgSpeedInit) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Initiation');
    subplot(2,1,2);
    histfit(rewardWavesSpeed,100,'kernel');
    xline(avgSpeedReward,'-r',{'Mean speed = ' num2str(avgSpeedReward) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Reward');

    figure('Name','Wave speeds in motion Initiation and reward');
    group = [ones(size(initWavesSpeed')); 2.*ones(size(rewardWavesSpeed'))];
    boxplot([initWavesSpeed';rewardWavesSpeed'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Reward'});
    ylabel('Wave speed in cm/s');
end


% Wavelength stats
avgWavelengthInit = mean(initWavesWavelength);
avgWavelengthReward = mean(rewardWavesWavelength);

[p,t] = ranksum(initWavesWavelength,rewardWavesWavelength);
disp('Wave Wavelength')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);

if plot == 1
    figure('Name','Histogram of Wavelengths in motion Initiation and Reward');
    subplot(2,1,1);
    histfit(initWavesWavelength,100,'kernel');
    xline(avgWavelengthInit,'-r',{'Mean Wavelength = ' num2str(avgWavelengthInit) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength : Motion Initiation');
    subplot(2,1,2);
    histfit(rewardWavesWavelength,100,'kernel');
    xline(avgWavelengthReward,'-r',{'Mean Wavelength = ' num2str(avgWavelengthReward) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wave Wavelength : Reward');

    figure('Name','Wavelengths in motion Initiation and Reward');
    group = [ones(size(initWavesWavelength')); 2.*ones(size(rewardWavesWavelength'))];
    boxplot([initWavesWavelength';rewardWavesWavelength'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Reward'});
    ylabel('Wave Wavelength in cm');
end


% Wave direction stats

avgDirInit = mean(initWavesDir);
avgDirReward = mean(rewardWavesDir);

[p, t] = ranksum(initWavesDir, rewardWavesDir);
% Print the results.
disp('Wave Direction')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Polar Histogram for wave direction in Motion Initiation and Reward');
    subplot(2,1,1);
    polarhistogram(initWavesDir,30);
    title('Wave Direction : Motion Initiation');
    subplot(2,1,2);
    polarhistogram(rewardWavesDir,30);
    title('Wave Direction : Reward');

    figure('Name','Histogram for wave direction in Motion Initiation and Reward');
    subplot(2,1,1);
    histogram(rad2deg(initWavesDir),72);
    title('Wave Direction : Motion Initiation');
    subplot(2,1,2);
    histogram(rad2deg(rewardWavesDir),72);
    title('Wave Direction : Reward');

    figure('Name','Wave Direction in motion Initiation and Reward');
    group = [ones(size(initWavesDir'));2.*ones(size(rewardWavesDir'))];
    boxplot([rad2deg(initWavesDir)';rad2deg(rewardWavesDir)'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Reward'});
end



wavesStat.initWavesIndx = initWavesIndx;
wavesStat.initWavesSpeed = initWavesSpeed;
wavesStat.initWavesDir = initWavesDir;
wavesStat.initWavesWavelength = initWavesWavelength;
wavesStat.initWavesSource = initWavesSource;

wavesStat.rewardWavesIndx = rewardWavesIndx;
wavesStat.rewardWavesSpeed = rewardWavesSpeed;
wavesStat.rewardWavesDir = rewardWavesDir;
wavesStat.rewardWavesWavelength = rewardWavesWavelength;
wavesStat.rewardWavesSource = rewardWavesSource;

wavesStat.nWavesInit = nWavesInit;
wavesStat.nWavesReward = nWavesReward;

wavesStat.avgSpeedInit = avgSpeedInit;
wavesStat.avgSpeedReward = avgSpeedReward;
wavesStat.avgWavelengthInit = avgWavelengthInit;
wavesStat.avgWavelengthReward = avgWavelengthReward;
wavesStat.avgDirInit = avgDirInit;
wavesStat.avgDirReward = avgDirReward;

end

