%% Binning Spikes 
if ~isfield('M1Spikes','spikes')
    disp('Bining M1 spikes ...')
    M1Spikes = binSpikes(M1Spikes,20000,1000);
end
if ~isfield('M2Spikes','spikes')
    disp('Bining M2 spikes ...')
    M2Spikes = binSpikes(M2Spikes,20000,1000);
end

if ~(M1Spikes.biningFs == M2Spikes.biningFs)
    disp('Bining Fs for M1 and M2 spikes is different')
end

%% Setting parameters for CCA Analysis 

% The units of the arguments are with respect to the binning window used to bin spikes.
argIn.BinWidth = 1000*(1/M1Spikes.biningFs);
argIn.MaxDelay = 0;  
argIn.TimeStep = 20;    
argIn.WindowLength = 80;
argIn.numCanonDim = 5;

argIn.NumWorkers = Inf; % Requires Parallel Processing Toolbox

%% Calculating CCA 
% Cue Hit
spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.cueHitTrace,0);
expCond = ones(size(spikes{1,1},3),1);

disp(argIn)
disp('Calculating CCA for cue Hit')
CCA.cueHit = ComputeCorrMap(spikes, expCond, argIn);

% Cue Miss
spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.cueMissTrace,0);
expCond = ones(size(spikes{1,1},3),1);

disp(argIn)
disp('Calculating CCA for cue Miss')
CCA.cueMiss = ComputeCorrMap(spikes, expCond, argIn);

% MI Hit
spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.MIHitTrace,0);
expCond = ones(size(spikes{1,1},3),1);

disp(argIn)
disp('Calculating CCA for MI Hit')
CCA.MIHit = ComputeCorrMap(spikes, expCond, argIn);

% MI FA
spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.MIFATrace,0);
expCond = ones(size(spikes{1,1},3),1);

disp(argIn)
disp('Calculating CCA for MI FA')
CCA.MIFA = ComputeCorrMap(spikes, expCond, argIn);

%% Plotting 
CANONICAL_PAIR_IDX = 1;
disp(["Plotting CCA for Dimension - " , CANONICAL_PAIR_IDX]);
mapDim = size(CCA.cueHit.CorrMap, 2);
delays = (-argIn.MaxDelay:argIn.MaxDelay); % Convert to ms
t = (1:argIn.TimeStep:argIn.TimeStep*mapDim); % Convert to ms

zeroDelayIndex = floor(numel(delays)/2) + 1;

figure,h1=plot(t,squeeze(CCA.cueHit.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX)),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
h2=plot(t,squeeze(CCA.cueMiss.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX)),'Color', [0.5 0.5 0.5],'LineWidth',2);
xline(1501,'--r','Cue');xline(1500+parameters.Fs*mean(IntanBehaviour.reactionTime,'all'),'--r','RT');
xlabel('Time (ms)'); ylabel('Population Correlation');
legend([h1 h2],'Hit','Miss','Location','best'); %ylim([5 15]);
title('CCA - Hit vs Miss')
xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');
drawnow;

figure,h1=plot(t,squeeze(CCA.MIHit.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX)),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
h2=plot(t,squeeze(CCA.MIFA.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX)),'Color', [0.9 0.1 0.1],'LineWidth',2);
xline(1501,'--r','MI');
xlabel('Time (ms)'); ylabel('Population Correlation');
legend([h1 h2],'Hit','FA','Location','best'); %ylim([5 15]);
title('CCA - MI Hit vs MI FA')
xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');
drawnow;

%% Ploting evoked CCA after post with respect to delay 
evokedPeriodStart = 1500;
evokedPeriodStop = 2000;
startIdx = find(t>=evokedPeriodStart,1);
stopIdx = find(t>=evokedPeriodStop,1);
CANONICAL_PAIR_IDX = 1;
figure;
subplot(2,1,1);
plot(delays,mean(CCA.cueHit.CorrMap(:,startIdx:stopIdx,CANONICAL_PAIR_IDX),2),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
xlabel('Delay (ms)');ylabel('Mean Evoked Population Correlation');
subplot(2,1,2);
plot(delays,max(CCA.cueHit.CorrMap(:,startIdx:stopIdx,CANONICAL_PAIR_IDX),[],2),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
xlabel('Delay (ms)');ylabel('Max Evoked Population Correlation')
drawnow;

%% Ploting CCA with delays 
CANONICAL_PAIR_IDX = 1;
mapDim = size(CCA.cueHit.CorrMap, 2);
delays = (-argIn.MaxDelay:argIn.MaxDelay); % Convert to ms
t = (1:argIn.TimeStep:argIn.TimeStep*mapDim); % Convert to ms

figure,
subplot(4,1,1)
imagesc(t,delays,CCA.cueHit.CorrMap(:,:,CANONICAL_PAIR_IDX));
ax = gca; ax.YDir = 'Normal';
ylabel('Delay');xlabel('Time')
subplot(4,1,2)
imagesc(t,delays,CCA.cueMiss.CorrMap(:,:,CANONICAL_PAIR_IDX));
ax = gca; ax.YDir = 'Normal';
ylabel('Delay');xlabel('Time')
subplot(4,1,3)
imagesc(t,delays,CCA.MIHit.CorrMap(:,:,CANONICAL_PAIR_IDX));
ax = gca; ax.YDir = 'Normal';
ylabel('Delay');xlabel('Time')
subplot(4,1,4)
imagesc(t,delays,CCA.MIFA.CorrMap(:,:,CANONICAL_PAIR_IDX));
ax = gca; ax.YDir = 'Normal';
ylabel('Delay');xlabel('Time')

%% Itterating over subset of trials
nItteration = 20;       % Number of interrations
nTrialFraction = 0.8;   % Fraction of selected trails for each itteration

% Parameters 
argIn.BinWidth = 1;
argIn.MaxDelay = 0;  
argIn.TimeStep = 40;    
argIn.WindowLength = 80;
argIn.NumWorkers = Inf; % Requires Parallel Processing Toolbox

t = (1:argIn.TimeStep:argIn.TimeStep*mapDim); % Convert to ms

% Cue Hit
CCAShuffled.cueHit.CorrMap = zeros(nItteration,numel(t),2);
disp(argIn)
for i=1:nItteration
    i
    trials = randperm(size(M1Spikes.CCASpikes.cueHit,3));
    spikes{1,1} = M1Spikes.CCASpikes.cueHit(:,:,trials(1:round(0.8*numel(trials))));
    spikes{1,2} = M2Spikes.CCASpikes.cueHit(:,:,trials(1:round(0.8*numel(trials))));
    expCond = ones(size(spikes{1,1},3),1);
    A = ComputeCorrMap(spikes, expCond, argIn);
    CCAShuffled.cueHit.CorrMap(i,:,:) = squeeze(A.CorrMap);
end

% Cue Miss
CCAShuffled.cueMiss.CorrMap = zeros(nItteration,numel(t),2);
for i=1:nItteration
    i
    trials = randperm(size(M1Spikes.CCASpikes.cueMiss,3));
    spikes{1,1} = M1Spikes.CCASpikes.cueMiss(:,:,trials(1:round(0.8*numel(trials))));
    spikes{1,2} = M2Spikes.CCASpikes.cueMiss(:,:,trials(1:round(0.8*numel(trials))));
    expCond = ones(size(spikes{1,1},3),1);
    disp(argIn)
    A = ComputeCorrMap(spikes, expCond, argIn);
    CCAShuffled.cueMiss.CorrMap(i,:,:) = squeeze(A.CorrMap);
end

% MI Hit
CCAShuffled.MIHit.CorrMap = zeros(nItteration,numel(t),2);
disp(argIn)
for i=1:nItteration
    i
    trials = randperm(size(M1Spikes.CCASpikes.MIHit,3));
    spikes{1,1} = M1Spikes.CCASpikes.MIHit(:,:,trials(1:round(0.8*numel(trials))));
    spikes{1,2} = M2Spikes.CCASpikes.MIHit(:,:,trials(1:round(0.8*numel(trials))));
    expCond = ones(size(spikes{1,1},3),1);
    A = ComputeCorrMap(spikes, expCond, argIn);
    CCAShuffled.MIHit.CorrMap(i,:,:) = squeeze(A.CorrMap);
end

% MI FA
CCAShuffled.MIFA.CorrMap = zeros(nItteration,numel(t),2);
disp(argIn)
for i=1:nItteration
    i
    trials = randperm(size(M1Spikes.CCASpikes.MIFA,3));
    spikes{1,1} = M1Spikes.CCASpikes.MIFA(:,:,trials(1:round(0.8*numel(trials))));
    spikes{1,2} = M2Spikes.CCASpikes.MIFA(:,:,trials(1:round(0.8*numel(trials))));
    expCond = ones(size(spikes{1,1},3),1);
    A = ComputeCorrMap(spikes, expCond, argIn);
    CCAShuffled.MIFA.CorrMap(i,:,:) = squeeze(A.CorrMap);
end

%% Plotting 
CANONICAL_PAIR_IDX = 1;
mapDim = size(CCA.cueHit.CorrMap, 2);
delays = (-argIn.MaxDelay:argIn.MaxDelay); % Convert to ms
t = (1:argIn.TimeStep:argIn.TimeStep*mapDim); % Convert to ms

zeroDelayIndex = floor(numel(delays)/2) + 1;
% Cue Hit vs Miss
figure,
avgCCA = mean(CCAShuffled.cueHit.CorrMap(:,:,CANONICAL_PAIR_IDX),1);
semCCA = std(CCAShuffled.cueHit.CorrMap(:,:,CANONICAL_PAIR_IDX),1)/sqrt(size(CCAShuffled.cueHit.CorrMap,1));

h1=plot(t,avgCCA,'Color', [0 0.1 0.8],'LineWidth',2); hold on;
plot(t,avgCCA-semCCA,'Color', [0 0.1 0.8],'LineWidth',1); hold on;
plot(t,avgCCA+semCCA,'Color', [0 0.1 0.8],'LineWidth',1); hold on;

avgCCA = mean(CCAShuffled.cueMiss.CorrMap(:,:,CANONICAL_PAIR_IDX),1);
semCCA = std(CCAShuffled.cueMiss.CorrMap(:,:,CANONICAL_PAIR_IDX),1)/sqrt(size(CCAShuffled.cueMiss.CorrMap,1));

h2=plot(t,avgCCA,'Color', [0.5 0.5 0.5],'LineWidth',2); hold on;
plot(t,avgCCA-semCCA,'Color', [0.5 0.5 0.5],'LineWidth',1); hold on;
plot(t,avgCCA+semCCA,'Color', [0.5 0.5 0.5],'LineWidth',1); hold on;

xline(1501,'--r','Cue');xline(1500+parameters.Fs*mean(IntanBehaviour.reactionTime,'all'),'--r','RT');
xlabel('Time (ms)'); ylabel('Population Correlation');
legend([h1 h2],'Hit','Miss','Location','best'); %ylim([5 15]);
title('CCA - Hit vs Miss')
xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');

% MI Hit vs MI FA
figure,
avgCCA = mean(CCAShuffled.MIHit.CorrMap(:,:,CANONICAL_PAIR_IDX),1);
semCCA = std(CCAShuffled.MIHit.CorrMap(:,:,CANONICAL_PAIR_IDX),1)/sqrt(size(CCAShuffled.MIHit.CorrMap,1));

h1=plot(t,avgCCA,'Color', [0 0.1 0.8],'LineWidth',2); hold on;
plot(t,avgCCA-semCCA,'Color', [0 0.1 0.8],'LineWidth',1); hold on;
plot(t,avgCCA+semCCA,'Color', [0 0.1 0.8],'LineWidth',1); hold on;

avgCCA = mean(CCAShuffled.MIFA.CorrMap(:,:,CANONICAL_PAIR_IDX),1);
semCCA = std(CCAShuffled.MIFA.CorrMap(:,:,CANONICAL_PAIR_IDX),1)/sqrt(size(CCAShuffled.MIFA.CorrMap,1));

h2=plot(t,avgCCA,'Color', [0.9 0.1 0.1],'LineWidth',2); hold on;
plot(t,avgCCA-semCCA,'Color', [0.9 0.1 0.1],'LineWidth',1); hold on;
plot(t,avgCCA+semCCA,'Color', [0.9 0.1 0.1],'LineWidth',1); hold on;

xline(1501,'--r','MI');
xlabel('Time (ms)'); ylabel('Population Correlation');
legend([h1 h2],'MI Hit','MI FA','Location','best'); %ylim([5 15]);
title('CCA - Hit vs FA')
xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');


