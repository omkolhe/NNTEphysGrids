%% Setting parameters for CCA Analysis 

% The units of the arguments are with respect to the binning window used to bin spikes.
argIn.BinWidth = 1;
argIn.MaxDelay = 0;  
argIn.TimeStep = 5;    
argIn.WindowLength = 50;

argIn.NumWorkers = Inf; % Requires Parallel Processing Toolbox

%% Calculating CCA 
maxLag = 40; % in ms
lagStep = 2; % in ms

% Cue Hit
disp(argIn)
disp('Calculating CCA for cue Hit')
i = 1;
for M2M1Lag = -maxLag:lagStep:maxLag
    spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.cueHitTrace,M2M1Lag);
    expCond = ones(size(spikes{1,1},3),1);
    A = ComputeCorrMap(spikes, expCond, argIn);
    CCADelay.cueHit(i,:,:) = squeeze(A.CorrMap);
    i=i+1;
end

% Cue Miss
disp(argIn)
disp('Calculating CCA for cue Miss')
i = 1;
for M2M1Lag = -maxLag:lagStep:maxLag
    spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.cueMissTrace,M2M1Lag);
    expCond = ones(size(spikes{1,1},3),1);
    A = ComputeCorrMap(spikes, expCond, argIn);
    CCADelay.cueMiss(i,:,:) = squeeze(A.CorrMap);
    i=i+1;
end

% MI Hit
disp(argIn)
disp('Calculating CCA for MI Hit')
i = 1;
for M2M1Lag = -maxLag:lagStep:maxLag
    spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.MIHitTrace,M2M1Lag);
    expCond = ones(size(spikes{1,1},3),1);
    A = ComputeCorrMap(spikes, expCond, argIn);
    CCADelay.MIHit(i,:,:) = squeeze(A.CorrMap);
    i=i+1;
end

% MI FA
disp(argIn)
disp('Calculating CCA for MI FA')
i = 1;
for M2M1Lag = -maxLag:lagStep:maxLag
    spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.MIFATrace,M2M1Lag);
    expCond = ones(size(spikes{1,1},3),1);
    A = ComputeCorrMap(spikes, expCond, argIn);
    CCADelay.MIFA(i,:,:) = squeeze(A.CorrMap);
    i=i+1;
end


%% Plotting
CANONICAL_PAIR_IDX = 1;
mapDim = size(CCADelay.cueHit, 2);
delays = -maxLag:lagStep:maxLag; % Convert to ms
t = (1:argIn.TimeStep:argIn.TimeStep*mapDim); % Convert to ms

figure;
ax1 = subplot(4,1,1);
imagesc(t,delays,CCADelay.cueHit(:,:,CANONICAL_PAIR_IDX));
ax = gca; ax.YDir = 'Normal';
ylabel('Delay');xlabel('Time');title('Cue Hit');
ax2 = subplot(4,1,2);
imagesc(t,delays,CCADelay.cueMiss(:,:,CANONICAL_PAIR_IDX));
ax = gca; ax.YDir = 'Normal';
ylabel('Delay');xlabel('Time');title('Cue Miss');
ax3 = subplot(4,1,3);
imagesc(t,delays,CCADelay.MIHit(:,:,CANONICAL_PAIR_IDX));
ax = gca; ax.YDir = 'Normal';
ylabel('Delay');xlabel('Time');title('MI Hit');
ax4 = subplot(4,1,4);
imagesc(t,delays,CCADelay.MIFA(:,:,CANONICAL_PAIR_IDX));
ax = gca; ax.YDir = 'Normal';
ylabel('Delay');xlabel('Time');title('MI FA');
linkaxes([ax1,ax2,ax3,ax4], 'xy');


%% Plotting evoked 

figure,
plot(delays,mean(CCADelay.cueHit(:,250:350,CANONICAL_PAIR_IDX),2));