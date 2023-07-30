function [BetaEvent] = detectBetaEvent(LFP,behaviourTrace,parameters)

% Reference - https://academic.oup.com/brain/article/140/11/2968/4430808#112576655

betaBurstThres = 2.5; % 2.5 * std dev 
betaAttenThres = 2.5; % -2.5 * std dev
minInterRippleInterval = 30; % 30ms
minBetaDuration = 20; % 50ms
maxBetaDuration = 250; % 200ms

%% Getting LFP for each trial and channels 
rows = size(LFP,1);
cols = size(LFP,2);
nlength = size(behaviourTrace(2).trace,2);
nTrials = size(behaviourTrace,2);

for i=1:nTrials
    BetaEvent(i).xamp = squeeze(mean(LFP(:,:,behaviourTrace(i).LFPIndex(1):behaviourTrace(i).LFPIndex(end)),[1,2]));
    BetaEvent(i).xpow = BetaEvent(i).xamp.^2;
end

%% Normalizing and getting thresholds
normByTrial = 1; % 1 - normalize power by noise in each trial, 0 - normalize power by noise across all trials

if normByTrial == 1
    for i=1:nTrials
        BetaEvent(i).xnorm = (BetaEvent(i).xpow - mean(BetaEvent(i).xpow,'all'))/std(BetaEvent(i).xpow);
    end
else
    meanPower = mean(horzcat(BetaEvent(1:end).xpow),'all');
    stdPower = std(horzcat(BetaEvent(1:end).xpow),'all');

end


%% Getting beta burst segments 
for i=1:nTrials
    BetaBurst(i).burstsegments = findThresSeg(xTrials(i).xamp,betaThres);
end

function normSig = normSignal(xpow)
    norm