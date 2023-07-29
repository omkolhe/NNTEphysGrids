function [BetaBurst] = detectBetaBurst(LFP,behaviourTrace,parameters)

% Reference - https://academic.oup.com/brain/article/140/11/2968/4430808#112576655

betaBurstThres = 75; % 75 percenticle 
xamp = abs(LFP.xgp);

%% Getting LFP for each trial and channels 
rows = size(LFP,1);
cols = size(LFP,2);
nlength = size(behaviourTrace.trace,1);
nTrials = size(behaviourTrace,2);

xTrials = zeros(rows,cols,nlength);

for i=1:nTrials
    xTrials(:,:,i) = xamp(:,:,behaviourTrace.LFPIndex(1):behaviourTrace.LFPIndex(end));
end

%% Getting thresholds


end

