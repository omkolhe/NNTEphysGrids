%% Checking for motion artificats 
meanLFPPower = arrayfun(@(s) s.xf, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
meanLFPPower = cat(3,meanLFPPower{:});
meanLFPPower = mean(meanLFPPower.^2,3);

nBadTrials = 0;
for i=1:size(IntanBehaviour.cueHitTrace,2)
    nMotionChs = 0;
    for j=1:parameters.rows*parameters.cols
        k = chToGrid(j,parameters);
        a = find((IntanBehaviour.cueHitTrace(i).xf(k(1),k(2),:)).^2 > 15*meanLFPPower(k(1),k(2)));
        if numel(a)>10
            nMotionChs = nMotionChs+1;
        end
    end
    if nMotionChs > 4
        IntanBehaviour.cueHitTrace(i).motion = 1;
        nBadTrials = nBadTrials+1;
    else
        IntanBehaviour.cueHitTrace(i).motion = 0;
    end
end