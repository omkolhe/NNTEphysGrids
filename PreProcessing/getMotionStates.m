function [Encoder] = getMotionStates(Encoder,LFPTime)

% 0 - Rest
% 1 - Initiation
% 2 - Running 
% 3 - Termination

Encoder.acc = zeros(1,numel(Encoder.time));
Encoder.acc(2:end) = diff(Encoder.vel)*Encoder.fs;

% Discretizing velocity into 0.2 cm/s 
maxVel = max(Encoder.vel);
segmentSize = 0.2;
fillFindingGap = 1*Encoder.fs;
velEdges = -segmentSize:segmentSize:maxVel+segmentSize;
Encoder.discVel = (discretize(Encoder.vel,velEdges)-1)*segmentSize;

VelTrigLow = 0.5;  % Lower threshold for rest in cm/s
VelTrigHigh = 8; % Higher threshold for running in cm/s

Encoder.state = zeros(size(Encoder.discVel,2),1);
Encoder.state(Encoder.discVel>VelTrigHigh) = 2;
Encoder.state(Encoder.discVel<=VelTrigHigh & Encoder.discVel >VelTrigLow) = 1;

runindx = find(Encoder.state == 2);
diffrunindx = zeros(numel(runindx),1);
diffrunindx(2:end) = diff(runindx);
rungapfill = find(diffrunindx<=fillFindingGap & diffrunindx>1);
for i=1:numel(rungapfill)
    Encoder.state(runindx(rungapfill(i)-1):runindx(rungapfill(i))) = 2;
end

% Segmenting state = 1 into initiation and termination
motionindx = find(Encoder.state == 1);
diffmotionindx = zeros(numel(motionindx),1);
diffmotionindx(2:end) = diff(motionindx);
motiongap = find(diffmotionindx > 1);
for i = 1:numel(motiongap)
    if i==1
        st = motionindx(1);
        sp = motionindx(motiongap(1)-1);
    else
        st = motionindx(motiongap(i-1));
        sp = motionindx(motiongap(i)-1);
    end
    beforeState = Encoder.state(st-1);
    afterState = Encoder.state(sp+1);
    
    if(beforeState == 0 && afterState == 2)
        Encoder.state(st:sp) = 1;
    elseif beforeState == 2 && afterState == 0
            Encoder.state(st:sp) = 3;
    elseif beforeState == 0 && afterState == 0
            Encoder.state(st:sp) = 0;
    end
end

% Getting trial times for each state
diffstate = zeros(numel(Encoder.state),1);
diffstate(2:end) = diff(Encoder.state);
diffstateindx = find(diffstate ~= 0);
Encoder.restTrialTime = [];
Encoder.initTrialTime = [];
Encoder.runTrialTime = [];
Encoder.termTrialTime = [];

for i=1:numel(diffstateindx)
    if i==1
        st = 1;
        sp = diffstateindx(1)-1;
    else
        st = diffstateindx(i-1);
        sp = diffstateindx(i)-1;
    end

    if(sp-st < 1*Encoder.fs)
        continue;
    end

    stLFP = max(find(LFPTime<=Encoder.time(st)));
    spLFP = max(find(LFPTime<=Encoder.time(sp)));

    stateestimate = mean(Encoder.state(st:sp),'all');

    if stateestimate <= 0.5
        Encoder.restTrialTime = [Encoder.restTrialTime; [st sp stLFP spLFP]];
    elseif stateestimate >= 0.5 && stateestimate <= 1.5 
        Encoder.initTrialTime = [Encoder.initTrialTime; [st sp stLFP spLFP]];
    elseif stateestimate >= 1.5 && stateestimate <= 2.5 
        Encoder.runTrialTime = [Encoder.runTrialTime; [st sp stLFP spLFP]];
    else
        Encoder.termTrialTime = [Encoder.termTrialTime; [st sp stLFP spLFP]];
    end
end

end

