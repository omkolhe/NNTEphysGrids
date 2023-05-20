function [Encoder] = getVelTrigTrials(Encoder,windowBeforeTrig,windowAfterTrig,LFPFs)

rejectedTrigs = [];

for i=1:size(Encoder.velTrig,2)
    startTime1 = Encoder.velTrig(1,i)-windowBeforeTrig*Encoder.fs;
    startTime2 = Encoder.velTrig(2,i)-windowBeforeTrig*LFPFs;
    if startTime2<0 
        %         Encoder.newTrig(i) = NaN;
        Encoder.trialTime(i,1) = NaN;
        Encoder.trialTime(i,2) = NaN;
        Encoder.trialTime(i,3) = NaN;
        Encoder.trialTime(i,4) = NaN;
        rejectedTrigs = [rejectedTrigs i];
        continue;

    end
    stopTime1 = Encoder.velTrig(1,i)+windowAfterTrig*Encoder.fs - 1;
    stopTime2 = Encoder.velTrig(2,i)+windowAfterTrig*LFPFs - 1;
    %     Encoder.newTrig(i) = Encoder.velTrig(i);
    if(mean(Encoder.vel(startTime1:stopTime1-windowAfterTrig*Encoder.fs))> 2)
        Encoder.trialTime(i,1) = NaN;
        Encoder.trialTime(i,2) = NaN;
        Encoder.trialTime(i,3) = NaN;
        Encoder.trialTime(i,4) = NaN;
        rejectedTrigs = [rejectedTrigs i];
    elseif(mean(Encoder.vel(stopTime1-windowAfterTrig*Encoder.fs:stopTime1+0.0*Encoder.fs))<4)
        Encoder.trialTime(i,1) = NaN;
        Encoder.trialTime(i,2) = NaN;
        Encoder.trialTime(i,3) = NaN;
        Encoder.trialTime(i,4) = NaN;
        rejectedTrigs = [rejectedTrigs i];
    else
        Encoder.trialTime(i,1) = startTime1;
        Encoder.trialTime(i,2) = stopTime1;
        Encoder.trialTime(i,3) = startTime2;
        Encoder.trialTime(i,4) = stopTime2;
    end
end
% Encoder.goodInd = Encoder.goodInd(sum(isnan(Encoder.trialTime),2)==0,:);
Encoder.trialTime = Encoder.trialTime(sum(isnan(Encoder.trialTime),2)==0,:);
Encoder.nTrig = size(Encoder.trialTime,1); 
Encoder.velTrig(:,rejectedTrigs) = [];
end


