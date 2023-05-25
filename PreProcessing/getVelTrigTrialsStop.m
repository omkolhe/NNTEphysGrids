function [Encoder] = getVelTrigTrialsStop(Encoder,windowBeforeTrig,windowAfterTrig,LFPFs)

rejectedTrigs = [];

for i=1:size(Encoder.velTrigStop,2)
    startTime1 = Encoder.velTrigStop(1,i)-windowBeforeTrig*Encoder.fs;
    startTime2 = Encoder.velTrigStop(2,i)-windowBeforeTrig*LFPFs;
    if startTime2<0 
        %         Encoder.newTrig(i) = NaN;
        Encoder.trialTimeStop(i,1) = NaN;
        Encoder.trialTimeStop(i,2) = NaN;
        Encoder.trialTimeStop(i,3) = NaN;
        Encoder.trialTimeStop(i,4) = NaN;
        rejectedTrigs = [rejectedTrigs i];
        continue;

    end
    stopTime1 = Encoder.velTrigStop(1,i)+windowAfterTrig*Encoder.fs - 1;
    stopTime2 = Encoder.velTrigStop(2,i)+windowAfterTrig*LFPFs - 1;
    %     Encoder.newTrig(i) = Encoder.velTrigStop(i);
    if(mean(Encoder.vel(startTime1:stopTime1-windowAfterTrig*Encoder.fs))< 5)
        Encoder.trialTimeStop(i,1) = NaN;
        Encoder.trialTimeStop(i,2) = NaN;
        Encoder.trialTimeStop(i,3) = NaN;
        Encoder.trialTimeStop(i,4) = NaN;
        rejectedTrigs = [rejectedTrigs i];
    elseif(mean(Encoder.vel(stopTime1-windowAfterTrig*Encoder.fs:stopTime1+0.0*Encoder.fs))>2)
        Encoder.trialTimeStop(i,1) = NaN;
        Encoder.trialTimeStop(i,2) = NaN;
        Encoder.trialTimeStop(i,3) = NaN;
        Encoder.trialTimeStop(i,4) = NaN;
        rejectedTrigs = [rejectedTrigs i];
    else
        Encoder.trialTimeStop(i,1) = startTime1;
        Encoder.trialTimeStop(i,2) = stopTime1;
        Encoder.trialTimeStop(i,3) = startTime2;
        Encoder.trialTimeStop(i,4) = stopTime2;
    end
end
% Encoder.goodInd = Encoder.goodInd(sum(isnan(Encoder.trialTimeStop),2)==0,:);
Encoder.trialTimeStop = Encoder.trialTimeStop(sum(isnan(Encoder.trialTimeStop),2)==0,:);
Encoder.nTrigStop = size(Encoder.trialTimeStop,1); 
Encoder.velTrigStop(:,rejectedTrigs) = [];
end


