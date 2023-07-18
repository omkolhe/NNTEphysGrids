function [PPL] = getPPL(xgp,behaviourTrace,parameters)

% Reference : Propagating waves mediate information transfer in the motor
% cortex, Nature Neuro, 2006 
   
totaltime = parameters.windowAfterPull + parameters.windowBeforePull + parameters.ts;
time = [parameters.ts:parameters.ts:totaltime]; 

nTrials = size(behaviourTrace,2);
allXGP = zeros(nTrials ,size(time,2));
for trialno=1:nTrials 
    allXGP(trialno,:) = angle(xgp(behaviourTrace(trialno).LFPIndex(1):behaviourTrace(trialno).LFPIndex(end)));
end

Hmax = log2(nTrials);

H = zeros(parameters.rows,parameters.cols,size(time,2));
PPL = zeros(parameters.rows,parameters.cols,size(time,2));

for t=1:size(time,2)
    for i=1:parameters.rows
        for j=1:parameters.cols
            H(i,j,t) = getShannonEntropy(squeeze(allXGP(:,t)),6,-pi,pi);
            PPL(i,j,t) = 100*(1-(H(i,j,t)/Hmax));
        end
    end
end

