function [PPL] = getPPL(xgp,behaviourTrace,parameters)

% Reference : Propagating waves mediate information transfer in the motor
% cortex, Nature Neuro, 2006 
   
totaltime = parameters.windowAfterCue + parameters.windowBeforeCue + parameters.ts;
time = [parameters.ts:parameters.ts:totaltime]; 

nTrials = size(behaviourTrace,2);
allXGP = zeros(nTrials ,size(time,2));
for trialno=1:nTrials 
    allXGP(trialno,:) = angle(xgp(behaviourTrace(trialno).LFPIndex(1):behaviourTrace(trialno).LFPIndex(end)));
end

N = ceil(exp(0.626 + (0.4*log(nTrials-1)))); % number of bins
% Ref - Comparison of Hilbert transform and wavelet methods for the
%       analysis of neuronal synchrony, JNeuroMethods, 2001

Hmax = log2(N);

H = zeros(parameters.rows,parameters.cols,size(time,2));
PPL = zeros(parameters.rows,parameters.cols,size(time,2));

for t=1:size(time,2)
    for i=1:parameters.rows
        for j=1:parameters.cols
            H(i,j,t) = getShannonEntropy(squeeze(allXGP(:,t)),N,-pi,pi);
            PPL(i,j,t) = 100*(1-(H(i,j,t)/Hmax));
        end
    end
end

