function [MI] = getMutualInformation(xHit,xMiss,parameters)
% Reference : Propagating waves mediate information transfer in the motor
% cortex, Nature Neuro, 2006 
   
totaltime = parameters.windowAfterCue + parameters.windowBeforeCue + parameters.ts;
time = [parameters.ts:parameters.ts:totaltime]; 

nHits = size(xHit,2);
nMiss = size(xMiss,2);

Nbins = ceil(exp(0.626 + (0.4*log(nHits+nMiss-1)))); % number of bins
% Ref - Comparison of Hilbert transform and wavelet methods for the
%       analysis of neuronal synchrony, JNeuroMethods, 2001

HX = zeros(parameters.rows,parameters.cols,size(time,2)); % Shannon entropy 
HXY = zeros(parameters.rows,parameters.cols,size(time,2)); % conditonal shannon entropy 
MI = zeros(parameters.rows,parameters.cols,size(time,2)); % Mutual information 

for t=1:size(time,2)
    for i=1:parameters.rows
        for j=1:parameters.cols
            xijtHit = cellfun(@(s) s(i,j,t),xHit);
            xijtMiss = cellfun(@(s) s(i,j,t),xMiss);
            xijtAll = horzcat(xijtHit,xijtMiss);
            HX(i,j,t) = getShannonEntropy(xijtAll,Nbins,-pi,pi);
            HXY(i,j,t) = getCondShannonEntropy(xijtHit,xijtMiss,Nbins,-pi,pi);
            MI(i,j,t) = HX(i,j,t) - HXY(i,j,t);
        end
    end
end
