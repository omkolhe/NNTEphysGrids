function Spikes = binSpikes(Spikes,spikeFs,biningFs)
%% Generate binned spikes from M1/M2Spikes and IntanBehaviour 
% Creating spiking matrix - N x t 
%      N - Number of neurons 
%      t - total time 

allSpikeTime = cell2mat((arrayfun(@(s) s.spikeTime, Spikes.Clusters,"UniformOutput",false))');
maxTime = max(allSpikeTime);
maxIndex = round(maxTime*spikeFs);
clear allSpikeTime;

% Binning data 
numBins = floor((maxIndex/spikeFs)*1000);
nPointsBin = spikeFs/biningFs;
nExtra = maxIndex - numBins*nPointsBin;
ind = reshape(repmat(1:1:numBins,nPointsBin,1),[],1)';
ind = [ind ones(1,round(nExtra))*(numBins(end)+1)];

% Binning M1 Spikes
Spikes.allSpikes = zeros(Spikes.nSpikes,maxIndex,'uint8');
for i=1:Spikes.nSpikes   % number of neurons
    Spikes.allSpikes(i,Spikes.Clusters(i).cluster) = 1;
end

Spikes.spikes = zeros(Spikes.nSpikes,ind(end),'uint8');
for i=1:Spikes.nSpikes
    Spikes.spikes(i,:) = accumarray(ind',Spikes.allSpikes(i,:)')';
end
Spikes = rmfield(Spikes,"allSpikes");
Spikes.biningFs = biningFs;

end