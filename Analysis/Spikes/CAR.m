function [Spikes] = CAR(Spikes)
%common average referencing 

    Spikes.commAvgRef = mean(Spikes.hpSpikes,1);
    Spikes.hpSpikes = Spikes.hpSpikes - Spikes.commAvgRef;
end

