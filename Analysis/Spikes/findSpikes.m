function [Spikes] = findSpikes(Spikes,thres)
%  An amplitude threshold window
% was set 3.5 standard deviations above and below the mean of
% the sample distribution. For each peak exceeding the threshold
% window, a 2.4 ms candidate waveform snippet centered on the
% absolute minimum of the waveform was removed from the
% recorded segment and stored
% Ref - K A Ludwig et al, 2006, Journal of Neural Engineering

Spikes.threshold = zeros(size(Spikes.hpSpikes,1),1);
Spikes.threshold = thres*std(Spikes.hpSpikes,0,2);

end

