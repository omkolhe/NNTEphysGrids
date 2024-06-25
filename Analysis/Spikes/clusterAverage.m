function [avgTempProp] = clusterAverage(spikeTemplate,spikeProperty)

templates = unique(spikeTemplate);

for i=1:numel(templates)
    indx = find(spikeTemplate==templates(i));
    avgTempProp(i) = mean(spikeProperty(indx));
end
