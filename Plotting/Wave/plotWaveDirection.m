function [dirComb, fittedVm] = plotWaveDirection(Waves,varargin)

if isempty(varargin{1})
    fit = 1;
    disp('Fitting mixture von mises distributions')
else
    fit = 0;
    disp('Fitted distribution passed')
end

dirComb = horzcat(Waves(1:end).waveDir);
angles = linspace(-pi, pi, 1000)';

if fit == 0
    fittedVm = varargin{1};
    nComponents = fittedVm.nComponents;
    likelihoods = fittedVm.pdf(angles);
else
    [~, fittedVm] = getClusterDirectionality2(Waves,0.95,3);
    nComponents = fittedVm.nComponents;
    likelihoods = fittedVm.pdf(angles);
end
polarhistogram(dirComb,72,'Normalization','pdf','FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.3,'EdgeAlpha',0.3);hold on;
polarplot(angles,likelihoods,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);%rlim([0 0.6]);
if nComponents ~= 0 polarplot(repmat(fittedVm.mu',2,1),[zeros(1,nComponents);fittedVm.pdf(fittedVm.mu)'],'-r','LineWidth',1.5), end
title('Wave Direction');box off;

end

% 
% blue - [153/255 153/255 255/255]
% green - [0.4660 0.6740 0.1880]