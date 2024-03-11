function [dirComb, fittedVm] = plotWaveDirection(Waves)

dirComb = horzcat(Waves(1:end).waveDir);

% Fit VonMisesMixture Model 
nComponents = 2;    
fittedVm = fitmvmdist(dirComb',nComponents,'MaxIter', 250);
angles = linspace(-pi, pi, 1000)';
likelihoods = fittedVm.pdf(angles);

polarhistogram(dirComb,36,'Normalization','pdf','FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.3,'EdgeAlpha',0.3);hold on;
polarplot(angles,likelihoods,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);%rlim([0 0.6]);
polarplot([fittedVm.mu';fittedVm.mu'],[0 0;fittedVm.pdf(fittedVm.mu)'],'-r','LineWidth',1.5);
title('Wave Direction');box off;

end

% 
% blue - [153/255 153/255 255/255]
% green - [0.4660 0.6740 0.1880]