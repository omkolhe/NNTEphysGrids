%% Plotting trajectories wrt time 

% Cue hit trajectories 
figure;
nDim = 3;
for dimension=1:nDim
    subplot(nDim,2,2*(dimension-1)+1)
    for i=1:size(neuralDynamics.hit.X,3)
        plot(neuralDynamics.time,squeeze(neuralDynamics.hit.X(dimension,:,i)),'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(neuralDynamics.time,squeeze(mean(neuralDynamics.hit.X(dimension,:,:),3)),'Color',[1 0 0 1],'LineWidth',2);
    xlabel('Time');ylabel(['Dimension' string(dimension)] );xline(1500,'--r','Cue','LabelVerticalAlignment','top');
    xline(1500+mean(IntanBehaviour.reactionTime,'all')*IntanBehaviour.parameters.Fs,'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    box off;
    subplot(nDim,2,2*dimension)
    for i=1:size(neuralDynamics.hit.speed.vel,3)
        plot(neuralDynamics.time(2:end),squeeze(neuralDynamics.hit.speed.vel(dimension,:,i)),'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(neuralDynamics.time(2:end),squeeze(mean(neuralDynamics.hit.speed.vel(dimension,:,:),3)),'Color',[1 0 0 1],'LineWidth',2);
    xlabel('Time');ylabel(['Velocity - Dimension' string(dimension)] );xline(1500,'--r','Cue','LabelVerticalAlignment','top');
    xline(1500+mean(IntanBehaviour.reactionTime,'all')*IntanBehaviour.parameters.Fs,'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    box off;
end

% MI hit trajectories 
figure;
nDim = 3;
for dimension=1:nDim
    subplot(nDim,2,2*(dimension-1)+1)
    for i=1:size(neuralDynamics.MIhit.X,3)
        plot(neuralDynamics.time,squeeze(neuralDynamics.MIhit.X(dimension,:,i)),'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(neuralDynamics.time,squeeze(mean(neuralDynamics.MIhit.X(dimension,:,:),3)),'Color',[1 0 0 1],'LineWidth',2);
    xlabel('Time');ylabel(['Dimension' string(dimension)] );xline(1500,'--r','MI','LabelVerticalAlignment','top');
    xline(1500+mean(IntanBehaviour.reactionTime,'all')*IntanBehaviour.parameters.Fs,'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    box off;
    subplot(nDim,2,2*dimension)
    for i=1:size(neuralDynamics.MIhit.speed.vel,3)
        plot(neuralDynamics.time(2:end),squeeze(neuralDynamics.MIhit.speed.vel(dimension,:,i)),'Color',[0 0 0 0.1],'LineWidth',1.5);
        hold on;
    end
    plot(neuralDynamics.time(2:end),squeeze(mean(neuralDynamics.MIhit.speed.vel(dimension,:,:),3)),'Color',[1 0 0 1],'LineWidth',2);
    xlabel('Time');ylabel(['Veloscity - Dimension' string(dimension)] );xline(1500,'--r','MI','LabelVerticalAlignment','top');
    xline(1500+mean(IntanBehaviour.reactionTime,'all')*IntanBehaviour.parameters.Fs,'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    box off;
end
%% Plotting trajectories in 3D space 
figure;
plot3d(neuralDynamics.hit.X(:,1))