plotFlag = 0;

for filenumber = 1:size(M1neuralDynamics,2)
    filenumber
    IntanBehaviour = M1neuralDynamics(filenumber).IntanBehaviour;
    neuralDynamics = M1neuralDynamics(filenumber).neuralDynamics;
    parameters = IntanBehaviour.parameters;
    rt = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false));
    if mean(rt) < 0
        warning('Bad file')
        M1neuralDynamics(filenumber).skip = 1;
        continue
    else
        M1neuralDynamics(filenumber).skip = 0;
    end
    %% Matching Neural Dynamics to the IntanBehaviour file 
    % Removing disengage miss trials, hit trials with RT>1.5 and RT<0,
    % overlaping Hit and FA trials
    if ~isfield(IntanBehaviour,'cleanFlag')
        IntanBehaviour = cleanIntanBehaviour(IntanBehaviour);
    end
    if ~isfield(neuralDynamics,'matchFlag')
        % Removing the same trials as above from the neural trajectory data
        neuralDynamics = matchBehaviourNeuralDynamics(neuralDynamics,IntanBehaviour);
    end
    % Calculating neural trajectory speed
    neuralDynamics = calTrajSpeed(neuralDynamics,IntanBehaviour,IntanBehaviour.parameters);
    
    if plotFlag == 1
        % Plotting Lever Trace and the neural trajectory speed - Cue alligned 
        figure;
        subplot(2,2,1);
        for i=1:size(neuralDynamics.hit.speed.speed,2)
            plot(neuralDynamics.time,neuralDynamics.hit.speed.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(neuralDynamics.time,mean(neuralDynamics.hit.speed.speed,2),'Color',[1 0 0 1],'LineWidth',2);
        xline(1500,'--r','Cue','LabelVerticalAlignment','top');
        xline(1500+mean(IntanBehaviour.reactionTime,'all')*IntanBehaviour.parameters.Fs,'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        xlabel('Time (in s)');box off;
        subplot(2,2,3);
        for i=1:size(IntanBehaviour.cueHitTrace,2)
            plot(IntanBehaviour.cueHitTrace(i).time,IntanBehaviour.cueHitTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(IntanBehaviour.cueHitTrace(1).time,mean(horzcat(IntanBehaviour.cueHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
        yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left');
        xline(0,'--r','Cue','LabelVerticalAlignment','top');
        xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Hits');box off;
        
        subplot(2,2,2);
        for i=1:size(neuralDynamics.miss.speed.speed,2)
            plot(neuralDynamics.time,neuralDynamics.miss.speed.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(neuralDynamics.time,mean(neuralDynamics.miss.speed.speed,2),'Color',[1 0 0 1],'LineWidth',2);
        xline(1500,'--r','Cue','LabelVerticalAlignment','top');
        xline(1500+mean(IntanBehaviour.reactionTime,'all')*IntanBehaviour.parameters.Fs,'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        xlabel('Time (in s)');box off;
        subplot(2,2,4);
        for i=1:size(IntanBehaviour.cueMissTrace,2)
            plot(IntanBehaviour.cueMissTrace(i).time,IntanBehaviour.cueMissTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(IntanBehaviour.cueMissTrace(1).time,mean(horzcat(IntanBehaviour.cueMissTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
        yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left');
        xline(0,'--r','Cue','LabelVerticalAlignment','top');
        xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Misss');box off;
        
        % Plotting Lever Trace and the neural trajectory speed - MI alligned 
        figure;
        subplot(2,2,1);
        for i=1:size(neuralDynamics.MIhit.speed.speed,2)
            plot(neuralDynamics.time,neuralDynamics.MIhit.speed.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(neuralDynamics.time,mean(neuralDynamics.MIhit.speed.speed,2),'Color',[1 0 0 1],'LineWidth',2);
        xline(1500,'--r','MI','LabelVerticalAlignment','top');
        xlabel('Time (in s)');box off;ylabel('Trajectory speed');
        subplot(2,2,3);
        for i=1:size(IntanBehaviour.MIHitTrace,2)
            plot(IntanBehaviour.MIHitTrace(i).time,IntanBehaviour.MIHitTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(IntanBehaviour.MIHitTrace(1).time,mean(horzcat(IntanBehaviour.MIHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
        yline(IntanBehaviour.MIcutoffHit,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 
        yline(IntanBehaviour.threshold,'--.b','Reward Threshold','LabelHorizontalAlignment','left'); 
        xline(0,'--r','MI','LabelVerticalAlignment','top');
        ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for MI Hits');box off;
        
        subplot(2,2,2);
        for i=1:size(neuralDynamics.MIFA.speed.speed,2)
            plot(neuralDynamics.time,neuralDynamics.MIFA.speed.speed(:,i),'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(neuralDynamics.time,mean(neuralDynamics.MIFA.speed.speed,2),'Color',[1 0 0 1],'LineWidth',2);
        xline(1500,'--r','MI','LabelVerticalAlignment','top');
        xlabel('Time (in s)');box off;ylabel('Trajectory speed');
        subplot(2,2,4);
        for i=1:size(IntanBehaviour.MIFATrace,2)
            plot(IntanBehaviour.MIFATrace(i).time,IntanBehaviour.MIFATrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
            hold on;
        end
        plot(IntanBehaviour.MIFATrace(1).time,mean(horzcat(IntanBehaviour.MIFATrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
        yline(IntanBehaviour.MIcutoffHit,'--.b','MI Threshold','LabelHorizontalAlignment','left'); 
        yline(IntanBehaviour.threshold,'--.b','Reward Threshold','LabelHorizontalAlignment','left'); 
        xline(0,'--r','MI','LabelVerticalAlignment','top');
        ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for MI FA');box off;
    end
    %% Correlating behaviuour to the trajectory speed
    % Getting the peak speed from cue to MI in cue aligned hit traces 
    % Cue index will be cueIndex 
    % MI index will be MIIndex
    for i=1:size(IntanBehaviour.cueHitTrace,2)
        [neuralDynamics.hit.speed.cueMIPeakSpeed(i),neuralDynamics.hit.speed.cueMIPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(neuralDynamics.cueIndex:neuralDynamics.MIIndex(i),i));
%         [neuralDynamics.hit.speed.cueMIPeakSpeed(i),neuralDynamics.hit.speed.cueMIPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(neuralDynamics.cueIndex:end,i));
%         [neuralDynamics.hit.speed.cueMIPeakSpeed(i),neuralDynamics.hit.speed.cueMIPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(neuralDynamics.cueIndex:neuralDynamics.RewardIndex(i),i));
    end
    
    % Getting the peak speed from cue to MI in MIHit traces 
    % Cue index will be cueIndex-(MIIndex-cueIndex) = 2*cueIndex-MIIndex
    % MI index will be cueIndex (midpoint)
    for i=1:size(IntanBehaviour.MIHitTrace,2)
        a = 2*neuralDynamics.cueIndex-neuralDynamics.MIIndex(i);
        if a < 1
            a = 1;
        end
        [neuralDynamics.MIhit.speed.cueMIPeakSpeed(i),neuralDynamics.MIhit.speed.cueMIPeakSpeedTime(i)] = max(neuralDynamics.hit.speed.speed(a:neuralDynamics.cueIndex,i));
    end
    
    % Movement time - Time from MI to the lever reaching the threshold
    neuralDynamics.movementTime = cell2mat(arrayfun(@(s) s.rewardIndex - s.LFPIndex(parameters.Fs*parameters.windowBeforePull), IntanBehaviour.MIHitTrace, 'UniformOutput', false));
    
    if plotFlag == 1
        % Fitting peak speed between cue and MI to reaction time 
        mdl = fitlm(neuralDynamics.hit.speed.cueMIPeakSpeed,IntanBehaviour.reactionTime)
        figure,plot(mdl);
        xlabel('Peak neural trajectory speed post cue');ylabel('Reaction Time');
        
        % Fitting peak speed between cue and MI to movement time 
        mdl = fitlm(neuralDynamics.MIhit.speed.cueMIPeakSpeed,neuralDynamics.movementTime)
        figure,plot(mdl);
        xlabel('Peak neural trajectory speed post cue');ylabel('Movement time');
        
        % Fitting time location of peak speed between cue and MI to reaction time 
        mdl = fitlm(neuralDynamics.hit.speed.cueMIPeakSpeedTime,IntanBehaviour.reactionTime)
        figure,plot(mdl);
        xlabel('Peak neural trajectory speed post cue');ylabel('Reaction Time');
    end
    
    %% Saving the IntanBehaviour and neuralDyanmics variable 
    M1neuralDynamics(filenumber).IntanBehaviour = IntanBehaviour;
    M1neuralDynamics(filenumber).neuralDynamics = neuralDynamics;
end

%% Pooling data from all animals and sessions 

% SHOULD I NORMALIZE FOR EACH MOUSE?
PooledData.hit.speed.cueMIPeakSpeedTime = [];
PooledData.hit.speed.cueMIPeakSpeed = [];
PooledData.RT = [];
for i = 1:size(M1neuralDynamics,2)
    if M1neuralDynamics(i).skip == 1
        continue
    else
        PooledData.RT = [PooledData.RT ,  M1neuralDynamics(i).IntanBehaviour.reactionTime];
        PooledData.hit.speed.cueMIPeakSpeedTime = ([PooledData.hit.speed.cueMIPeakSpeedTime , M1neuralDynamics(i).neuralDynamics.hit.speed.cueMIPeakSpeedTime ]);
        PooledData.hit.speed.cueMIPeakSpeed = ([PooledData.hit.speed.cueMIPeakSpeed , M1neuralDynamics(i).neuralDynamics.hit.speed.cueMIPeakSpeed ]);
    end
end
% Fitting time location of peak speed between cue and MI to reaction time 
neuralTrajectoryTs = 20; % in ms
mdl = fitlm(neuralTrajectoryTs*PooledData.hit.speed.cueMIPeakSpeedTime,PooledData.RT)
figure,plot(mdl);
xlabel('Time of peak neural trajectory speed from cue');ylabel('Reaction Time');

% Fitting time location of peak speed between cue and MI to reaction time 
mdl = fitlm(PooledData.hit.speed.cueMIPeakSpeed,PooledData.RT)
figure,plot(mdl);
xlabel('Peak neural trajectory speed from cue');ylabel('Reaction Time');

[ ~,outlierIndex] = rmoutliers(PooledData.RT);
PooledData.hit.speed.cueMIPeakSpeedRMOutlier = neuralTrajectoryTs*PooledData.hit.speed.cueMIPeakSpeed;
PooledData.hit.speed.cueMIPeakSpeedRMOutlier(outlierIndex) = [];
PooledData.RTRMOutlier = PooledData.RT; PooledData.RTRMOutlier(outlierIndex)=[];
mdl = fitlm(PooledData.hit.speed.cueMIPeakSpeedRMOutlier,PooledData.RTRMOutlier)
figure,plot(mdl);
xlabel('Peak neural trajectory speed from cue');ylabel('Reaction Time');




