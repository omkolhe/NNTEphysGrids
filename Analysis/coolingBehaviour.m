function coolingBehaviour(IntanBehaviour,parameters)

%% Comparing Reaction Time
baselineRT = [];
coolRT = [];
recoveryRT = [];
recoveryFlag = 0; % Run before loop; recoveryFlag = 0 - Baseline, = 1 - Recovery
% only works for 1 baseline, cool and recovery cycle
for i=1:size(IntanBehaviour.cueHitTrace,2)  
    if IntanBehaviour.cueHitTrace(i).temp <= 20 
        coolRT = [coolRT;IntanBehaviour.cueHitTrace(i).reactionTime];
        recoveryFlag = 1; 
    elseif IntanBehaviour.cueHitTrace(i).temp >= 25 && recoveryFlag == 0
        baselineRT= [baselineRT;IntanBehaviour.cueHitTrace(i).reactionTime];
    elseif IntanBehaviour.cueHitTrace(i).temp >= 25 && recoveryFlag == 1
        recoveryRT= [recoveryRT;IntanBehaviour.cueHitTrace(i).reactionTime];
    end
end

ranksum(coolRT,baselineRT)
ranksum(coolRT,recoveryRT)
ranksum(baselineRT,recoveryRT)

data = baselineRT;
% manually copy noOptoRT into data

figure,customBoxplot(data);

%% Analog control of RT wrt to temp

temp = cell2mat(arrayfun(@(s) s.temp, IntanBehaviour.cueHitTrace, 'UniformOutput', false));
RT = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false));

mdl = fitlm(temp,RT)
figure,plot(mdl);
