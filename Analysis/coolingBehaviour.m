function optoBehaviour(IntanBehaviour,parameters)

%% Comparing Reaction Time 
coolRT = [];
noCoolRT = [];

for i=1:size(IntanBehaviour.cueHitTrace,2)  
    if IntanBehaviour.cueHitTrace(i).temp <= 20 
        coolRT = [coolRT;IntanBehaviour.cueHitTrace(i).reactionTime];
    elseif IntanBehaviour.cueHitTrace(i).temp >= 25
        noCoolRT= [noCoolRT;IntanBehaviour.cueHitTrace(i).reactionTime];
    end
end

ranksum(coolRT,noCoolRT)

data = coolRT;
% manually copy noOptoRT into data

figure,customBoxplot(data);

%% Analog control of RT wrt to temp

temp = cell2mat(arrayfun(@(s) s.temp, IntanBehaviour.cueHitTrace, 'UniformOutput', false));
RT = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false));

mdl = fitlm(temp,RT)
figure,plot(mdl);
