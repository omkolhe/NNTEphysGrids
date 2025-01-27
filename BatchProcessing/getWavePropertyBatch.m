function [wavesHit,wavesMiss,wavesFA,wavesHitReward,wavesMIHit,wavesMIFA] = getWavePropertyBatch(filename,filepath,fileID,waveProperties,wavesHit,wavesMiss,wavesFA,wavesHitReward,wavesMIHit,wavesMIFA )

%% Checking if variable exists
disp(['Reading from '+filepath+'\'+filename]);
matObj = matfile(filepath+'\'+filename);
variables = who(matObj);
skip = 0;
skipMissReject = 0;
if ismember('Waves', variables)
    disp("Loading variables - paramters and Waves");
    parameters = load(filepath+'\'+filename,"parameters");
    parameters = parameters.parameters;
    Waves = load(filepath+'\'+filename,"Waves");
    Waves = Waves.Waves;
elseif ismember('WavesBaseline',variables)
    disp("Loading variables - paramters and WavesBaseline");
    parameters = load(filepath+'\'+filename,"parameters");
    parameters = parameters.parameters;
    Waves = load(filepath+'\'+filename,"WavesBaseline");
    Waves = Waves.WavesBaseline;
else
    disp("Could not find Waves or WavesBaseline");
    skip = 1;
end

if ismember('IntanBehaviour', variables)
    disp("Loading variables - IntanBehaviour");
    IntanBehaviour = load(filepath+'\'+filename,"IntanBehaviour");
    IntanBehaviour = IntanBehaviour.IntanBehaviour;
elseif ismember('IntanBehaviourBaseline',variables)
    disp("Loading variables - IntanBehaviour");
    IntanBehaviour = load(filepath+'\'+filename,"IntanBehaviourBaseline");
    IntanBehaviour = IntanBehaviour.IntanBehaviourBaseline;
else
    disp("Could not find IntanBehaviour or IntanBehaviourBaseline");
    skipMissReject = 1;
end

%% Rejecting all misses after the last hit + 40 misses
if skipMissReject == 0
    disp("Rejecting extra miss trials at the end");
    extraMisses = 40;
    lastHitIndex = IntanBehaviour.cueHitTrace(end).LFPIndex(end);
    lastMiss = 1;
    for i=1:size(IntanBehaviour.cueMissTrace,2)
        if IntanBehaviour.cueMissTrace(i).LFPIndex(end) > lastHitIndex
            lastMiss = i;
            break;
        end
    end
    
    if lastMiss+extraMisses > size(IntanBehaviour.cueMissTrace,2)
        lastMiss = size(IntanBehaviour.cueMissTrace,2);
    else
        lastMiss = lastMiss+extraMisses;
    end
    
    disp(['Rejected last ', string(size(IntanBehaviour.cueMissTrace,2)-lastMiss)]);
end

%% Getting wave properties
if skip == 0
    % Hits 
    nWavesCombinedHits = size( wavesHit,2);
    for i=1:size(Waves.wavesHit,2)
        for j=1:length(waveProperties)
            [ wavesHit(nWavesCombinedHits+i).(waveProperties{j})] = Waves.wavesHit(i).(waveProperties{j});
        end
         wavesHit(nWavesCombinedHits+i).fileID = fileID;
         wavesHit(nWavesCombinedHits+i).RT = IntanBehaviour.cueHitTrace(i).reactionTime;
    end
    
    % Miss 
    nWavesCombinedMiss = size( wavesMiss,2);
    for i=1:lastMiss
        for j=1:length(waveProperties)
            [ wavesMiss(nWavesCombinedMiss+i).(waveProperties{j})] = Waves.wavesMiss(i).(waveProperties{j});
        end
         wavesMiss(nWavesCombinedMiss+i).fileID = fileID;
    end
    
    % FA 
    nWavesCombinedFA = size( wavesFA,2);
    for i=1:size(Waves.wavesFA,2)
        for j=1:length(waveProperties)
            [ wavesFA(nWavesCombinedFA+i).(waveProperties{j})] = Waves.wavesFA(i).(waveProperties{j});
        end
         wavesFA(nWavesCombinedFA+i).fileID = fileID;
    end
    
    % Hits Reward aligned
    nWavesCombinedHitReward = size( wavesHitReward,2);
    for i=1:size(Waves.wavesHitReward,2)
        for j=1:length(waveProperties)
            [ wavesHitReward(nWavesCombinedHitReward+i).(waveProperties{j})] = Waves.wavesHitReward(i).(waveProperties{j});
        end
         wavesHitReward(nWavesCombinedHitReward+i).fileID = fileID;
    end
    
    % Hits MI
    nWavesCombinedMIHit = size( wavesMIHit,2);
    for i=1:size(Waves.wavesMIHit,2)
        for j=1:length(waveProperties)
            [ wavesMIHit(nWavesCombinedMIHit+i).(waveProperties{j})] = Waves.wavesMIHit(i).(waveProperties{j});
        end
         wavesMIHit(nWavesCombinedMIHit+i).fileID = fileID;
    end
    
    % Hits MI
    nWavesCombinedMIFA = size( wavesMIFA,2);
    for i=1:size(Waves.wavesMIFA,2)
        for j=1:length(waveProperties)
            [ wavesMIFA(nWavesCombinedMIFA+i).(waveProperties{j})] = Waves.wavesMIFA(i).(waveProperties{j});
        end
         wavesMIFA(nWavesCombinedMIFA+i).fileID = fileID;
    end

    % Clearing loaded variables 
    clear Waves IntanBehaviour;
end

end

