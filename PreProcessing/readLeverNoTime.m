function [Behaviour] = readLeverNoTime(parameters)

% if ~exist('parameters.experiment','var')
%     parameters.experiment = 'self';
%     disp('No experiment argument passed. Experiment type set to self initiated');
% end

if strcmp(parameters.experiment,'cue')
    cue = 1;
    disp('Experiment type set to cue initiated. . . ')
else
    cue = 0;
    disp('Experiment type set to self initiated. . . ')
end

%% Reading file from arduino 
[enfile,enpath] = uigetfile('*.csv');
if isequal(enfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(enpath,enfile)]);
end

resting_position = 241;
flip = 1;
nlengthBeforePull = round(parameters.windowBeforePull/parameters.ts);
nlength = round(parameters.windowBeforePull/parameters.ts + parameters.windowAfterPull/parameters.ts + 1);
nlengthBeforeCue = round(parameters.windowBeforeCue/parameters.ts);
nlengthCue = round(parameters.windowBeforeCue/parameters.ts + parameters.windowAfterCue/parameters.ts + 1);

B = readmatrix([enpath,'/',enfile]);
Behaviour.leverTrace = (B(2:end,1) - resting_position)*flip;
Behaviour.time = (B(2:end,2) - B(2,2))/1e6; % time in seconds
Behaviour.nHit = B(end,3);
Behaviour.nMiss = B(end,4);
if cue==1 
    Behaviour.nCue = B(end,5); 
    Behaviour.nCueHit = Behaviour.nHit;
    Behaviour.nCueMiss = Behaviour.nCue-Behaviour.nCueHit;
end   
Behaviour.B = B(2:end,:);

%% Getting hit and miss timings
hitIndex = find(diff(B(:,3)) == 1) + 1;
hitTime = Behaviour.time(hitIndex);
Behaviour.hit = [hitIndex hitTime];

missIndex = find(diff(B(:,4)) == 1) + 1;
missTime = Behaviour.time(missIndex);
Behaviour.miss = [missIndex missTime ];

%% Getting cues, hit cues and miss cues
if cue == 1
    cueIndex = find(Behaviour.B(:,end) == 1);
    cueTime = Behaviour.time(cueIndex);
    Behaviour.cue = [cueIndex cueTime];
    
    % Getting cueHits and cueMisses 
    cueHitIndex = zeros(Behaviour.nHit,1);
    cueHitTime = zeros(Behaviour.nHit,1);
    cueHitPullIndex = zeros(Behaviour.nHit,1);
    cueHitPullTime = zeros(Behaviour.nHit,1);
    
    a = zeros(Behaviour.nHit,1);
    
    for i=1:Behaviour.nHit
        a(i) = max(find(cueTime<hitTime(i))); % No need to check for reaction time. If there is a hit, there is a cue
        cueHitIndex(i) = cueIndex(a(i));
        cueHitTime(i) = cueTime(a(i));
        cueHitPullIndex(i) = hitIndex(i);
        cueHitPullTime(i) = hitTime(i);
    end
    
    Behaviour.cueHit = [cueHitIndex cueHitTime cueHitPullIndex cueHitPullTime];
    
    Behaviour.reactionTime = cueHitPullTime - cueHitTime;
    Behaviour.meanReactionTime = mean(Behaviour.reactionTime,'all');
    
    % For cue Miss trials
    Behaviour.cueMiss = Behaviour.cue;
    Behaviour.cueMiss(a,:) = [];
end

%% get lever traces for hits and miss 
st_hit1 = max(find(Behaviour.time < Behaviour.hit(1,2)-parameters.windowBeforePull));
if isempty(st_hit1)
    disp('First hit rejected');
    Behaviour.nHit = Behaviour.nHit-1;
    Behaviour.hit(1,:) = [];
end 
sp_hitend =  max(find(Behaviour.time < Behaviour.hit(end,2)+parameters.windowAfterPull));
if isempty(sp_hitend) 
    disp('Last hit rejected')
    Behaviour.nHit = Behaviour.nHit-1;
    Behaviour.hit(end,:) = [];
end

for i=1:Behaviour.nHit
    Behaviour.hitTrace(i).i1 = max(find(Behaviour.time < Behaviour.hit(i,2)-parameters.windowBeforePull));
    Behaviour.hitTrace(i).i0 = Behaviour.hit(i,1);
    Behaviour.hitTrace(i).i2 = max(find(Behaviour.time < Behaviour.hit(i,2)+parameters.windowAfterPull));
    Behaviour.hitTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2);
    Behaviour.hitTrace(i).rawtime = Behaviour.time(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2) - Behaviour.time(Behaviour.hitTrace(i).i1);
    Behaviour.hitTrace(i).time1 = Behaviour.time(Behaviour.hitTrace(i).i1:Behaviour.hitTrace(i).i2);
    Behaviour.hitTrace(i).t1 = Behaviour.time(Behaviour.hitTrace(i).i1);
    Behaviour.hitTrace(i).t0 = Behaviour.hit(i,2);
    Behaviour.hitTrace(i).t2 = Behaviour.time(Behaviour.hitTrace(i).i2);
end

st_miss1 = max(find(Behaviour.time < Behaviour.miss(1,2)-parameters.windowBeforePull));
if isempty(st_miss1)
    disp('First miss rejected');
    Behaviour.nMiss = Behaviour.nMiss-1;
    Behaviour.miss(1,:) = [];
end 
sp_missend =  max(find(Behaviour.time < Behaviour.miss(end,2)+parameters.windowAfterPull));
if isempty(sp_missend)
    disp('Last miss rejected')
    Behaviour.nMiss = Behaviour.nMiss-1;
    Behaviour.miss(end,:) = [];
end

for i=1:Behaviour.nMiss
    Behaviour.missTrace(i).i1 = max(find(Behaviour.time < Behaviour.miss(i,2)-parameters.windowBeforePull));
    Behaviour.missTrace(i).i0 = Behaviour.miss(i,1);
    Behaviour.missTrace(i).i2 = max(find(Behaviour.time < Behaviour.miss(i,2)+parameters.windowAfterPull));
    Behaviour.missTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2);
    Behaviour.missTrace(i).rawtime = Behaviour.time(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2) - Behaviour.time(Behaviour.missTrace(i).i1);
    Behaviour.missTrace(i).time1 = Behaviour.time(Behaviour.missTrace(i).i1:Behaviour.missTrace(i).i2);
    Behaviour.missTrace(i).t1 = Behaviour.time(Behaviour.missTrace(i).i1);
    Behaviour.missTrace(i).t0 = Behaviour.miss(i,2);
    Behaviour.missTrace(i).t2 = Behaviour.time(Behaviour.missTrace(i).i2);
end

if cue == 1 
    st_cuehit1 = max(find(Behaviour.time < Behaviour.cueHit(1,2)-parameters.windowBeforeCue));
    if isempty(st_cuehit1)
        disp('First cue hit rejected');
        Behaviour.nCueHit = Behaviour.nCueHit-1;
        Behaviour.cueHit(1,:) = [];
    end 
    sp_cuehitend =  max(find(Behaviour.time < Behaviour.cueHit(end,2)+parameters.windowAfterCue));
    if isempty(sp_cuehitend)
        disp('Last cue hit rejected')
        Behaviour.nCueHit = Behaviour.nCueHit-1;
        Behaviour.cueHit(end,:) = [];
    end
    for i=1:Behaviour.nCueHit
        Behaviour.cueHitTrace(i).i1 = max(find(Behaviour.time < Behaviour.cueHit(i,2)-parameters.windowBeforeCue));
        Behaviour.cueHitTrace(i).i0 = Behaviour.cueHit(i,1);
        Behaviour.cueHitTrace(i).i2 = max(find(Behaviour.time < Behaviour.cueHit(i,2)+parameters.windowAfterCue));
        Behaviour.cueHitTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.cueHitTrace(i).i1:Behaviour.cueHitTrace(i).i2);
        Behaviour.cueHitTrace(i).rawtime = Behaviour.time(Behaviour.cueHitTrace(i).i1:Behaviour.cueHitTrace(i).i2) - Behaviour.time(Behaviour.cueHitTrace(i).i1);
        Behaviour.cueHitTrace(i).time1 = Behaviour.time(Behaviour.cueHitTrace(i).i1:Behaviour.cueHitTrace(i).i2);
        Behaviour.cueHitTrace(i).t1 = Behaviour.time(Behaviour.cueHitTrace(i).i1);
        Behaviour.cueHitTrace(i).t0 = Behaviour.cueHit(i,2);
        Behaviour.cueHitTrace(i).t2 = Behaviour.time(Behaviour.cueHitTrace(i).i2);
    end
    
    % Getting cueMiss traces
    if Behaviour.nCueMiss > 0
        st_cuemiss1 = max(find(Behaviour.time < Behaviour.cueMiss(1,2)-parameters.windowBeforeCue));
        if isempty(st_cuemiss1)
            disp('First cue miss rejected');
            Behaviour.nCueMiss = Behaviour.nCueMiss-1;
            Behaviour.cueMiss(1,:) = [];
        end 
        sp_cuemissend =  max(find(Behaviour.time < Behaviour.cueMiss(end,2)+parameters.windowAfterCue));
        if isempty(sp_cuemissend)
            disp('Last cue hit rejected')
            Behaviour.nCueMiss = Behaviour.nCueMiss-1;
            Behaviour.cueMiss(end,:) = [];
        end
        for i=1:Behaviour.nCueMiss
            Behaviour.cueMissTrace(i).i1 = max(find(Behaviour.time < Behaviour.cueMiss(i,2)-parameters.windowBeforeCue));
            Behaviour.cueMissTrace(i).i0 = Behaviour.cueMiss(i,1);
            Behaviour.cueMissTrace(i).i2 = max(find(Behaviour.time < Behaviour.cueMiss(i,2)+parameters.windowAfterCue));
            Behaviour.cueMissTrace(i).rawtrace = Behaviour.leverTrace(Behaviour.cueMissTrace(i).i1:Behaviour.cueMissTrace(i).i2);
            Behaviour.cueMissTrace(i).rawtime = Behaviour.time(Behaviour.cueMissTrace(i).i1:Behaviour.cueMissTrace(i).i2) - Behaviour.time(Behaviour.cueMissTrace(i).i1);
            Behaviour.cueMissTrace(i).time1 = Behaviour.time(Behaviour.cueMissTrace(i).i1:Behaviour.cueMissTrace(i).i2);
            Behaviour.cueMissTrace(i).t1 = Behaviour.time(Behaviour.cueMissTrace(i).i1);
            Behaviour.cueMissTrace(i).t0 = Behaviour.cueMiss(i,2);
            Behaviour.cueMissTrace(i).t2 = Behaviour.time(Behaviour.cueMissTrace(i).i2);
        end
    end
end

