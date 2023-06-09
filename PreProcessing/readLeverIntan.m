function [IntanBehaviour] = readLeverIntan(parameters,lfpTime,leverTrace,rewardTrace,Behaviour)
intanFs = 5000;

resting_position = 242*5/1024;
flip = 1;
nlengthBeforePull = round(parameters.windowAfterPull/parameters.ts);
nlength = round(parameters.windowAfterPull/parameters.ts + parameters.windowAfterPull/parameters.ts + 1);

IntanBehaviour.leverTrace = resample(((double(leverTrace)- resting_position)*flip),parameters.Fs,intanFs);
IntanBehaviour.time = lfpTime; % time in seconds
IntanBehaviour.rewardTrace = downsample(rewardTrace,round(intanFs/parameters.Fs),2); 

rewardIndex = find(diff(IntanBehaviour.rewardTrace)==1)+1;

IntanBehaviour.nHit = size(rewardIndex,2);
% IntanBehaviour.nMiss = B(end,4);

% Estimating the threshold for reward
IntanBehaviour.threshold = mean(IntanBehaviour.leverTrace(rewardIndex),'all');

%% Getting hit traces and timings

for i=1:IntanBehaviour.nHit
%     IntanBehaviour.hit(i) = [rewardIndex(i) lfpTime(rewardIndex(i)) rewardIndex(i) lfpTime(rewardIndex(i))];
    IntanBehaviour.hitTrace(i).trace = IntanBehaviour.leverTrace(rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs)';
    IntanBehaviour.hitTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.hitTrace(i).trace,1)-1)*1/parameters.Fs)';
    IntanBehaviour.hitTrace(i).LFPIndex = ([rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:1:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs])';
    IntanBehaviour.hitTrace(i).LFPtime = IntanBehaviour.time(rewardIndex(i)-parameters.windowBeforePull*parameters.Fs:rewardIndex(i)+parameters.windowAfterPull*parameters.Fs)';
end


%% Getting miss traces and timings

% figure();plot(1:1:10001,(5/1024)*Behaviour.leverTrace(Behaviour.miss(15,1)-5000:Behaviour.miss(15,1)+5000));hold on;
% plot(0.1:0.1:10000.1,IntanBehaviour.leverTrace(Behaviour.miss(15,3)-50000:Behaviour.miss(15,3)+50000));
correctionWindow = 600; % in number of points in LFPFs
tol = 0.004;
nDiffSlope = 10;
disp('Finding miss trials in the Intan data ...');
for i=1:Behaviour.nMiss
    missIndexAr = Behaviour.miss(i,3);
    trace1 = IntanBehaviour.leverTrace(missIndexAr-correctionWindow:missIndexAr+correctionWindow);
    misstrigs1 = find(trace1 <IntanBehaviour.threshold+tol & trace1>IntanBehaviour.threshold-tol);  
    % Checking slope 
    for j=1:size(misstrigs1,2)
        % Checking edge cases, rejects all the edge cases 
        if((misstrigs1(j)+nDiffSlope >= (correctionWindow*2+1)) || (misstrigs1(j)-nDiffSlope <= 0)) 
            misstrigs1(j) = NaN;
            continue
        end
        % Checking slope, reject all negative slope 
        slope = mean(trace1(misstrigs1(j):misstrigs1(j)+nDiffSlope)) - mean(trace1(misstrigs1(j)-nDiffSlope:misstrigs1(j)));
        if slope < 0
            misstrigs1(j) = NaN;
        end
    end
    misstrigs1 = misstrigs1(~isnan(misstrigs1));
    missIndex(i) =  missIndexAr - correctionWindow + min(misstrigs1);
end

IntanBehaviour.nMiss = Behaviour.nMiss;

for i=1:IntanBehaviour.nMiss
%     IntanBehaviour.hit(i) = [rewardIndex(i) lfpTime(rewardIndex(i)) rewardIndex(i) lfpTime(rewardIndex(i))];
    IntanBehaviour.missTrace(i).trace = IntanBehaviour.leverTrace(missIndex(i)-parameters.windowBeforePull*parameters.Fs:missIndex(i)+parameters.windowAfterPull*parameters.Fs)';
    IntanBehaviour.missTrace(i).time = (0:1/parameters.Fs:(size(IntanBehaviour.hitTrace(i).trace,1)-1)*1/parameters.Fs)';
    IntanBehaviour.missTrace(i).LFPIndex = ([missIndex(i)-parameters.windowBeforePull*parameters.Fs:1:missIndex(i)+parameters.windowAfterPull*parameters.Fs])';
    IntanBehaviour.missTrace(i).LFPtime = IntanBehaviour.time(missIndex(i)-parameters.windowBeforePull*parameters.Fs:missIndex(i)+parameters.windowAfterPull*parameters.Fs)';
end


