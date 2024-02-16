function [Encoder] = readPos(parameters)
tic

reward = parameters.rewardTrials;
timecorrection = 50e-3; 

[enfile,enpath] = uigetfile('*.csv');
if isequal(enfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(enpath,enfile)]);
end

A = readmatrix([enpath,'/',enfile]);

if reward == 0
    pos = A(:,1);
    time = A(:,2)/1e3;       % time from arduino converted to seconds
    fsAr = 1/(time(2)-time(1));
    vel = zeros(1,numel(time));
    vel(2:end) = movmean(diff(pos)*fsAr,fsAr);
else
    pos = A(:,2);
    time1 = A(:,3)/1e3;       % time from arduino converted to seconds
    fsAr = 1/(time1(2)-time1(1));
    rewardIndX = find(time1 == max(time1));
    time = time1;
    for i=1:numel(rewardIndX)
        time(rewardIndX(i)) = time(rewardIndX(i)-1) + timecorrection;
        time(rewardIndX(i)+1:end) = time1(rewardIndX(i)+1:end) + time(rewardIndX(i)); 
    end
    Encoder.rewardIndx = rewardIndX-1;
    vel = zeros(1,numel(time));
    vel(2:end) = movmean(diff(pos)*fsAr,fsAr);
end
Encoder.pos = pos;
Encoder.vel = vel;
Encoder.time = time;
Encoder.fs = fsAr;

toc
end


%% Previous code

% vel = zeros(1,size(tLFP,2));
% time = tLFP/1000;
% preInd = 1;
% for i=1:size(posAr,1)
%     ind = max(find(tLFP<=tAr(i)));
%     if(i==1)
%         pos(preInd:ind) = 0;
%     else
%         pos(preInd:ind) = posAr(i);
%     end
%     preInd = ind+1;
%     if mod(i,1000)==0 i, end
% end