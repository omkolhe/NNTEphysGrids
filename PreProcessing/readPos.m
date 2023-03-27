function [pos,vel,time] = readPos(tLFP)
tic
tLFP = tLFP*1000;

[enfile,enpath] = uigetfile('*.csv');
if isequal(enfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(enpath,enfile)]);
end

A = readmatrix([enpath,'/',enfile]);
cngPos = A(:,1);    % change in pos from arduino
posAr = cumsum(cngPos);
tAr = A(:,2);       % time from arduino in ms

pos = zeros(1,size(tLFP,2));
vel = zeros(1,size(tLFP,2));
time = tLFP/1000;
preInd = 1;
for i=1:size(posAr,1)
    ind = max(find(tLFP<=tAr(i)));
    if(i==1)
        pos(preInd:ind) = 0;
    else
        pos(preInd:ind) = posAr(i);
    end
    preInd = ind+1;
    if mod(i,1000)==0 i, end
end
vel(2:end) = movmean(diff(pos)*1000/(tLFP(2)-tLFP(1)),500);
toc
