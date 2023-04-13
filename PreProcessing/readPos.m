function [pos,vel,time,fsAr] = readPos()
tic

[enfile,enpath] = uigetfile('*.csv');
if isequal(enfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(enpath,enfile)]);
end

A = readmatrix([enpath,'/',enfile]);
pos = A(:,2);
time = A(:,3)/1e3;       % time from arduino converted to seconds
fsAr = 1/(time(2)-time(1));

vel = zeros(1,numel(time));

vel(2:end) = movmean(diff(pos)*fsAr,fsAr);
toc


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