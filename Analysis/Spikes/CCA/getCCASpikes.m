function spikes = getCCASpikes(M1Spikes,M2Spikes,behaviourTrace,delay)

%% Positive delay -> M2Spikes lag behind M1Spikes

nTrials = size(behaviourTrace,2);

spikes{1,1} = zeros(M1Spikes.nSpikes,size(behaviourTrace(1).trace,1),nTrials);
spikes{1,2} = zeros(M2Spikes.nSpikes,size(behaviourTrace(1).trace,1),nTrials);

for i=1:nTrials
    spikes{1,1}(:,:,i) = M1Spikes.spikes(:,behaviourTrace(i).LFPIndex(1):behaviourTrace(i).LFPIndex(end));
    spikes{1,2}(:,:,i) = M2Spikes.spikes(:,behaviourTrace(i).LFPIndex(1)-delay:behaviourTrace(i).LFPIndex(end)-delay);
end

end
% 
% %% Creating spike strucutre for CCA Analysis for Region 1 (M1)
% % Cue Hits 
% M1Spikes.CCASpikes.cueHit = zeros(size(M1Spikes.PSTH.hit.spks',1),size(M1Spikes.PSTH.hit.spks{1,1},2),size(M1Spikes.PSTH.hit.spks{1,1},1));
% for i=1:size(M1Spikes.CCASpikes.cueHit,1) % Number of neurons
%     for j=1:size(M1Spikes.CCASpikes.cueHit,3) % Number of trial
%         M1Spikes.CCASpikes.cueHit(i,:,j) = M1Spikes.PSTH.hit.spks{1,i}(j,:);
%     end
% end
% % Cue Miss
% M1Spikes.CCASpikes.cueMiss = zeros(size(M1Spikes.PSTH.miss.spks',1),size(M1Spikes.PSTH.miss.spks{1,1},2),size(M1Spikes.PSTH.miss.spks{1,1},1));
% for i=1:size(M1Spikes.CCASpikes.cueMiss,1) % Number of neurons
%     for j=1:size(M1Spikes.CCASpikes.cueMiss,3) % Number of trial
%         M1Spikes.CCASpikes.cueMiss(i,:,j) = M1Spikes.PSTH.miss.spks{1,i}(j,:);
%     end
% end
% % MI Hits
% M1Spikes.CCASpikes.MIHit = zeros(size(M1Spikes.PSTH.MIHit.spks',1),size(M1Spikes.PSTH.MIHit.spks{1,1},2),size(M1Spikes.PSTH.MIHit.spks{1,1},1));
% for i=1:size(M1Spikes.CCASpikes.MIHit,1) % Number of neurons
%     for j=1:size(M1Spikes.CCASpikes.MIHit,3) % Number of trial
%         M1Spikes.CCASpikes.MIHit(i,:,j) = M1Spikes.PSTH.MIHit.spks{1,i}(j,:);
%     end
% end
% % MI FA
% M1Spikes.CCASpikes.MIFA = zeros(size(M1Spikes.PSTH.MIFA.spks',1),size(M1Spikes.PSTH.MIFA.spks{1,1},2),size(M1Spikes.PSTH.MIFA.spks{1,1},1));
% for i=1:size(M1Spikes.CCASpikes.MIFA,1) % Number of neurons
%     for j=1:size(M1Spikes.CCASpikes.MIFA,3) % Number of trial
%         M1Spikes.CCASpikes.MIFA(i,:,j) = M1Spikes.PSTH.MIFA.spks{1,i}(j,:);
%     end
% end
% 
% %% Creating spike strucutre for CCA Analysis for Region 2 (M2)
% % Cue Hits 
% M2Spikes.CCASpikes.cueHit = zeros(size(M2Spikes.PSTH.hit.spks',1),size(M2Spikes.PSTH.hit.spks{1,1},2),size(M2Spikes.PSTH.hit.spks{1,1},1));
% for i=1:size(M2Spikes.CCASpikes.cueHit,1) % Number of neurons
%     for j=1:size(M2Spikes.CCASpikes.cueHit,3) % Number of trial
%         M2Spikes.CCASpikes.cueHit(i,:,j) = M2Spikes.PSTH.hit.spks{1,i}(j,:);
%     end
% end
% % Cue Miss
% M2Spikes.CCASpikes.cueMiss = zeros(size(M2Spikes.PSTH.miss.spks',1),size(M2Spikes.PSTH.miss.spks{1,1},2),size(M2Spikes.PSTH.miss.spks{1,1},1));
% for i=1:size(M2Spikes.CCASpikes.cueMiss,1) % Number of neurons
%     for j=1:size(M2Spikes.CCASpikes.cueMiss,3) % Number of trial
%         M2Spikes.CCASpikes.cueMiss(i,:,j) = M2Spikes.PSTH.miss.spks{1,i}(j,:);
%     end
% end
% % MI Hits
% M2Spikes.CCASpikes.MIHit = zeros(size(M2Spikes.PSTH.MIHit.spks',1),size(M2Spikes.PSTH.MIHit.spks{1,1},2),size(M2Spikes.PSTH.MIHit.spks{1,1},1));
% for i=1:size(M2Spikes.CCASpikes.MIHit,1) % Number of neurons
%     for j=1:size(M2Spikes.CCASpikes.MIHit,3) % Number of trial
%         M2Spikes.CCASpikes.MIHit(i,:,j) = M2Spikes.PSTH.MIHit.spks{1,i}(j,:);
%     end
% end
% % MI FA
% M2Spikes.CCASpikes.MIFA = zeros(size(M2Spikes.PSTH.MIFA.spks',1),size(M2Spikes.PSTH.MIFA.spks{1,1},2),size(M2Spikes.PSTH.MIFA.spks{1,1},1));
% for i=1:size(M2Spikes.CCASpikes.MIFA,1) % Number of neurons
%     for j=1:size(M2Spikes.CCASpikes.MIFA,3) % Number of trial
%         M2Spikes.CCASpikes.MIFA(i,:,j) = M2Spikes.PSTH.MIFA.spks{1,i}(j,:);
%     end
% end
