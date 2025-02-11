function [neuralDynamics] = matchBehaviourNeuralDynamics(neuralDynamics,IntanBehaviour)
% Rejecting neural trajectories according to IntanBehaviour 
% Hit
neuralDynamics.hit.X(:,:,IntanBehaviour.rejectTrials.Hit) = [];
% Miss
neuralDynamics.miss.X(:,:,IntanBehaviour.rejectTrials.cueMiss) = [];
% MIHit
neuralDynamics.MIhit.X(:,:,IntanBehaviour.rejectTrials.Hit) = [];
% MIFA
neuralDynamics.MIFA.X(:,:,IntanBehaviour.rejectTrials.FA) = [];

neuralDynamics.matchFlag = 1;
end
