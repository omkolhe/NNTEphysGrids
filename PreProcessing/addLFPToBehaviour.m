function IntanBehaviour = addLFPToBehaviour(IntanBehaviour,LFP,parameters)

if strcmp(parameters.experiment,'cue')
    % Adding LFP an GP traces to cueHitTrace variable 
    for i=1:size(IntanBehaviour.cueHitTrace,2)
        IntanBehaviour.cueHitTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        IntanBehaviour.cueHitTrace(i).xf = LFP.xf(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        IntanBehaviour.cueHitTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        IntanBehaviour.cueHitTrace(i).wt = LFP.wt(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
    end
    
    if isfield(IntanBehaviour,'cueMissTrace')
        % Adding LFP an GP traces to cueMissTrace variable 
        for i=1:size(IntanBehaviour.cueMissTrace,2)
            IntanBehaviour.cueMissTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            IntanBehaviour.cueMissTrace(i).xf = LFP.xf(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            IntanBehaviour.cueMissTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            IntanBehaviour.cueMissTrace(i).wt = LFP.wt(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
        end
    end
end

if isfield(IntanBehaviour,'missTrace')
    % Adding LFP an GP traces to missTrace variable 
    for i=1:size(IntanBehaviour.missTrace,2)
        IntanBehaviour.missTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
        IntanBehaviour.missTrace(i).xf = LFP.xf(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
        IntanBehaviour.missTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
        IntanBehaviour.missTrace(i).wt = LFP.wt(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
    end
end

if isfield(IntanBehaviour,'hitTrace')
    % Adding LFP an GP traces to hitTrace variable 
    for i=1:size(IntanBehaviour. hitTrace,2)
        IntanBehaviour. hitTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        IntanBehaviour. hitTrace(i).xf = LFP.xf(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        IntanBehaviour. hitTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        IntanBehaviour. hitTrace(i).wt = LFP.wt(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
    end
end

if isfield(IntanBehaviour,'MIHitTrace')
    % Adding LFP an GP traces to hitTrace variable 
    for i=1:size(IntanBehaviour.  MIHitTrace,2)
        IntanBehaviour.  MIHitTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        IntanBehaviour.  MIHitTrace(i).xf = LFP.xf(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        IntanBehaviour.  MIHitTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        IntanBehaviour.  MIHitTrace(i).wt = LFP.wt(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
    end
end

if isfield(IntanBehaviour,'MIFATrace')
    % Adding LFP an GP traces to  MIFATrace variable 
    for i=1:size(IntanBehaviour.  MIFATrace,2)
        IntanBehaviour.  MIFATrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        IntanBehaviour.  MIFATrace(i).xf = LFP.xf(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        IntanBehaviour.  MIFATrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        IntanBehaviour.  MIFATrace(i).wt = LFP.wt(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
    end
end