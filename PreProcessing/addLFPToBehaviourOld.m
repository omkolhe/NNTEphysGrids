function IntanBehaviour = addLFPToBehaviour(IntanBehaviour,LFP,parameters)

betaPresent = 0;
gammaPresent = 0;

if strcmp(parameters.experiment,'cue')
    % Adding LFP an GP traces to cueHitTrace variable 
    for i=1:size(IntanBehaviour.cueHitTrace,2)
        IntanBehaviour.cueHitTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        IntanBehaviour.cueHitTrace(i).xf = LFP.xf(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        IntanBehaviour.cueHitTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        IntanBehaviour.cueHitTrace(i).wt = LFP.wt(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        if parameters.shank == 1
            IntanBehaviour.cueHitTrace(i).rawLFPProbe = LFP.LFPprobe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
            IntanBehaviour.cueHitTrace(i).xfProbe = LFP.xfProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
            IntanBehaviour.cueHitTrace(i).xgpProbe = LFP.xgpProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
            IntanBehaviour.cueHitTrace(i).wtProbe = LFP.wtProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        end
        if betaPresent == 1
            IntanBehaviour.cueHitTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
            IntanBehaviour.cueHitTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
            IntanBehaviour.cueHitTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        end
%         IntanBehaviour.cueHitTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%         IntanBehaviour.cueHitTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        if gammaPresent == 1
            IntanBehaviour.cueHitTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
            IntanBehaviour.cueHitTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
            IntanBehaviour.cueHitTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        end
        if parameters.shank == 1
%             IntanBehaviour.cueHitTrace(i).xfbetaProbe = LFP.xfbetaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%             IntanBehaviour.cueHitTrace(i).xgpbetaProbe = LFP.xgpbetaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%             IntanBehaviour.cueHitTrace(i).wtbetaProbe = LFP.wtbetaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%             IntanBehaviour.cueHitTrace(i).xfthetaProbe = LFP.xfthetaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%             IntanBehaviour.cueHitTrace(i).xgpthetaProbe = LFP.xgpthetaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%             IntanBehaviour.cueHitTrace(i).wtthetaProbe = LFP.wtthetaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%             IntanBehaviour.cueHitTrace(i).xfgammaProbe = LFP.xfgammaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%             IntanBehaviour.cueHitTrace(i).xgpgammaProbe = LFP.xgpgammaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
%             IntanBehaviour.cueHitTrace(i).wtgammaProbe = LFP.wtgammaProbe(:,:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
        end
    end

    if isfield(IntanBehaviour,'cueMissTrace')
        % Adding LFP an GP traces to cueMissTrace variable
        for i=1:size(IntanBehaviour.cueMissTrace,2)
            IntanBehaviour.cueMissTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            IntanBehaviour.cueMissTrace(i).xf = LFP.xf(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            IntanBehaviour.cueMissTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            IntanBehaviour.cueMissTrace(i).wt = LFP.wt(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            if parameters.shank == 1
                IntanBehaviour.cueMissTrace(i).rawLFPProbe = LFP.LFPprobe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
                IntanBehaviour.cueMissTrace(i).xfProbe = LFP.xfProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
                IntanBehaviour.cueMissTrace(i).xgpProbe = LFP.xgpProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
                IntanBehaviour.cueMissTrace(i).wtProbe = LFP.wtProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            end
            if betaPresent == 1
                IntanBehaviour.cueMissTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
                IntanBehaviour.cueMissTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
                IntanBehaviour.cueMissTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            end
%             IntanBehaviour.cueMissTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%             IntanBehaviour.cueMissTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            if gammaPresent == 1
                IntanBehaviour.cueMissTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
                IntanBehaviour.cueMissTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
                IntanBehaviour.cueMissTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            end
            if parameters.shank == 1
%                 IntanBehaviour.cueMissTrace(i).xfbetaProbe = LFP.xfbetaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%                 IntanBehaviour.cueMissTrace(i).xgpbetaProbe = LFP.xgpbetaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%                 IntanBehaviour.cueMissTrace(i).wtbetaProbe = LFP.wtbetaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%                 IntanBehaviour.cueMissTrace(i).xfthetaProbe = LFP.xfthetaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%                 IntanBehaviour.cueMissTrace(i).xgpthetaProbe = LFP.xgpthetaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%                 IntanBehaviour.cueMissTrace(i).wtthetaProbe = LFP.wtthetaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%                 IntanBehaviour.cueMissTrace(i).xfgammaProbe = LFP.xfgammaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%                 IntanBehaviour.cueMissTrace(i).xgpgammaProbe = LFP.xgpgammaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
%                 IntanBehaviour.cueMissTrace(i).wtgammaProbe = LFP.wtgammaProbe(:,:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
            end
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
        if parameters.shank == 1
            IntanBehaviour.missTrace(i).rawLFPProbe = LFP.LFPprobe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
            IntanBehaviour.missTrace(i).xfProbe = LFP.xfProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
            IntanBehaviour.missTrace(i).xgpProbe = LFP.xgpProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
            IntanBehaviour.missTrace(i).wtProbe = LFP.wtProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
        end
        if betaPresent == 1
            IntanBehaviour.missTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
            IntanBehaviour.missTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
            IntanBehaviour.missTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
        end
%         IntanBehaviour.missTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%         IntanBehaviour.missTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
        if gammaPresent == 1
            IntanBehaviour.missTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
            IntanBehaviour.missTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
            IntanBehaviour.missTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
        end
        if parameters.shank == 1
%                 IntanBehaviour.missTrace(i).xfbetaProbe = LFP.xfbetaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%                 IntanBehaviour.missTrace(i).xgpbetaProbe = LFP.xgpbetaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%                 IntanBehaviour.missTrace(i).wtbetaProbe = LFP.wtbetaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%                 IntanBehaviour.missTrace(i).xfthetaProbe = LFP.xfthetaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%                 IntanBehaviour.missTrace(i).xgpthetaProbe = LFP.xgpthetaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%                 IntanBehaviour.missTrace(i).wtthetaProbe = LFP.wtthetaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%                 IntanBehaviour.missTrace(i).xfgammaProbe = LFP.xfgammaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%                 IntanBehaviour.missTrace(i).xgpgammaProbe = LFP.xgpgammaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
%                 IntanBehaviour.missTrace(i).wtgammaProbe = LFP.wtgammaProbe(:,:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
        end
    end
end

if isfield(IntanBehaviour,'hitTrace')
    % Adding LFP an GP traces to hitTrace variable 
    for i=1:size(IntanBehaviour. hitTrace,2)
        IntanBehaviour. hitTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        IntanBehaviour. hitTrace(i).xf = LFP.xf(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        IntanBehaviour. hitTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        IntanBehaviour. hitTrace(i).wt = LFP.wt(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        if parameters.shank == 1
            IntanBehaviour. hitTrace(i).rawLFPProbe = LFP.LFPprobe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
            IntanBehaviour. hitTrace(i).xfProbe = LFP.xfProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
            IntanBehaviour. hitTrace(i).xgpProbe = LFP.xgpProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
            IntanBehaviour. hitTrace(i).wtProbe = LFP.wtProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        end
        if betaPresent == 1
            IntanBehaviour. hitTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
            IntanBehaviour. hitTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
            IntanBehaviour. hitTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        end
%         IntanBehaviour. hitTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
%         IntanBehaviour. hitTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        if gammaPresent == 1
            IntanBehaviour. hitTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
            IntanBehaviour. hitTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
            IntanBehaviour. hitTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour. hitTrace(i).LFPIndex(1):IntanBehaviour. hitTrace(i).LFPIndex(end));
        end
        if parameters.shank == 1
%                 IntanBehaviour.hitTrace(i).xfbetaProbe = LFP.xfbetaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
%                 IntanBehaviour.hitTrace(i).xgpbetaProbe = LFP.xgpbetaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
%                 IntanBehaviour.hitTrace(i).wtbetaProbe = LFP.wtbetaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
%                 IntanBehaviour.hitTrace(i).xfthetaProbe = LFP.xfthetaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
%                 IntanBehaviour.hitTrace(i).xgpthetaProbe = LFP.xgpthetaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
%                 IntanBehaviour.hitTrace(i).wtthetaProbe = LFP.wtthetaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
%                 IntanBehaviour.hitTrace(i).xfgammaProbe = LFP.xfgammaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
%                 IntanBehaviour.hitTrace(i).xgpgammaProbe = LFP.xgpgammaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
%                 IntanBehaviour.hitTrace(i).wtgammaProbe = LFP.wtgammaProbe(:,:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
        end
    end
end

if isfield(IntanBehaviour,'MIHitTrace')
    % Adding LFP an GP traces to hitTrace variable 
    for i=1:size(IntanBehaviour.  MIHitTrace,2)
        IntanBehaviour.  MIHitTrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        IntanBehaviour.  MIHitTrace(i).xf = LFP.xf(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        IntanBehaviour.  MIHitTrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        IntanBehaviour.  MIHitTrace(i).wt = LFP.wt(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        if parameters.shank == 1
            IntanBehaviour.  MIHitTrace(i).rawLFPProbe = LFP.LFPprobe(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
            IntanBehaviour.  MIHitTrace(i).xfProbe = LFP.xfProbe(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
            IntanBehaviour.  MIHitTrace(i).xgpProbe = LFP.xgpProbe(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
            IntanBehaviour.  MIHitTrace(i).wtProbe = LFP.wtProbe(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        end
        if betaPresent == 1
            IntanBehaviour.  MIHitTrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
            IntanBehaviour.  MIHitTrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
            IntanBehaviour.  MIHitTrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        end
%         IntanBehaviour.  MIHitTrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
%         IntanBehaviour.  MIHitTrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        if gammaPresent == 1
            IntanBehaviour.  MIHitTrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
            IntanBehaviour.  MIHitTrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
            IntanBehaviour.  MIHitTrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.  MIHitTrace(i).LFPIndex(1):IntanBehaviour.  MIHitTrace(i).LFPIndex(end));
        end
        if parameters.shank == 1
%             IntanBehaviour.MIHitTrace(i).xfbetaProbe = LFP.xfbetaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
%             IntanBehaviour.MIHitTrace(i).xgpbetaProbe = LFP.xgpbetaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
%             IntanBehaviour.MIHitTrace(i).wtbetaProbe = LFP.wtbetaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
%             IntanBehaviour.MIHitTrace(i).xfthetaProbe = LFP.xfthetaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
%             IntanBehaviour.MIHitTrace(i).xgpthetaProbe = LFP.xgpthetaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
%             IntanBehaviour.MIHitTrace(i).wtthetaProbe = LFP.wtthetaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
%             IntanBehaviour.MIHitTrace(i).xfgammaProbe = LFP.xfgammaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
%             IntanBehaviour.MIHitTrace(i).xgpgammaProbe = LFP.xgpgammaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
%             IntanBehaviour.MIHitTrace(i).wtgammaProbe = LFP.wtgammaProbe(:,:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
        end
    end
end

if isfield(IntanBehaviour,'MIFATrace')
    % Adding LFP an GP traces to  MIFATrace variable 
    for i=1:size(IntanBehaviour.  MIFATrace,2)
        IntanBehaviour.  MIFATrace(i).rawLFP = LFP.LFPdatacube(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        IntanBehaviour.  MIFATrace(i).xf = LFP.xf(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        IntanBehaviour.  MIFATrace(i).xgp = LFP.xgp(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        IntanBehaviour.  MIFATrace(i).wt = LFP.wt(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        if parameters.shank == 1
            IntanBehaviour.  MIFATrace(i).rawLFPProbe = LFP.LFPprobe(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
            IntanBehaviour.  MIFATrace(i).xfProbe = LFP.xfProbe(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
            IntanBehaviour.  MIFATrace(i).xgpProbe = LFP.xgpProbe(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
            IntanBehaviour.  MIFATrace(i).wtProbe = LFP.wtProbe(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        end
        if betaPresent == 1
            IntanBehaviour.  MIFATrace(i).xfbeta = LFP.xfbeta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
            IntanBehaviour.  MIFATrace(i).xgpbeta = LFP.xgpbeta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
            IntanBehaviour.  MIFATrace(i).wtbeta = LFP.wtbeta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        end
%         IntanBehaviour.  MIFATrace(i).xftheta = LFP.xftheta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).xgptheta = LFP.xgptheta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
%         IntanBehaviour.  MIFATrace(i).wttheta = LFP.wttheta(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        if gammaPresent == 1
            IntanBehaviour.  MIFATrace(i).xfgamma = LFP.xfgamma(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
            IntanBehaviour.  MIFATrace(i).xgpgamma = LFP.xgpgamma(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
            IntanBehaviour.  MIFATrace(i).wtgamma = LFP.wtgamma(:,:,IntanBehaviour.  MIFATrace(i).LFPIndex(1):IntanBehaviour.  MIFATrace(i).LFPIndex(end));
        end
        if parameters.shank == 1
%             IntanBehaviour.MIFATrace(i).xfbetaProbe = LFP.xfbetaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
%             IntanBehaviour.MIFATrace(i).xgpbetaProbe = LFP.xgpbetaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
%             IntanBehaviour.MIFATrace(i).wtbetaProbe = LFP.wtbetaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
%             IntanBehaviour.MIFATrace(i).xfthetaProbe = LFP.xfthetaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
%             IntanBehaviour.MIFATrace(i).xgpthetaProbe = LFP.xgpthetaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
%             IntanBehaviour.MIFATrace(i).wtthetaProbe = LFP.wtthetaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
%             IntanBehaviour.MIFATrace(i).xfgammaProbe = LFP.xfgammaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
%             IntanBehaviour.MIFATrace(i).xgpgammaProbe = LFP.xgpgammaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
%             IntanBehaviour.MIFATrace(i).wtgammaProbe = LFP.wtgammaProbe(:,:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
        end
    end
end