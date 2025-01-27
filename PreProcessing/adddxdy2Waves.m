function Waves = adddxdy2Waves(behaviourTrace,xgp,wt1,Waves,parameters)
    for ii=1:size(behaviourTrace,2)    
        p = xgp{1,ii};
        wt = wt1{1,ii};
        if sum(isnan(p),'all') >= 1
            for ll=1:size(p,3)
        %         xf{1,ii}(:,:,ll) = inpaint_nans(xf{1,ii}(:,:,ll),3);
                p(:,:,ll) = inpaint_nans(p(:,:,ll),3);
                wt(:,:,ll) = inpaint_nans(wt(:,:,ll),3);
            end
        end
        [pm,pd,dx,dy] = phase_gradient_complex_multiplication(p, parameters.xspacing, parameters.yspacing );
        Waves(ii).dx = dx;
        Waves(ii).dy = dy;
    end
end

% xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.cueHitTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourBaseline.cueHitTrace, 'UniformOutput', false);
% WavesBaseline.wavesHit = adddxdy2Waves(IntanBehaviourBaseline.cueHitTrace,xgp,wt,WavesBaseline.wavesHit,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.cueMissTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourBaseline.cueMissTrace, 'UniformOutput', false);
% WavesBaseline.wavesMiss = adddxdy2Waves(IntanBehaviourBaseline.cueMissTrace,xgp,wt,WavesBaseline.wavesMiss,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.hitTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourBaseline.hitTrace, 'UniformOutput', false);
% WavesBaseline.wavesHitReward = adddxdy2Waves(IntanBehaviourBaseline.hitTrace,xgp,wt,WavesBaseline.wavesHitReward,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.missTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourBaseline.missTrace, 'UniformOutput', false);
% WavesBaseline.wavesFA = adddxdy2Waves(IntanBehaviourBaseline.missTrace,xgp,wt,WavesBaseline.wavesFA,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.MIHitTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourBaseline.MIHitTrace, 'UniformOutput', false);
% WavesBaseline.wavesMIHit = adddxdy2Waves(IntanBehaviourBaseline.MIHitTrace,xgp,wt,WavesBaseline.wavesMIHit,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.MIFATrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourBaseline.MIFATrace, 'UniformOutput', false);
% WavesBaseline.wavesMIFA = adddxdy2Waves(IntanBehaviourBaseline.MIFATrace,xgp,wt,WavesBaseline.wavesMIFA,parameters);
% 
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourOpto.cueHitTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourOpto.cueHitTrace, 'UniformOutput', false);
% WavesOpto.wavesHit = adddxdy2Waves(IntanBehaviourOpto.cueHitTrace,xgp,wt,WavesOpto.wavesHit,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourOpto.cueMissTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourOpto.cueMissTrace, 'UniformOutput', false);
% WavesOpto.wavesMiss = adddxdy2Waves(IntanBehaviourOpto.cueMissTrace,xgp,wt,WavesOpto.wavesMiss,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourOpto.hitTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourOpto.hitTrace, 'UniformOutput', false);
% WavesOpto.wavesHitReward = adddxdy2Waves(IntanBehaviourOpto.hitTrace,xgp,wt,WavesOpto.wavesHitReward,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourOpto.missTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourOpto.missTrace, 'UniformOutput', false);
% WavesOpto.wavesFA = adddxdy2Waves(IntanBehaviourOpto.missTrace,xgp,wt,WavesOpto.wavesFA,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourOpto.MIHitTrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourOpto.MIHitTrace, 'UniformOutput', false);
% WavesOpto.wavesMIHit = adddxdy2Waves(IntanBehaviourOpto.MIHitTrace,xgp,wt,WavesOpto.wavesMIHit,parameters);
% 
% xgp = arrayfun(@(s) s.xgp, IntanBehaviourOpto.MIFATrace, 'UniformOutput', false);
% wt = arrayfun(@(s) s.wt, IntanBehaviourOpto.MIFATrace, 'UniformOutput', false);
% WavesOpto.wavesMIFA = adddxdy2Waves(IntanBehaviourOpto.MIFATrace,xgp,wt,WavesOpto.wavesMIFA,parameters);