function [PGDfreqHit,PGDfreqMiss] = getPGDFreqBand(LFPdatacube,IntanBehaviour,parameters)

nMax = 19;
PGDfreqHit = zeros(nMax+1,1);
PGDfreqMiss = zeros(nMax+1,1);

for i=1:nMax
    xf = bandpass_filter(LFPdatacube,i*5,(i+1)*5,4,1000);
    [xgp, ~] = generalized_phase(xf,1000,0);
    PGDfreqHit1 = getPGD(xgp,IntanBehaviour.hitTrace,parameters);
    PGDfreqMiss1 = getPGD(xgp,IntanBehaviour.missTrace,parameters);
    PGDfreqHit(i) = mean(PGDfreqHit1,'all');
    PGDfreqMiss(i) = mean(PGDfreqMiss1,'all');
end

