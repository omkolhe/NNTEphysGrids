function [PGD] = getPGD(xgp,behaviourTrace,parameters)

spacing = parameters.spacing;

totaltime = parameters.windowAfterPull + parameters.windowBeforePull + parameters.ts;
time = [parameters.ts:parameters.ts:totaltime]; 

PGD = zeros(size(behaviourTrace,2),size(time,2));

for ii=1:size(behaviourTrace,2)
    p = xgp(:,:,behaviourTrace(ii).LFPIndex(1):behaviourTrace(ii).LFPIndex(end));
    [pm,~,dx,dy] = phase_gradient_complex_multiplication(p, spacing );
    % Phase gradient directionality 
    PGD(ii,:) = phase_gradient_directionality(pm,dx,dy); 
end
  
end

