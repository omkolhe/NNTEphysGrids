function [rhoThres] = getRhoThreshold(xgp,behaviourTrace,parameters,nShuffle,trialno,threshold)

% warningid = 'MATLAB:smoothn:SUpperBound';
warning('off','all');

spacing = parameters.spacing;
X = parameters.X;
Y = parameters.Y;

ii = trialno;
p = xgp(:,:,behaviourTrace(ii).LFPIndex(1):behaviourTrace(ii).LFPIndex(end));
evaluationPoints = find_evaluation_points(p,pi,spacing);
rho = zeros( nShuffle, length(evaluationPoints) );
for kk=1:nShuffle
    if kk==1
        pShuffle = p;
    else
        pShuffle = shuffle_channels(p);
    end
    [~,~,dx,dy] = phase_gradient_complex_multiplication( pShuffle, spacing );
    source = find_source_points( evaluationPoints, X, Y, dx, dy );
    for jj = 1:length(evaluationPoints)
        ph = angle( pShuffle(:,:,evaluationPoints(jj)) );
        rho(kk,jj) = phase_correlation_distance( ph, source(:,jj), spacing );
    end
end

% Plotting rho distribution
nEvalPoints = size(rho,2);
rho1 = reshape(rho',[1,nShuffle*nEvalPoints]);
rhoThres = prctile(rho1,threshold);
figure();histogram(rho1(1:nEvalPoints),'FaceColor','r'); hold on;
histogram(rho1(nEvalPoints+1:end),'FaceColor','b'); 
xline(rhoThres,'-r',{num2str(threshold) ' Percentile'});
xlabel(rho);
ylabel('Frequency');
title('Histogram of rho values for spatial shuffled electrodes');

warning('on','all');

end

