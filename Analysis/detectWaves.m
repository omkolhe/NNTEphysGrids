function [Waves] = detectWaves(LFP,Encoder,paramerters)

spacing = paramerters.spacing;
rhoThres = paramerters.rhoThres;
X = paramerters.X;
Y = paramerters.Y;
for ii=1:size(Encoder.velTrig,2)
    Waves(ii).p = LFP.xgp(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4));
    Waves(ii).wt = LFP.wt(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4));
    %p = arrayfun(@(jj) inpaint_nans(p(:,:,jj)),1:size(p,3));
    Waves(ii).evaluationPoints = find_evaluation_points(Waves(ii).p,pi,0.2);
    %plot_evaluation_points( p, evaluationPoints );
    [Waves(ii).pm,Waves(ii).pd,Waves(ii).dx,Waves(ii).dy] = phase_gradient_complex_multiplication( Waves(ii).p, spacing );
    % Phase gradient directionality 
    [Waves(ii).PGD] = phase_gradient_directionality(Waves(ii).pm,Waves(ii).dx,Waves(ii).dy);
    % Wavelength
    Waves(ii).wl = 1./abs(Waves(ii).pm);
    % Instantaneous speed
    Waves(ii).insts = instantaneous_speed(Waves(ii).wt,Waves(ii).pm);
    Waves(ii).s = speedSpatial(Waves(ii).wt,Waves(ii).pm);
    % Velcity direction unit vector
    [Waves(ii).vx, Waves(ii).vy] = wavefront_direction(Waves(ii).pd,Waves(ii).insts);
    Waves(ii).velDir = atan2(Waves(ii).vy,Waves(ii).vx);
    % divergence calculation
    Waves(ii).source = find_source_points( Waves(ii).evaluationPoints, X, Y, Waves(ii).dx, Waves(ii).dy );
    % phase correlation with distance (\rho_{\phi,d} measure)
    Waves(ii).rho = zeros( 1, length(Waves(ii).evaluationPoints) );
    for jj = 1:length(Waves(ii).evaluationPoints)
        Waves(ii).ph = angle( Waves(ii).p(:,:,Waves(ii).evaluationPoints(jj)) );
        [Waves(ii).rho(jj),~,Waves(ii).D] = phase_correlation_distance( Waves(ii).ph, Waves(ii).source(:,jj), spacing );
    end
    % Thresholding using rhoThres. Reject all waves with rho less than
    % rhoThres
    indxBadWave = find(Waves(ii).rho<rhoThres);
    Waves(ii).evaluationPoints(indxBadWave) = [];
    Waves(ii).source(:,indxBadWave) = [];
    Waves(ii).rho(indxBadWave) = [];
    Waves(ii).D(indxBadWave) = [];
    Waves(ii).nWaves = size(Waves(ii).evaluationPoints,2);

    % Calculating the speed and wave direction of the detected waves in
    % that trial window
    for kk = 1:Waves(ii).nWaves
        st = Waves(ii).evaluationPoints(kk)-2;
        if st<0, st = 0; end  % adjusting is waves is detected at the start of window
        sp = Waves(ii).evaluationPoints(kk)+2;
        if sp>size(Waves(ii).p,3), sp = Waves(ii).evaluationPoints(kk); end % adjusting is waves is detected at the end of window
        %Waves(ii).speed(kk) = mean(abs(Waves(ii).insts(:,:,st:sp)),[1 2 3]); % speed in cm/s
        Waves(ii).speed(kk) = mean(abs(Waves(ii).s(st:sp)),'all'); % speed in cm/s
        Waves(ii).waveDir(kk) = mean(Waves(ii).velDir(st:sp),'all'); 
    end
    %Waves(ii).speedpdg = pgdMean(Waves(ii).PGD,Waves(ii).s,0.51);
    %Waves(ii).dirpdg = pgdMean(Waves(ii).PGD,Waves(ii).velDir,0.51);
    %plot_vector_field( exp( 1i .* Waves(ii).pd(:,:,100) ), 0 );
    %plot_wave_examples( LFP.xf(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4)), options, ii, Waves(ii).evaluationPoints, Waves(ii).source, Waves(ii).rho );
end

end

