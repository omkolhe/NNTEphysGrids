parameters.Fs = 1000;
parameters.timeBeforeStim = 1.5;
parameters.timeAfterStim = 1.5;
[parameters.X,parameters.Y] = meshgrid( 1:parameters.cols, 1:parameters.rows );

for i=1:50%size(stimData,1)
    IntanBehaviour.stimControl(i).LFP = (downsample(squeeze(stimData(i,:,:))',4))';
    IntanBehaviour.stimControl(i).LFPdatacube = createDataCube2(IntanBehaviour.stimControl(i).LFP,parameters.rows,parameters.cols,[]); % Creating datacube);
    IntanBehaviour.stimControl(i).xf = bandpass_filter(IntanBehaviour.stimControl(i).LFPdatacube,5,40,4,1000);
    [IntanBehaviour.stimControl(i).xgp, IntanBehaviour.stimControl(i).wt] = generalized_phase(IntanBehaviour.stimControl(i).xf,1000,0);
end

%% Plotting waves
figure,stack_plot(reshape(IntanBehaviour.stimControl(34).xf,[],3001),0,4,parameters.Fs);
xline(1.5);

%% 
parameters.rhoThres = 0.65;


xf = arrayfun(@(s) s.xf, IntanBehaviour.stimControl, 'UniformOutput', false);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.stimControl, 'UniformOutput', false);
wt = arrayfun(@(s) s.wt, IntanBehaviour.stimControl, 'UniformOutput', false);
Waves.wavesHit = detectWaves2(xf,xgp,wt,IntanBehaviour.stimControl,parameters,parameters.rhoThres);

speedComb = horzcat(Waves.wavesHit(1:end).speed);
avgSpeed = mean(speedComb);
figure('Name','Histogram of wave speeds in Hits and Misses');
histfit(speedComb,100,'kernel');
xline(avgSpeed,'-r',{'Mean speed = ' num2str(avgSpeed) ' cm/s'});
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Hits');box off;


wavesHitPresent = vertcat(Waves.wavesHit.wavePresent);
wavesHitStart = vertcat(Waves.wavesHit.waveStart);

figure();
subplot(2,1,1);
title('Waves During Hit Trials')
rasterPlot(wavesHitPresent);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Stim','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);%ylim([1 20]);
set(gca,'TickDir','out','fontsize',14'); box off;
subplot(2,1,2)
bar((sum(wavesHitPresent,1)/size(IntanBehaviour.stimControl,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Stim','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);
ylabel('Wave probability');xlabel('Time (in ms)');set(gca,'TickDir','out','fontsize',14'); box off;


animateWaves(34,IntanBehaviour.stimControl,Waves.wavesHit,0,6);


function [LFPdatacube] = createDataCube2(LFP,rows,cols,chMap)
    
LFPdatacube = NaN(rows,cols,size(LFP,2));

for i=1:rows
    for j=1:cols
        if(~ismember((cols*(i-1)+j),chMap))
            LFPdatacube(i,j,:) = LFP(cols*(i-1)+j,:);
        end
    end
end

end


function [Waves] = detectWaves2(xf,xgp,wt1,behaviourTrace,parameters,threshold)

rhoThres = threshold;
dirThres = 15; %in degrees
X = parameters.X;
Y = parameters.Y;

waveTimeWindow = 40; % in points

for ii=1:50%size(behaviourTrace,2)    
    p = xgp{1,ii};
    wt = wt1{1,ii};
    Waves(ii).wavePresent = zeros(1,size(xf{1,ii},3));
    Waves(ii).waveStart = zeros(1,size(xf{1,ii},3));
%     for ll=1:size(Waves(ii).p,3)
%         xf{1,ii}(:,:,ll) = inpaint_nans(xf{1,ii}(:,:,ll),3);
%         Waves(ii).p(:,:,ll) = inpaint_nans(Waves(ii).p(:,:,ll),3);
%         Waves(ii).wt(:,:,ll) = inpaint_nans(Waves(ii).wt(:,:,ll),3);
%     end
%     p = arrayfun(@(jj) inpaint_nans(p(:,:,jj)),1:size(p,3));
    Waves(ii).evaluationPoints = find_evaluation_points(p,-pi/2,0.2);
%     plot_evaluation_points( Waves(ii).p, Waves(ii).evaluationPoints );
    [pm,pd,dx,dy] = phase_gradient_complex_multiplication(p, parameters.xspacing, parameters.yspacing );
    % Phase gradient directionality 
    [Waves(ii).PGD] = phase_gradient_directionality(pm,dx,dy);
    % Wavelength
    wl = 1./abs(pm);
    % Instantaneous speed
    insts = instantaneous_speed(wt,pm);
    s = speedSpatial(wt,pm);
    l = wavelengthSpatial(pm);
    % Velcity direction unit vector
    [vx, vy] = wavefront_direction(pd,insts);
    Waves(ii).vx = vx;
    Waves(ii).vy = vy;
    velDir = atan2(vy,vx);
    % divergence calculation
%     Waves(ii).source = find_source_points( Waves(ii).evaluationPoints, X, Y, Waves(ii).dx, Waves(ii).dy );
    % phase correlation with distance (\rho_{\phi,d} measure)
%     Waves(ii).rho = zeros( 1, length(Waves(ii).evaluationPoints) );
    Waves(ii).rho = {};
    for jj = 1:length(Waves(ii).evaluationPoints)
        st = Waves(ii).evaluationPoints(jj)-waveTimeWindow;
        if st<=0, st = 1; end  % adjusting is waves is detected at the start of window
        sp = Waves(ii).evaluationPoints(jj)+waveTimeWindow;
        if sp>size(p,3), sp = Waves(ii).evaluationPoints(jj); end % adjusting is waves is detected at the end of window
        rho = zeros(sp-st+1,1);
        sourcepoint = zeros(2,sp-st+1);
        dir = zeros(sp-st+1,1);
        for kk=1:(sp-st+1)
            ph = angle( p(:,:,(st+kk-1)));
            sourcepoint(:,kk) = find_source_points(st+kk-1,X,Y,dx, dy );
            [rho(kk,1),~,~] = phase_correlation_distance( ph,sourcepoint(:,kk), parameters.xspacing,parameters.yspacing );
            dir(kk) = rad2deg(velDir(st+kk-1));
        end
        waveClusters = findThresSeg(rho,rhoThres);
        % Rejecting segements less than 9ms
        clusterSize = diff(waveClusters,1,2);
        rejectIndex = find(clusterSize<9);
        waveClusters(rejectIndex,:) = [];
        waveClusters = removeNaNRows(waveClusters);
        % Selecting the wavecluster with maximum duration 
        % This is done to have just one wave at each evaluation point
        if isempty(waveClusters)
            Waves(ii).evaluationPoints(jj) = NaN;
%             Waves(ii).rho(jj) = NaN;
            Waves(ii).waveTime(jj,:) = [NaN,NaN];
            Waves(ii).source(jj,:) = [NaN,NaN];
            continue
        end
        clusterSize = diff(waveClusters,1,2);
        [~, clusterMaxIndx] = max(clusterSize);
        Waves(ii).waveTime(jj,:) = [st-1+waveClusters(clusterMaxIndx,1),st-1+waveClusters(clusterMaxIndx,2)];
        Waves(ii).evaluationPoints(jj) = Waves(ii).waveTime(jj,1);
        Waves(ii).source(jj,1) = round(mode(sourcepoint(1,waveClusters(clusterMaxIndx,1):waveClusters(clusterMaxIndx,2))));
        Waves(ii).source(jj,2) = round(mode(sourcepoint(1,waveClusters(clusterMaxIndx,1):waveClusters(clusterMaxIndx,2))));
        Waves(ii).rho{jj} = rho(waveClusters(clusterMaxIndx,1):waveClusters(clusterMaxIndx,2));
    end
    % Removing evaluation points in which no wave was detected 
    Waves(ii).evaluationPoints = Waves(ii).evaluationPoints(~isnan(Waves(ii).evaluationPoints));
    Waves(ii).waveTime = removeNaNRows(Waves(ii).waveTime);
    Waves(ii).source = removeNaNRows(Waves(ii).source);
    Waves(ii).rho = Waves(ii).rho(~cellfun(@isempty,Waves(ii).rho));
    % Removing duplicate detected waves
%     Waves(ii).evaluationPoints = unique((Waves(ii).evaluationPoints)','rows','stable');
%     Waves(ii).waveTime = unique(Waves(ii).waveTime,'rows','stable');
    % Merging repetetive waves 
    Waves(ii).nWaves = size(Waves(ii).waveTime,1);
    kk = 1;
    while kk < Waves(ii).nWaves
        if(Waves(ii).waveTime(kk,2) > Waves(ii).waveTime(kk+1,1)) % check there is a overlap
            if (Waves(ii).waveTime(kk,2) < Waves(ii).waveTime(kk+1,2))
                Waves(ii).rho{kk} = [Waves(ii).rho{kk};Waves(ii).rho{kk+1}(Waves(ii).waveTime(kk,2)-Waves(ii).waveTime(kk+1,1)+1:end)];
                Waves(ii).waveTime(kk,2) = Waves(ii).waveTime(kk+1,2); % Merge the two waves
            end
            Waves(ii).evaluationPoints(kk) = Waves(ii).waveTime(kk,1);
            Waves(ii).source(kk,:) = round(mean(Waves(ii).source(kk:kk+1,:),1));
            Waves(ii).waveTime(kk+1,:) = [];
            Waves(ii).source(kk+1,:) = [];
            Waves(ii).evaluationPoints(kk+1) = [];
            Waves(ii).nWaves = Waves(ii).nWaves-1;
            Waves(ii).rho(kk+1) = [];
        else
            kk = kk+1;
        end
    end
    
    %plot_evaluation_points( Waves(ii).p, Waves(ii).evaluationPoints );
    % Calculating the speed, amplitude and wave direction of the detected waves in
    % that trial window
    for kk = 1:Waves(ii).nWaves
        %Waves(ii).speed(kk) = mean(abs(Waves(ii).insts(:,:,st:sp)),[1 2 3]); % speed in cm/s
        Waves(ii).wavePresent(Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2)) = 1;
        Waves(ii).waveStart(Waves(ii).waveTime(kk,1)) = 1;
        Waves(ii).speed(kk) = mean(abs(s(Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2))),'all','omitnan'); % speed in cm/s
        Waves(ii).waveDir(kk) = atan2(mean(vy(Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2)),"all",'omitnan'),mean(vx(Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2)),"all",'omitnan'));
        Waves(ii).wavelength(kk) = mean(l(Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2)),'all','omitnan');
        Waves(ii).waveDuration(kk) = Waves(ii).waveTime(kk,2)-Waves(ii).waveTime(kk,1)+1;
        Waves(ii).waveAmp(kk) = max(abs(p(:,:,Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2))),[],"all",'omitnan')-min(abs(p(:,:,Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2))),[],"all",'omitnan');
    end
    %Waves(ii).speedpdg = pgdMean(Waves(ii).PGD,Waves(ii).s,0.51);
    %Waves(ii).dirpdg = pgdMean(Waves(ii).PGD,Waves(ii).velDir,0.51);
    %plot_vector_field( exp( 1i .* Waves(ii).pd(:,:,100) ), 0 );
    %plot_wave_examples( LFP.xf(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4)), options, ii, Waves(ii).evaluationPoints, Waves(ii).source, Waves(ii).rho );
end
  
end
