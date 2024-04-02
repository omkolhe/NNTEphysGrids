Spikes = SpikesBaseline;
IntanBehaviour = IntanBehaviourBaseline;
%% Removing clusters that are not correct
Spikes = rejectSpikes(Spikes,0.75,1,parameters);

%% Plotting spikes raster and FR to check 
figure,hold on
for n = 53
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.hit.spks{n}),axis off
subplot(2,1,[2]),bar(-1500:1500,smoothdata(Spikes.PSTH.hit.spkRates(n,:)),'FaceColor',[28/255 117/255 188/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out');sgtitle(['Neuron - ' num2str(n)])
set(gca,'fontsize',16)
xline(0)
% xline(mean(IntanBehaviour.reactionTime)*1000)
end
%% Segmenting spiking data according to behaviour
FsScaling = (Spikes.Clusters(1).cluster(1)/Spikes.Clusters(1).spikeTime(1))/parameters.Fs;
for i=1:size(Spikes.Clusters,2)
    Spikes.Clusters(i).spikeIndexLFP = round(Spikes.Clusters(i).cluster/FsScaling);
end
% Making spikes rasters 
Spikes.spikeRaster = zeros(size(Spikes.Clusters,2),size(IntanBehaviour.time,2));
for i=1:size(Spikes.Clusters,2)
    Spikes.spikeRaster(i,Spikes.Clusters(i).spikeIndexLFP) = 1;
end

% Adding spike raster to IntanBehaviour
for i = 1:size(IntanBehaviour.cueHitTrace,2)
    IntanBehaviour.cueHitTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.cueMissTrace,2)
    IntanBehaviour.cueMissTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.hitTrace,2)
    IntanBehaviour.hitTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.missTrace,2)
    IntanBehaviour.missTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.MIHitTrace,2)
    IntanBehaviour.MIHitTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.MIFATrace,2)
    IntanBehaviour.MIFATrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
end

%% Getting average phase map for each spike
for i=1:Spikes.nSpikes
    Spikes.spikeTrigPhase(i).rawPhaseMap = circ_mean(angle(LFP.xgp(:,:,Spikes.Clusters(i).spikeIndexLFP)),[],3);
    a = 1*exp(1i*Spikes.spikeTrigPhase(i).rawPhaseMap);
    Spikes.spikeTrigPhase(i).phaseMap = angle(inpaint_nans(a,3));
    [Spikes.spikeTrigPhase(i).pm,Spikes.spikeTrigPhase(i).pd,Spikes.spikeTrigPhase(i).dx,Spikes.spikeTrigPhase(i).dy] = getPhaseGradient( Spikes.spikeTrigPhase(i).phaseMap, parameters.xspacing,parameters.yspacing);
    Spikes.spikeTrigPhase(i).sourcePoint = find_source_points2(1,parameters.X,parameters.Y,Spikes.spikeTrigPhase(i).dx,Spikes.spikeTrigPhase(i).dy);
    Spikes.spikeTrigPhase(i).PGD = phase_gradient_directionality(Spikes.spikeTrigPhase(i).pm,Spikes.spikeTrigPhase(i).dx,Spikes.spikeTrigPhase(i).dy);
    [Spikes.spikeTrigPhase(i).rho,~,~] = phase_correlation_distance( Spikes.spikeTrigPhase(i).phaseMap,Spikes.spikeTrigPhase(i).sourcePoint, parameters.xspacing,parameters.yspacing );
    [Spikes.spikeTrigPhase(i).vx, Spikes.spikeTrigPhase(i).vy] = wavefront_direction(Spikes.spikeTrigPhase(i).pd,[]);
    Spikes.spikeTrigPhase(i).waveDir = atan2(Spikes.spikeTrigPhase(i).vy,Spikes.spikeTrigPhase(i).vx);
end
%% Plotting spike triggered phase maps
figure,
n = 10;
imagesc(Spikes.spikeTrigPhase(n).phaseMap);
map = colorcet( 'C2' );
map = circshift(map,1);
colormap(map)
c = colorbar;
hold on;
[XX,YY] = meshgrid( 1:size(Spikes.spikeTrigPhase(n).dx+1i*Spikes.spikeTrigPhase(n).dy,2), 1:size(Spikes.spikeTrigPhase(n).dx+1i*Spikes.spikeTrigPhase(n).dy,1) );
M = real( exp( 1i * angle(Spikes.spikeTrigPhase(n).dx+1i*Spikes.spikeTrigPhase(n).dy) ) ); N = imag( exp( 1i * angle(Spikes.spikeTrigPhase(n).dx+1i*Spikes.spikeTrigPhase(n).dy) ) );
quiver( XX, YY, M, N, 0.25, 'r' );
scatter(Spikes.spikeTrigPhase(n).sourcePoint(1),Spikes.spikeTrigPhase(n).sourcePoint(2),'filled');
%% Plotting Phase Gradient directions for the spike triggered phase maps
PMGDirection = cell2mat(arrayfun(@(s) s.waveDir, Spikes.spikeTrigPhase,'UniformOutput',false));
figure,plotDirectionHistogram(PMGDirection,36,[]);
phaseBoundary1 = pi/2;
PMGDirectionRotated = angle(exp(1i*PMGDirection)*exp(1i*-phaseBoundary1));
figure,plotDirectionHistogram(PMGDirectionRotated,36,[]);
for i=1:Spikes.nSpikes
    if (PMGDirectionRotated(i)>0)
        Spikes.spikeTrigPhase(i).directionCluster = 1;
    else
        Spikes.spikeTrigPhase(i).directionCluster = 2;
    end
end