%% Random code
% %% Plotting video of phase 
% figure(); hold on;
% map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) );axis off;
% x = repmat(1:cols,rows,1); % generate x-coordinates
% y = repmat(1:rows,cols,1)'; % generate y-coordinates
% Ts = 1/1024;
% tvel = ((Ts:Ts:300*Ts))*1000;
% for i=1:numel(goodTrial.xgp)
%     ang = inpaint_nans(angle(goodTrial.xgp(:,:,i)));
% %     ang = inpaint_nans(goodTrial.xf(:,:,i));    
%     pause(0.1);
%     t= num2cell(rad2deg(ang));
%     t = cellfun(@num2str, t, 'UniformOutput',false);
% %     subplot(2,1,1);
%     imagesc(ang); title('Phase map for Electrode Map');
%     text(x(:),y(:),t,'HorizontalAlignment','Center');
% %     ax2 = axes; set( ax2, 'position', [0.02116    0.80976    0.0484    0.1000] ); axis image
% %     [x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
% %     h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )
% %     t1 = text( 0, 0, 'GP' );
% %     set( t1, 'fontname', 'arial', 'fontsize', 20, 'fontweight', 'bold', 'horizontalalignment', 'center' )
%     
% %     subplot(2,1,2)
% %     plot(Encoder.vel(1+Encoder.velTrig(137)-350:i+Encoder.velTrig(137)-350),'-b','LineWidth',2);xlim([tvel(1) tvel(end)]); ylim([-1 4]);ylabel('Velocity (in cm/s)');xlabel('Time (in ms)');
% end
% 
% %% Plotting video of LFP amplitude
% figure();
% x = repmat(1:cols,rows,1); % generate x-coordinates
% y = repmat(1:rows,cols,1)'; % generate y-coordinates
% tvel = Encoder.timeWindow1;
% for i=1:numel(goodTrial.xgp)
%     amp = inpaint_nans(goodTrial.xf(:,:,i));    
%     pause(0.1);
% %     subplot(2,1,1);
%     imagesc(amp); title('Phase map for Electrode Map');colorbar;
% %     if mod(i,2)==0, continue; end
% %     subplot(2,1,2)
% %     plot(tvel(1:floor(i/2)),Encoder.vel(1,Encoder.trialTime(trialno,1):Encoder.trialTime(trialno,1)+floor(i/2)),'-b','LineWidth',2)
% %     plot(tvel(1:floor(i/2)),'-b','LineWidth',2);
% %     xlim([tvel(1) tvel(end)]); ylim([-1 4]);ylabel('Velocity (in cm/s)');xlabel('Time (in s)');
% end
% %% Finding source points 
% % plot_wave_examples( goodTrial.xf, options, jj, evaluationPoints, source, rho );
% 
% figure();
% for i=1:numel(goodTrial.relTime)
%     quiver(dx(:,:,i),dy(:,:,i));
%     pause(0.1);
%     ylim([0 9]);xlim([0 5]);
% end
% 
% %% 
% 
% dt = 1 / 1024; T = size(LFP.LFP,2) / Fs; time = dt:dt:T;
time = LFP.times(Encoder.trialTime(4,3):Encoder.trialTime(4,4));
xw = squeeze(LFP.xf(4,3,Encoder.trialTime(4,3):Encoder.trialTime(4,4)));
xwRaw = squeeze(LFP.LFPdatacube(4,3,Encoder.trialTime(4,3):Encoder.trialTime(4,4)));
xgp1 = squeeze(LFP.xgp(4,3,Encoder.trialTime(4,3):Encoder.trialTime(4,4)));
% main figure
fg1 = figure; hold on; ax1 = gca; 
plot( time, xwRaw, 'linewidth', 1, 'color', 'k' ); h4 = cline( time, xw, [], angle(xgp1) );
xlim([166 166.5])
set( h4, 'linestyle', '-', 'linewidth', 2  ), axis off
l1 = line( [.1 .2], [-125 -125] ); set( l1, 'linewidth', 4, 'color', 'k' )
l2 = line( [.1 .1], [-125 -75] ); set( l2, 'linewidth', 4, 'color', 'k' )

% inset
map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )
ax2 = axes; set( ax2, 'position', [0.2116    0.6976    0.0884    0.2000] ); axis image
[x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )
%%
% text labels
t1 = text( 0, 0, 'GP' );
set( t1, 'fontname', 'arial', 'fontsize', 28, 'fontweight', 'bold', 'horizontalalignment', 'center' )
set( gcf, 'currentaxes', ax1 )
t2 = text( 0.1260, -146.8832, '100 ms' );
set( t2, 'fontname', 'arial', 'fontsize', 24, 'fontweight', 'bold' )
t2 = text( 0.0852, -130.2651, '50 \muV' );
set( t2, 'fontname', 'arial', 'fontsize', 24, 'fontweight', 'bold', 'rotation', 90 )
% 
% %% 
% 
% waveden = zeros(1,size(LFP.LFP,2));
% waveden(evaluationPoints) = 1;
% waveden1 = downsample(waveden,2);
% figure(),plot(waveden1);
% 
% 
% figure();
% plot(Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
% hold on;
% plot(waveden1);
% 
% 
% 
% %meanPhase is the mean phase map over a 5ms window starting at the wave
% %initiaion time
% [rho(jj),~,D,pl] = phase_correlation_distance(meanPhase,source(:,jj), parameters.pixel_spacing);
% %fit the phase (p1) vs distance (D) to a line
% pfit = polyfit(D,pl,1);
% %use phase vs distance to calculate wave speed
% k = abs(pfit(1));%rad/mm, WAVE NUMBER, slope of the lin fit (spatial derivative k=dPhase/dx)
% IF = nanmean(IF,3);%Hz, INSTANTANEOUS FREQUENCY for each electrode in the small time window. f = dPhase/dt
% IF = nanmean(IF(:));%Hz, AVG INSTANTANEOUS FREQUENCY across the grid
% w = IF*2*pi;%(1/s)*(rad) = rad/s, AVG INSTANTANEOUS ANGULAR FREQUENCY across the grid. w = 2pi*f
% speed(jj) = (w/k).*(1/1000); %(rad/s)/(rad/mm) = (mm/s)*(1m/1000mm) = m/s, speed
% 


