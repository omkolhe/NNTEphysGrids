function animateWaves(trial, Waves,saveOption,waveID)
% *SPONTANEOUS WAVES DEMO*
%
% PLOT WAVE EXAMPLES     plot specific examples of spontaneous waves based
%                          on calculated values of \rho_{\phi,d}
%
% INPUT
% x - datacube (rows,cols,timepts)
% trial - trial number
% evaluation_points - evaluation time points (cf. find_evaluation_points)
% source - putative source point
%
% OUTPUT
% animated spatiotemporal plot
%

x = Waves(trial).xf;
evaluation_points = Waves(trial).evaluationPoints;
source = Waves(trial).source;
vx = Waves(trial).vx;
vy = Waves(trial).vy;

% parameters
plot_pre_time = 5; pause_length = 0.4; 

% init
M = load( 'myMap.mat' );

if saveOption == 1
    fn = 'filename';
    writerobj = VideoWriter([fn '.avi'],'Uncompressed AVI'); % Initialize movie file
    writerobj.FrameRate = 5;
    open(writerobj);
end

ctr = 1; % wave detections counter
for jj = 1:length(evaluation_points)
    % get start and stop time, truncating if necessary
    st = Waves(trial).waveTime(jj,1) - plot_pre_time; sp = Waves(trial).waveTime(jj,2) + plot_pre_time;
    if ( st < 1 ), st = 1; end; if ( sp > size(x,3) ), sp = size(x,3); end

    % get data to plot, shuffle if option is chosen
    x_plot = x(:,:,st:sp);
    vx_plot = vx(1,st:sp);
    vy_plot = vy(1,st:sp);

    % create plot
    figure; title( sprintf( 'trial %d, wave example %d, 0 of %d ms', trial, ctr, size(x_plot,3) ) );
    color_range = [ min(reshape(x_plot,[],1)) max(reshape(x_plot,[],1)) ];
    h = imagesc( x_plot(:,:,1) ); hold on; axis image;
    plot( source(jj,1), source(jj,2), '.', 'markersize', 35, 'color', [.7 .7 .7] );
    h2 = quiver(source(jj,1), source(jj,2),vx_plot(1),vy_plot(1));
    h2.Color = 'White';
    h2.LineWidth = 2;
    set( gca, 'linewidth', 3, 'xtick', [], 'ytick', [], 'fontname', 'arial', 'fontsize', 16, 'ydir', 'reverse' );
    colormap( M.myMap ); box on; xlabel( 'electrodes' ); ylabel( 'electrodes' ); caxis( color_range )

    % create colorbar
    cb = colorbar();
    set( cb, 'location', 'southoutside' )
    set( cb, 'position', [0.6661    0.1674    0.2429    0.0588] );
    set( get(cb,'ylabel'), 'string', 'Amplitude (\muV)' ); set( cb, 'linewidth', 2 )
    if (saveOption == 1 && jj == waveID)
        writeVideo(writerobj,getframe(gcf)); %grabs current fig frame
    end

    % animate plot
    for kk = 1:size(x_plot,3)
        set( h, 'cdata', x_plot(:,:,kk) );
        set(h2, 'udata',vx_plot(1,kk),'vdata',vy_plot(1,kk));
        set( get(gca,'title'), 'string', ...
            sprintf( 'trial %d, wave example %d, %d of %d ms', trial, ctr, kk, size(x_plot,3) ) )
        pause(pause_length);
        if (saveOption == 1 && jj == waveID)
            writeVideo(writerobj,getframe(gcf)); %grabs current fig frame
        end
    end

    % increment counter
    ctr = ctr + 1;
    if (saveOption == 1 && jj == waveID)
        close(writerobj)
        disp('Video saved to current directory')
    end
end
