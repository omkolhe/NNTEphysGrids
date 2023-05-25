function [globalAvgSpectrogram, avgSpectrogramCWT,globalAvgVel,fwt] = getAvgSpectogram(xf,LFP,Encoder,trialTime,velTrial,parameters,flimit)


for trialno = 1:size(trialTime,1)
    xf1 = xf(:,:,trialTime(trialno,3):trialTime(trialno,4));
    relTime = Encoder.timeWindow2;
% %     figure();
%     for i=1:32
% %         subplot(parameters.rows,parameters.cols,i);
% %         if ismember(i,Intan.badChMap), continue; end
%         calCWTSpectogram(squeeze(xf1(floor((i-1)/parameters.cols)+1,mod(i-1,parameters.cols)+1,:)),relTime,LFP.Fs,20,flimit,0);
%     %     [s,f,t,ps,fc,tc] = spectrogram(squeeze(xf1(floor((i-1)/cols)+1,mod(i-1,cols)+1,:)),25,16,1024,1024,'yaxis','onesided');
%     %     imagesc(t,f(10:45),abs(ps(10:45,:)));hold on; xline(301/1024,'-r', 'LineWidth',2);set(gca,'YDir','normal') 
%     end
    
    % Average spectrogram across all channels
    for i=1:parameters.rows
        for j=1:parameters.cols
            [spectrogramCh((i-1)*parameters.cols + j,:,:) ,fwt] = calCWTSpectogram(squeeze(xf1(parameters.rows,parameters.cols,:)),relTime,LFP.Fs,20,flimit,0);
        end
    end
    avgSpectrogramCWT(trialno,:,:) = mean(spectrogramCh,1);
end
figure();
globalAvgSpectrogram = mean(avgSpectrogramCWT,1);
globalAvgVel = interp(mean(velTrial,1),LFP.Fs/Encoder.fs);
imagesc(relTime,fwt,squeeze(globalAvgSpectrogram));colormap('jet');set(gca,'YDir','normal');title('Wavelet based Average Spectogram');ylabel('Frequency (Hz)');xlabel('Time (s)');
c=colorbar;ylabel(c, 'Relative Power to white noise','FontSize',10);
hold on; yyaxis right; box off;
plot(Encoder.timeWindow2,globalAvgVel,'-w','LineWidth',2.5);
ylabel('Velocity (cm/s)');
end

