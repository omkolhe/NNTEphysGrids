function [globalAvgSpectrogram, avgSpectrogramCWT,globalAvgBehaviour,fwt] = getAvgSpectogram(xf,LFPFs,behaviourTrace,parameters,flimit)

%% Inputs
% xf - filtered signal in datacube format
% LFPFs - sampling rate of LFP/xf
% behaviourTrace - the struct containing lever traces and the timing
% details arranged as trials
% parameters - paramaters files
% flimit - the frequency limits for the spectogram

voicesPerOctave = 20;

for trialno = 1:size(behaviourTrace,2)
    xf1 = xf(:,:,behaviourTrace(trialno).LFPIndex(1):behaviourTrace(trialno).LFPIndex(end));
    relTime = behaviourTrace(trialno).time;
    xf1Avg = squeeze(mean(xf1,[1 2]));
    [avgSpectrogramCWT(trialno,:,:),fwt] = calCWTSpectogram(xf1Avg,relTime,LFPFs,voicesPerOctave,flimit,0);
    % Average spectrogram across all channels
%     for i=1:parameters.rows
%         for j=1:parameters.cols
%             [spectrogramCh((i-1)*parameters.cols + j,:,:) ,fwt] = calCWTSpectogram(squeeze(xf1(parameters.rows,parameters.cols,:)),relTime,LFPFs,voicesPerOctave,flimit,0);
%         end
%     end
%     avgSpectrogramCWT(trialno,:,:) = mean(spectrogramCh,1);
end
figure();
globalAvgSpectrogram = mean(avgSpectrogramCWT,1);
globalAvgBehaviour = mean(horzcat(behaviourTrace(1:end).trace),2);
plotSpectrogram(squeeze(globalAvgSpectrogram),relTime,fwt);
hold on; yyaxis right; box off;
plot(relTime,globalAvgBehaviour,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)');
end

