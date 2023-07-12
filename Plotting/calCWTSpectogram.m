function [powerCWT, fwt] = calCWTSpectogram(x,t,Fs,VoicesperOctave,flimit,plot)

%Ref - https://www.cosmos.esa.int/documents/1655127/1655136/Torrence_1998_Wavelet_Guide_BAMS.pdf/001d8327-b255-3024-a2f0-ce02e33ac98f

assert( numel(x)==numel(t), 'Number of elements in x and t dont match' );

fb = cwtfilterbank(SignalLength=numel(x),SamplingFrequency=Fs,Wavelet='amor',VoicesPerOctave=VoicesperOctave,FrequencyLimits=flimit);
[wt,fwt] = cwt(x,Filterbank=fb); % Calculate the CWT
powerCWT = (abs(wt).^2)/abs(var(x,1)); % Estimate power from CWT 
% The power calculated is normalzied/relative to the white noise power 
if plot==1 
    plotSpectrogram(powerCWT,t,fwt,'Wavelet Based Spectrogram','Time (s)','Frequency (Hz)')
end