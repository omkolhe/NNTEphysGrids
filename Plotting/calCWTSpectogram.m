function [powerCWT, fwt] = calCWTSpectogram(x,t,Fs,VoicesperOctave,flimit,plot,baselinesub)

%Ref - https://www.cosmos.esa.int/documents/1655127/1655136/Torrence_1998_Wavelet_Guide_BAMS.pdf/001d8327-b255-3024-a2f0-ce02e33ac98f

if ~exist('baselinesub','var')
    BLsub = 0;
else
    BLsub = baselinesub;
    nBLSub = 0.25*Fs; % in number of points
end

assert( numel(x)==numel(t), 'Number of elements in x and t dont match' );

fb = cwtfilterbank(SignalLength=numel(x),SamplingFrequency=Fs,Wavelet='amor',VoicesPerOctave=VoicesperOctave,FrequencyLimits=flimit);
[wt,fwt] = cwt(x,Filterbank=fb); % Calculate the CWT
powerCWT = (abs(wt).^2);%/abs(var(x,1)); % Estimate power from CWT 

if BLsub == 1
    BLPower = mean(powerCWT(:,1:nBLSub),2);
    powerCWTBLsub = bsxfun(@rdivide,powerCWT,BLPower);
    powerCWT = powerCWTBLsub;
end


if (plot==1 )
    plotSpectrogram(10*log10(powerCWT),t,fwt,'Wavelet Based Spectrogram','Time (s)','Frequency (Hz)')
end


