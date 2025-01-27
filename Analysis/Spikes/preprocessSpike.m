function output = preprocessSpike(data,Fs)
Fc = [500 2000];
Wn = Fc./(Fs/2);
b = fir1(5000,Wn,'bandpass');


% Fc = [250];
% Wn = Fc./(Fs/2);
% b = fir1(5000,Wn,'high');
a = 1;
% Define filter specifications
% order = 2;             % Filter order (choose based on desired roll-off)
% F_low = 100;            % Lower cutoff frequency in Hz
% F_high = 6000;          % Upper cutoff frequency in Hz

% % Normalize cutoff frequencies to Nyquist frequency
% nyquist = Fs / 2;
% % Wn = [F_low, F_high] / nyquist;
% %
% % % Design bandpass Bessel filter
% % [b, a] = besself(order, Wn, 'bandpass');
% 
% Wn = [F_low] / nyquist;
% % Design bandpass Bessel filter
% [b, a] = butter(order, Wn, 'high');


rawspikeTrace = filtfilt(b,a,double(data)');
rawspikeTrace = rawspikeTrace';
commonModeAvg = rawspikeTrace-mean(rawspikeTrace,1,"omitnan");
whitenedSpikeTrace = commonModeAvg;
output.rawspikeTrace = rawspikeTrace;
output.whitenedSpikeTrace = whitenedSpikeTrace - mean(whitenedSpikeTrace,2);
end