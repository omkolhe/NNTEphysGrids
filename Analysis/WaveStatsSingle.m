function [] = WaveStatsSingle(Waves,parameters,plot)

avgSpeed = mean(Waves.speed,'all');
if plot == 1
    figure('Name','Histogram of wave speeds');
    histfit(Waves.speed,100,'kernel');
    xline(avgSpeed,'-r',{'Mean speed = ' num2str(avgSpeed) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed');
end

avgl = mean(Waves.wavelength,'all');
if plot == 1
    figure('Name','Histogram of wavelengths')
    histfit(Waves.wavelength,100,'kernel');
    xline(avgl,'-r',{'Mean wavelength = ' num2str(avgSpeed) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength');
end

% Wave direction stats

if plot == 1
    figure('Name','Polar Histogram for wave direction');
    polarhistogram(Waves.waveDir,30);
    title('Wave Direction');
end

end

