%%Time specifications:
Fs = 1000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 0.7;             % seconds
t1 = (0:dt:StopTime-dt)';     % seconds
Anoise = 0;
%%Sine wave:
Fc1 = 15;% hertz
Fc2 = 65;% hertz
noise = Anoise*(randn(length(t1),1));
x1 = 2*cos(2*pi*Fc1*t1) + 3*cos(2*pi*Fc2*t1) + noise ;
% x1 = 2*cos(2*pi*Fc1*t1) + noise ;

StopTime2 = 2;             % seconds
t2 = (StopTime:dt:StopTime2-dt)';     % seconds
Fc1 = 11;% hertz
Fc2 = 76;% hertz
noise = Anoise*(randn(length(t2),1));
x2 = 2.5*cos(2*pi*Fc1*t2) + 3*cos(2*pi*Fc2*t2) + noise ;
% x2 = 2*cos(2*pi*Fc1*t2) + noise ;

t = [t1;t2];
x = [x1;x2];
% Plot the signal versus time:
figure;
plot(t,x);
xlabel('time (in seconds)');
title('Signal versus Time');
zoom xon;

figure();
[avgSpectrogramCWT,fwt] = calCWTSpectogram(x,t,Fs,20,[5 90],1);
