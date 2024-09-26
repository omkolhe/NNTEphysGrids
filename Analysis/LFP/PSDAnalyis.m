[PSD.ChHit , PSD.f] = getAvgPSD(IntanBehaviour.cueHitTrace,parameters);
[PSD.ChMiss , PSD.f] = getAvgPSD(IntanBehaviour.cueMissTrace,parameters);
PSD.trialPSDHit = squeeze(10*log10(mean(PSD.ChHit,2,"omitnan")));
PSD.trialPSDMiss = squeeze(10*log10(mean(PSD.ChMiss,2,"omitnan")));

[PSD.ChHitReward , PSD.f] = getAvgPSD(IntanBehaviour.hitTrace,parameters);
[PSD.ChFA , PSD.f] = getAvgPSD(IntanBehaviour.missTrace,parameters);
PSD.trialPSDHitReward = squeeze(10*log10(mean(PSD.ChHitReward,2,"omitnan")));
PSD.trialPSDFA = squeeze(10*log10(mean(PSD.ChFA,2,"omitnan")));

[PSD.ChHitMI , PSD.f] = getAvgPSD(IntanBehaviour.MIHitTrace,parameters);
[PSD.ChFAMI , PSD.f] = getAvgPSD(IntanBehaviour.MIFATrace,parameters);
PSD.trialPSDHitMI = squeeze(10*log10(mean(PSD.ChHitMI,2,"omitnan")));
PSD.trialPSDFAMI = squeeze(10*log10(mean(PSD.ChFAMI,2,"omitnan")));

remove_artifact = 0;
if remove_artifact == 1
    % Removing trials with artifacts
    PSD.rejectThres = 35; % in db
    PSD.rejectFreq = 3; %  not Hz but index in PSD.f
    PSD.artifactTrialIndex = find(PSD.trialPSDHit(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of Hit trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDHit(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.cueHitTrace(PSD.artifactTrialIndex)=[];
    PSD.ChHit(PSD.artifactTrialIndex,:,:) = [];
    IntanBehaviour.reactionTime(PSD.artifactTrialIndex,:,:) = [];
    
    PSD.artifactTrialIndex = find(PSD.trialPSDMiss(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of Miss trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDMiss(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.cueMissTrace(PSD.artifactTrialIndex)=[];
    PSD.ChMiss(PSD.artifactTrialIndex,:,:) = [];
    
    PSD.artifactTrialIndex = find(PSD.trialPSDHitReward(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of Hit trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDHitReward(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.hitTrace(PSD.artifactTrialIndex)=[];
    PSD.ChHitReward(PSD.artifactTrialIndex,:,:) = [];
    
    PSD.artifactTrialIndex = find(PSD.trialPSDFA(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of FA trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDFA(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.missTrace(PSD.artifactTrialIndex)=[];
    PSD.ChFA(PSD.artifactTrialIndex,:,:) = [];

    PSD.artifactTrialIndex = find(PSD.trialPSDHitMI(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of Hit MI trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDHitMI(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.MIHitTrace(PSD.artifactTrialIndex)=[];
    PSD.ChHitMI(PSD.artifactTrialIndex,:,:) = [];

    PSD.artifactTrialIndex = find(PSD.trialPSDFAMI(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of FA MI trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDFAMI(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.MIFATrace(PSD.artifactTrialIndex)=[];
    PSD.ChFAMI(PSD.artifactTrialIndex,:,:) = [];
end

PSD.avgPSDHit = squeeze(10*log10(mean(PSD.ChHit,[1 2],"omitnan")));
PSD.avgPSDMiss = squeeze(10*log10(mean(PSD.ChMiss,[1 2],"omitnan")));
PSD.avgPSDHitReward = squeeze(10*log10(mean(PSD.ChHitReward,[1 2],"omitnan")));
PSD.avgPSDFA = squeeze(10*log10(mean(PSD.ChFA,[1 2],"omitnan")));
PSD.avgPSDHitMI = squeeze(10*log10(mean(PSD.ChHitMI,[1 2],"omitnan")));
PSD.avgPSDFAMI = squeeze(10*log10(mean(PSD.ChFAMI,[1 2],"omitnan")));

figure();
subplot(1,2,1);
plot(PSD.f(1:81),PSD.trialPSDHit(:,1:81),'Color', [0 0 1 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDHit(1:81),'Color', [0 0 1 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average PSD for Hit Trials');
box off;

subplot(1,2,2);
plot(PSD.f(1:81),PSD.trialPSDMiss(:,1:81),'Color', [1 0 0 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDMiss(1:81),'Color', [1 0 0 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average PSD for Miss Trials');
box off;

figure();
subplot(1,2,1);
plot(PSD.f(1:81),PSD.trialPSDHitReward(:,1:81),'Color', [0 0 1 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDHitReward(1:81),'Color', [0 0 1 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average PSD for Hit Trials');
box off;

subplot(1,2,2);
plot(PSD.f(1:81),PSD.trialPSDFA(:,1:81),'Color', [1 0 0 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDFA(1:81),'Color', [1 0 0 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average PSD for FA Trials');
box off;

figure();
subplot(1,2,1);
plot(PSD.f(1:81),PSD.trialPSDHitMI(:,1:81),'Color', [0 0 1 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDHitMI(1:81),'Color', [0 0 1 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average PSD for Hit MI Trials');
box off;

subplot(1,2,2);
plot(PSD.f(1:81),PSD.trialPSDFAMI(:,1:81),'Color', [1 0 0 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDFAMI(1:81),'Color', [1 0 0 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average PSD for FA MI Trials');
box off;