%% Loading variables from mat file 
if ~exist('IntanBehaviourBaseline','var')
    [enfile,enpath] = uigetfile('*.mat','Select mat file for session');
    if isequal(enfile,0)
       disp('User selected Cancel');
    else
       disp(['User selected ', fullfile(enpath,enfile)]);
    end
    IntanBehaviourBaseline = load(fullfile(enpath,enfile),"IntanBehaviourBaseline");
    IntanBehaviourBaseline = IntanBehaviourBaseline.IntanBehaviourBaseline;
    IntanBehaviourOpto = load(fullfile(enpath,enfile),"IntanBehaviourOpto");
    IntanBehaviourOpto = IntanBehaviourOpto.IntanBehaviourOpto;
    parameters = load(fullfile(enpath,enfile),"parameters");
    parameters = parameters.parameters;
end

%% Comparing Reaction Time 
outlierFlag = 1;
if outlierFlag == 0
    Behaviour.baselineRT = IntanBehaviourBaseline.reactionTime;
    Behaviour.optoRT = IntanBehaviourOpto.reactionTime;
else
    Behaviour.baselineRT = rmoutliers(IntanBehaviourBaseline.reactionTime);
    Behaviour.optoRT = rmoutliers(IntanBehaviourOpto.reactionTime);
end

[p,h] = ranksum(Behaviour.baselineRT,Behaviour.optoRT);
figure,plotBox2(Behaviour.baselineRT,Behaviour.optoRT);
xL=xlim;yL=ylim;
text(0.995*xL(2),0.995*yL(2),['p-val = ' num2str(p)],'HorizontalAlignment','right','VerticalAlignment','top')
ylabel('Reaction Time (s)'); title('M2 -> M1 Opto');subtitle(['p-val = ' num2str(p)]);
xtix = {'Baseline','Opto'}; xtixloc = [1 2]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);set(gca,'TickDir','out','fontsize',14');
set(gca,'TickDir','out','fontsize',14');

%% Hit rate and FA rate 

Behaviour.hitTimeBaseline = cell2mat(arrayfun(@(s) s.LFPtime(parameters.windowBeforePull*parameters.Fs+1), IntanBehaviourBaseline.cueHitTrace, 'UniformOutput', false));
Behaviour.hitTimeOpto = cell2mat(arrayfun(@(s) s.LFPtime(parameters.windowBeforePull*parameters.Fs+1), IntanBehaviourOpto.cueHitTrace, 'UniformOutput', false));
Behaviour.FATimeBaseline = cell2mat(arrayfun(@(s) s.LFPtime(parameters.windowBeforePull*parameters.Fs+1), IntanBehaviourBaseline.MIFATrace, 'UniformOutput', false));
Behaviour.FATimeOpto = cell2mat(arrayfun(@(s) s.LFPtime(parameters.windowBeforePull*parameters.Fs+1), IntanBehaviourOpto.MIFATrace, 'UniformOutput', false));
Behaviour.hitTime = [Behaviour.hitTimeBaseline Behaviour.hitTimeOpto];
Behaviour.FATime = [Behaviour.FATimeBaseline Behaviour.FATimeOpto];
Behaviour.timeResolution = 120 ; % in seconds
Behaviour.time = [0:Behaviour.timeResolution:floor(IntanBehaviourBaseline.time(end)/Behaviour.timeResolution)*Behaviour.timeResolution];
Behaviour.optoONTime = find(diff(IntanBehaviourBaseline.optoTrace,1,2)==1)/parameters.Fs; % in seconds

for i=1:size(Behaviour.time,2)-1
    Behaviour.hitrate(i) = numel(find(Behaviour.hitTime<Behaviour.time(i+1) & Behaviour.hitTime>Behaviour.time(i)))/(Behaviour.timeResolution/60); % per minute
    Behaviour.FArate(i) = numel(find(Behaviour.FATime<Behaviour.time(i+1) & Behaviour.FATime>Behaviour.time(i)))/(Behaviour.timeResolution/60); % per minute
end


figure();
Behaviour.plotTime = (0.5 * (Behaviour.time(1:end-1) + Behaviour.time(2:end)))/60;
h1=plot(Behaviour.plotTime,Behaviour.hitrate,'Color', [0 0.1 0.8],'LineWidth',2); hold on;
h2=plot(Behaviour.plotTime,Behaviour.FArate,'Color', [0.9 0.1 0.1],'LineWidth',2); hold on;
xline(Behaviour.optoONTime/60,'--r');
xlabel('Time (min)'); ylabel('Hit/FA rate (per min)');
legend([h1 h2],'Hit','FA','Location','best'); %ylim([5 15]);
title('Hit and FA rate for Opto')
box off;set(gca,'TickDir','out','fontsize',14');
drawnow;

%% Saving Behaviour 
% savepath = uigetdir(path);
sessionName = [enpath,'/','Behaviour.mat'];
save(sessionName,"Behaviour","savepath","parameters","-v7.3"); %