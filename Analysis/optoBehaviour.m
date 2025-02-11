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
outlierFlag = 0;
if outlierFlag == 0
    baselineRT = IntanBehaviourBaseline.reactionTime;
    optoRT = IntanBehaviourOpto.reactionTime;
else
    baselineRT = rmoutliers(IntanBehaviourBaseline.reactionTime);
    optoRT = rmoutliers(IntanBehaviourOpto.reactionTime);
end

[p,h] = ranksum(baselineRT,optoRT);
figure,plotBox2(baselineRT,optoRT);
xL=xlim;yL=ylim;
text(0.995*xL(2),0.995*yL(2),['p-val = ' num2str(p)],'HorizontalAlignment','right','VerticalAlignment','top')
ylabel('Reaction Time (s)'); title('M2 -> Th Opto');subtitle(['p-val = ' num2str(p)]);
xtix = {'Baseline','Opto'}; xtixloc = [1 2]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);set(gca,'TickDir','out','fontsize',14');
set(gca,'TickDir','out','fontsize',14');

%% 
