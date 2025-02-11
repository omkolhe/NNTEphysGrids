%% Getting the link to the folder containing all CCA Mat files
% folderPath = uigetdir(path);
matFiles = dir(fullfile(folderPath, '*.mat'));

combinedCCA = [];

%% Loop through each .mat file
for file = 1:length(matFiles)
    % Get the full file path
    filePath = fullfile(matFiles(file).folder, matFiles(file).name);
    fprintf('Processing file: %s\n', matFiles(file).name);

    disp("Loading variables - IntanBehaviour");
    IntanBehaviour = load(filePath,"IntanBehaviour");
    IntanBehaviour = IntanBehaviour.IntanBehaviour;

    disp("Loading variables - M1Spikes");
    M1Spikes = load(filePath,"M1Spikes");
    M1Spikes = M1Spikes.M1Spikes;
    
    disp("Loading variables - M2Spikes");
    M2Spikes = load(filePath,"M2Spikes");
    M2Spikes = M2Spikes.M2Spikes;  

    parameters = IntanBehaviour.parameters;

    combinedCCA(file).IntanBehaviour = IntanBehaviour;
    combinedCCA(file).CCA = [];
    combinedCCA(file).argIn = [];
    combinedCCA(file).filePath = filePath;

    rt = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false));
    if mean(rt) < 0
        warning('Bad file')
        continue
    end
    
    % Removing disengage miss trials, hit trials with RT>1.5 and RT<0,
    % overlaping Hit and FA trials
    if ~isfield(IntanBehaviour,'cleanFlag')
        IntanBehaviour = cleanIntanBehaviour(IntanBehaviour);
    end

    %% Binning Spikes 
    if ~isfield('M1Spikes','spikes')
        disp('Bining M1 spikes ...')
        M1Spikes = binSpikes(M1Spikes,20000,1000);
    end
    if ~isfield('M2Spikes','spikes')
        disp('Bining M2 spikes ...')
        M2Spikes = binSpikes(M2Spikes,20000,1000);
    end
    
    if ~(M1Spikes.biningFs == M2Spikes.biningFs)
        disp('Bining Fs for M1 and M2 spikes is different')
        continue
    end
    
    %% Setting parameters for CCA Analysis 
    
    % The units of the arguments are with respect to the binning window used to bin spikes.
    argIn.BinWidth = 1000*(1/M1Spikes.biningFs);
    argIn.MaxDelay = 50;  
    argIn.TimeStep = 20;    
    argIn.WindowLength = 80;
    argIn.numCanonDim = 3;
    
    argIn.NumWorkers = Inf; % Requires Parallel Processing Toolbox

    combinedCCA(file).argIn = argIn;
    
    %% Calculating CCA 
    % Cue Hit
    spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.cueHitTrace,0);
    expCond = ones(size(spikes{1,1},3),1);
    
    disp(argIn)
    disp('Calculating CCA for cue Hit')
    combinedCCA(file).CCA.cueHitSpikes = spikes;
    combinedCCA(file).CCA.cueHit = ComputeCorrMap(spikes, expCond, argIn);
    
    % Cue Miss
    spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.cueMissTrace,0);
    expCond = ones(size(spikes{1,1},3),1);
    
    disp(argIn)
    disp('Calculating CCA for cue Miss')
    combinedCCA(file).CCA.cueMissSpikes = spikes;
    combinedCCA(file).CCA.cueMiss = ComputeCorrMap(spikes, expCond, argIn);
    
    % MI Hit
    spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.MIHitTrace,0);
    expCond = ones(size(spikes{1,1},3),1);
    
    disp(argIn)
    disp('Calculating CCA for MI Hit')
    combinedCCA(file).CCA.MIHitSpikes = spikes;
    combinedCCA(file).CCA.MIHit = ComputeCorrMap(spikes, expCond, argIn);
    
    % MI FA
    spikes = getCCASpikes(M1Spikes,M2Spikes,IntanBehaviour.MIFATrace,0);
    expCond = ones(size(spikes{1,1},3),1);
    
    disp(argIn)
    disp('Calculating CCA for MI FA')
    combinedCCA(file).CCA.MIFASpikes = spikes;
    combinedCCA(file).CCA.MIFA = ComputeCorrMap(spikes, expCond, argIn);

end

%% Saving the combined CCA 
sessionName = [folderPath,'/CombinedFiles/','CCACombined.mat'];
save(sessionName,"combinedCCA","folderPath","-v7.3"); 

