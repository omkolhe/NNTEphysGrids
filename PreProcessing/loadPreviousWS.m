[wsfile,wspath] = uigetfile('*.mat');
if isequal(wsfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(wspath,wsfile)]);
end

load(fullfile(wspath,wsfile));