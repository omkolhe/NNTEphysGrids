% Concatenate trials

addpath(genpath('main'));
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.rhd')); %Parses RHD files
count = 1;
downsampleRate = 4;
targetedFs = 10000;
L = length(directory);
for idx = 1:L
    file = directory(idx).folder;
    path = directory(idx).name;
    Intan = read_Intan_RHD2000_file(file,path); 
    Fs =  Intan.frequency_parameters.amplifier_sample_rate;
    allIntan{count} = resample(Intan.amplifier_data',targetedFs,Fs);
    if ~isempty(Intan.board_adc_data)
        analog_adc_data{count} = resample(Intan.board_adc_data',targetedFs,Fs);
    else
        analog_adc_data = [];
    end
    if ~isempty(Intan.board_dig_in_data)
        dig_in_data{count} = downsample(Intan.board_dig_in_data',round(Fs/targetedFs),1);
    else
        dig_in_data = [];
    end
    count = count+1;
end % load Intan files
Fs = Intan.frequency_parameters.amplifier_sample_rate;
% Concatenate intan files for the whole session
disp('Combining...')
Intan.allIntan = vertcat(allIntan{:})';

if ~isempty(analog_adc_data)
    Intan.analog_adc_data = vertcat(analog_adc_data{:})';
else
    Intan.analog_adc_data = [];
end

if ~isempty(dig_in_data)
    Intan.dig_in_data = vertcat(dig_in_data{:})';
else
    Intan.dig_in_data = [];
end

disp('Compressing...')
Intan.allIntan = single(Intan.allIntan);
Intan.analog_adc_data = single(Intan.analog_adc_data);
Intan.dig_in_data = single(Intan.dig_in_data);
% Adjust electrode order by depth
% UCLA_probe_map %legacy file call
Intan.allIntan  = Intan.allIntan(electrode_map,:);
% Fix recording offset
Intan.offset = 1; % second
Intan.offsetSample = targetedFs*Intan.offset;
disp(['Adjusting for ' num2str(Intan.offset) ' second offset']);
Intan.allIntan = Intan.allIntan(:,Intan.offsetSample:(size(Intan.allIntan,2)-Intan.offsetSample));
if ~isempty(Intan.analog_adc_data)
    Intan.analog_adc_data = Intan.analog_adc_data(:,Intan.offsetSample:(size(Intan.analog_adc_data,2)-Intan.offsetSample));
else
    Intan.analog_adc_data = [];
end
if ~isempty(Intan.dig_in_data)
    Intan.dig_in_data = Intan.dig_in_data(:,Intan.offsetSample:(size(Intan.dig_in_data,2)-Intan.offsetSample));
else
    Intan.dig_in_data = [];
end
clear amplifier_data t_amplifier frequncy_parameters notes aux_input_channels...
    aux_input_data board_dig_in_channels board_dig_in_data amplifier_channels...
    board_adc_data board_adc_channels t_board_adc t_dig t_aux_input analog_adc_data dig_in_data allIntan


