function [VocTypeHD]=VocTypeHD_Def_4Matfiles(h5Path, matfilePath)
%% This function run through stimuli of a given cell and authorize the user to input the type of syllable for songs
%% Get the environment to figure out on which machine/cluster we are
fprintf(1,'The environment is: %s\n',getenv('HOSTNAME'))

if ~isempty(strfind(getenv('HOSTNAME'),'ln')) || ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))%savio Cluster
    Savio=1;
    Me=0;
    fprintf(1, 'We are on savio!\n')
    addpath(genpath('/global/home/users/jelie/CODE/SingleUnitModels'));
    addpath(genpath('/global/home/users/jelie/CODE/SingleUnitDataMining'));
    addpath(genpath('/global/home/users/jelie/CODE/GeneralCode'));
    addpath(genpath('/global/home/users/jelie/CODE/density_estimation'));
    addpath(genpath('/global/home/users/jelie/CODE/tlab/src'));
    addpath(genpath('/global/home/users/jelie/CODE/strflab/trunk'));
    StimDir = '/global/scratch/jelie/Stims/';
elseif ismac()
    Savio=0;
    Me = 1;
    fprintf(1, 'We are on my Mac Book Pro!\n')
%     addpath(genpath('/Users/elie/Documents/CODE/SingleUnitModels'));
%     addpath(genpath('/Users/elie/Documents/CODE/density_estimation'));
%     addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'));
%     addpath(genpath('/Users/elie/Documents/CODE/STRFLab/trunk'));
%     addpath(genpath('/Users/elie/Documents/CODE/tlab/src'));
%     StimDir = '/Users/elie/Documents/CODE/data/Stims/';
    addpath(genpath('/Users/ellenzippi/Documents/MATLAB/tlab/src'));
    addpath(genpath('/Users/ellenzippi/Documents/MATLAB/SingleUnitDataMining'));
    addpath(genpath('/Users/ellenzippi/Desktop/TheunissenLab/data')); 
    StimDir = '/Users/ellenzippi/Desktop/TheunissenLab/data/Stims/'; 
else %we are on strfinator or a cluster machine
    Savio = 0;
    Me = 0;
    fprintf(1, 'Hum.. We must be on one of the lab machine\n')
    addpath(genpath('/auto/fhome/julie/Code/SingleUnitModels'));
    addpath(genpath('/auto/fhome/julie/Code/density_estimation'));
    addpath(genpath('/auto/fhome/julie/Code/GeneralCode'));
    addpath(genpath('/auto/fhome/julie/Code/strflab'));
    addpath(genpath('/auto/fhome/julie/Code/tlab/src'));
    StimDir = '/auto/tdrive/fdata/julie/Stims';
end



%% Read h5 input data from data base
if Savio %savio Cluster
    Dir_local='/global/scratch/jelie/MatFiles/';
    fprintf('load file %s\n', h5Path);
    unit=loadfromTdrive_savio(h5Path, Dir_local);
elseif Me
    unit = read_unit_h5file(h5Path, 'r');
elseif (~Me && ~Savio)
    unit = read_unit_h5file(h5Path, 'r');
end

%% Read Matfile old output data
Res=load(matfilePath);

%% Select the protocol (SelX, Callx, Maskx, STRFx...)
%number of different protocols run for that unit, identify if the one
%you are asking for is existing in this unit and selecting this one in
%responses
nclasses = length(unit.classes);

classId = 0;
for ic=1:nclasses
    prot=unit.classes{ic};
    stimType='Call';
    [~, Name, ~]=fileparts(unit.source_directory);
    if strcmp('WhiWhi4522M', Name) && strcmp('Site1', unit.site)
        if strcmp('Call1c', prot(1:6))
            classId = ic;
            break;
        end
    elseif strcmp('Call', prot(1:4))
        classId = ic;
        break;
    end
end

if classId ~= 0 
    %% Get responses
    responses = unit.class_responses.(prot);

    %% Iterate through sound files and add input
    % This is the number of sound files played
    nfiles = length(responses);
    
    % set up output cell arrays
    Res.VocTypeHD = cell(nfiles,1); 

    for isound = 1:nfiles
        fprintf('sound %d/%d\n', isound, nfiles)
        response=responses{isound};
        stim_name=response.tdt_wavfile;

        disp(Res.VocType{isound})
        if strcmp(Res.VocType{isound},'song') == 1
            % Read the stim wave files on the cluster on a local mac machine.
            IndS = strfind(stim_name,'Stims');
            stim_name = strcat(StimDir, stim_name((IndS+6):end));
            [sound_in, samprate] = audioread(stim_name);
            info = audioinfo(stim_name); 

            % create on as a vector of zeros same length as sound_in
            on = zeros(size(sound_in)); 

            % identify the indices that correspond to selected sound by
            % replacing zero values by 0.8 using WaveIndices variable
            on = on(Res.WavIndices{isound}(1):Res.WavIndices{isound}(2))+0.8; 

            % play sound
            soundsc(sound_in, samprate, info.BitsPerSample)
%             audioplayer(sound_in, samprate, info.BitsPerSample); 
            plot_sound_spike_selection_h5(response, sound_in, samprate, on);

            % little input line that ask for what type of syllable if it is a
            % song (check VocType)
            Res.VocTypeHD{isound} = input('Introductory (in) or motif (mo)?', 's'); 
        else
            Res.VocTypeHD{isound} = Res.VocType{isound}; 
        end
        VocTypeHD = Res.VocTypeHD; 
    end
%% Saving data
     save(matfilePath, '-struct', 'Res');
    fprintf('saved data under: %s\n', matfilePath);
else
    fprintf(1, 'Warning: could not find stimType %s in h5 file %s\n', stimType, h5Path);
    fprintf(1, '\tAvailable options are:\n');
    for ic=1:nclasses
        fprintf(1,'%s\n', unit.classes{ic});
    end
end
end
    