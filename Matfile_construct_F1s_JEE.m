function [filename]=Matfile_construct_F1s_JEE(h5Path, OldWav, pl)
%% This function cut tdt stim to isolate the begining of the stimulus which
...might contains several vocalizations in the first second
    
if nargin<2
    OldWav=[]; %If you don't have access to the name of the wave in the vocalization banck, don't try to look for them
end
if nargin<3
    pl=0; %set to 0 for no graph; 1 if you want to see final results per stim; 2 if you to see both steps graph
end

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
elseif ismac()
    Savio=0;
    Me = 1;
    fprintf(1, 'We are on my Mac Book Pro!\n')
    addpath(genpath('/Users/elie/Documents/CODE/SingleUnitModels'));
    addpath(genpath('/Users/elie/Documents/CODE/density_estimation'));
    addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'));
    addpath(genpath('/Users/elie/Documents/CODE/STRFLab/trunk'));
    addpath(genpath('/Users/elie/Documents/CODE/tlab/src'));
else %we are on strfinator or a cluster machine
    Savio = 0;
    Me = 0;
    fprintf(1, 'Hum.. We must be on one of the lab machine\n')
    addpath(genpath('/auto/fhome/julie/Code/SingleUnitModels'));
    addpath(genpath('/auto/fhome/julie/Code/density_estimation'));
    addpath(genpath('/auto/fhome/julie/Code/GeneralCode'));
    addpath(genpath('/auto/fhome/julie/Code/strflab'));
    addpath(genpath('/auto/fhome/julie/Code/tlab/src'));
end

%%  I choose 1000 ms as a window to cut the voc of each stim
... to values obtain with soundduration_invest.m
    Win=1;
Response_samprate=1000;%I choose to analyse the spike patterns with a 1kHz resolution

%% Read input data from data base
%% Load the unit matfile
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/GreBlu9508M/ZS_Site2_L1100R1450_21.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L1100R1450_e21_s0_ss1.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L2000R1600_e27_s1_ss1.mat')
if Savio %savio Cluster
    Dir_local='/global/scratch/jelie/MatFiles/';
    fprintf('load file %s\n', h5Path);
    unit=loadfromTdrive_savio(h5Path, Dir_local);
elseif Me
    unit = read_unit_h5file(h5Path, 'r');
elseif (~Me && ~Savio)
    unit = read_unit_h5file(h5Path, 'r');
end

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
    %% Get the responses for that unit
    responses = unit.class_responses.(prot);
    
    %% Now run cut the neural responses and stims at Win seconds after the begining of the stim
    % This is the number of sound files played
    nfiles = length(responses);
    
    % Create the structure that will contain the results
    Res = struct();
    
    % indicate info on the site
    [pathh5, nameh5, ~]=fileparts(h5Path);
    Res.subject=pathh5(end-10:end);
    Res.Site=nameh5;
    
    % prepare cell arrays
    VocType=cell(nfiles,1);
    TDT_wavfiles=cell(nfiles,1);
    VocBank_Wavfiles=cell(nfiles,1);
    Original_wavfiles=cell(nfiles,1);
    SectionWave=cell(nfiles,1);
    WavIndices=cell(nfiles,1);
    SectionLength = zeros(nfiles,1);
    ESex=cell(nfiles,1);
    Eage=cell(nfiles,1);
    Erelated=cell(nfiles,1);
    Trials=cell(nfiles,1);
    JackKnife_KDE_Filtered = cell(nfiles,1);
%    JKinput_Trialsfiltered = cell(nfiles,1);
%    Trials_GaussFiltered=cell(nfiles,1);
    PSTH=cell(nfiles,1);
    PSTH_KDE_Filtered=cell(nfiles,1);
    Spike_array = cell(nfiles,1);
%     Kth_Neigh = cell(nfiles,1);
%     Kth_Neigh_JK = cell(nfiles,1);
%     HwidthSpikes_PSTH = cell(nfiles,1);
%     HwidthSpikes_JK = cell(nfiles,1);
    Rate_BG=nan(nfiles,1);
    Spectro=cell(nfiles,1);
    Spectroto = cell(nfiles,1);
    Spectrofo = cell(nfiles,1);
    Ntrials = nan(nfiles,1);
    Samprate = nan(nfiles,1);
    
    %% Configure Parallel computing
    if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
        MyParPool = parpool(str2num(getenv('SLURM_CPUS_ON_NODE')),'IdleTimeout', Inf);
        system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
        [~,JobID] = system('echo $SLURM_JOB_ID');
        parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];
    end
    
    %% Loop through stims
    %parfor
    for isound = 1:nfiles
        fprintf(1,'sound %d/%d\n',isound, nfiles);
        response=responses{isound};
        stim_name=response.tdt_wavfile;
        stim_number=str2double(response.number);
        if ~isempty(OldWav)
            VocBank_idx= find(cell2mat(OldWav(:,1))==stim_number);
        end
        
        % Read and save the callid.
        if strcmp(response.stim_type, 'call')
            VocType{isound}=response.callid;
        else
            VocType{isound}=response.stim_type; %song or ml noise
        end
        
        % Read the stim wave files on the cluster or on a local mac machine.
        if Me
            stim_name = strcat('/Users/elie/Documents/CODE/data', stim_name(18:end));
            [sound_in, Samprate(isound)] = audioread(stim_name);
        elseif Savio
            stim_name = ['/global/scratch/jelie/' stim_name(strfind(stim_name,'Stims'):end)];
            [sound_in, Samprate(isound)] = audioread(stim_name);
        else
            stim_name = ['/auto/tdrive/' stim_name(strfind(stim_name,'fdata'):end)];
            [sound_in, Samprate(isound)] = audioread(stim_name);
        end
        
        
        
        
        % Plot the sound if asked
        if pl>0
            fprintf(1,'Plot the sound isolated and the neural response\n');
            on=ones(size(sound_in));
            if length(sound_in)>Win*Samprate(isound)
                on(Win*Samprate(isound)+1 :end)=0;%Only keep the first response
            end
            on=on*0.8;
            plot_sound_spike_selection_h5(response, sound_in, Samprate(isound), on);
            pause(1);
        end
        
        % Save the wave file and wave indices
        if length(sound_in)<Win*Samprate(isound)
            EndIndex_sound = length(sound_in);
        else
            EndIndex_sound=Win*Samprate(isound);
        end
        WavIndices{isound} = [1; EndIndex_sound];
        SectionWave{isound} = sound_in(WavIndices{isound}(1): WavIndices{isound}(2));
        SectionLength(isound) = length(SectionWave{isound})/Samprate(isound)*1000;
        
        %Determine if there is enough background response after the end of
        %stim to complete the neural response up to Win ms
        Time_needed = Win-SectionLength(isound)/1000;
        if Time_needed>0
            ntrials = length(response.trials);
            post_response_time=nan(1,ntrials);
            for it = 1:ntrials
                post_response_time(it) = response.trials{it}.post_time; %time left after the stim in s
            end
            post_response_time=min(post_response_time);
            if post_response_time>=Time_needed
                EndIndex = EndIndex_sound + Time_needed*Samprate(isound);
            else
                EndIndex = EndIndex_sound + post_response_time*Samprate(isound);
            end
        else
            EndIndex = EndIndex_sound;
        end
        
        
        % Save average background activity calculated in the previous time
        % window before stimulus presentation
        Rates=nan(length(response.trials),1);
        for it =1:length(response.trials)
            Rates(it) = response.trials{it}.bg_rate;
        end
        Rate_BG(isound) = mean(Rates);
        
        %% calculate and store spectro
        % Parameters for the Spectrogram
        nstd = 6;
        fband = 50;
        twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
        winLength = fix(twindow*Samprate(isound)/1000.0);  % Window length in number of points
        winLength = fix(winLength/2)*2;            % Enforce even window length
        increment = fix(0.001*Samprate(isound));           % Sampling rate of spectrogram in number of points - set at 1 kHz
        %calculate spectro
        [s, to, fo, ~] = GaussianSpectrum(SectionWave{isound}, increment, winLength, Samprate(isound));
        
        %reshape and store spectro
        D=length(to);
        F=length(fo);
        Spectro{isound}=reshape(abs(s),1,F*D);
        Spectroto{isound}= to;
        Spectrofo{isound}= fo;
        
        %% Isolate spikes that relate to the section and...
        ...calculate average (psth) for this section.
        Ntrials(isound) = length(response.trials);
% The following code was investigating the effect of various degree of nearest neighbor     
%        Kth_Neigh{isound} =  Ntrials(isound) * flip((1:Ntrials(isound)).^-1); 
%        Kth_Neigh_JK{isound} =  (Ntrials(isound)-1) * flip((1:(Ntrials(isound)-1)).^-1); 
        
%         Kth_Neigh{isound} =  Ntrials(isound)/3;
%         Kth_Neigh_JK{isound} = (Ntrials(isound)-1)/3;
%        [Trials{isound},Trials_GaussFiltered{isound},PSTH{isound},PSTH_GaussFiltered{isound},JackKnife_GaussFiltered{isound},JKinput_Trialsfiltered{isound},~,~,~,HwidthSpikes{isound}] = spikeTimes_psth_gaussfilter_cal(1, EndIndex, Samprate(isound),response,Rate_BG(isound), Win,Response_samprate,Kth_Neigh, pl);
%        [Trials{isound},Trials_GaussFiltered{isound},PSTH{isound},PSTH_GaussFiltered{isound},JackKnife_GaussFiltered{isound},~,~,~,HwidthSpikes{isound}] = spikeTimes_psth_gaussfilter_cal(1, EndIndex, Samprate(isound),response,Rate_BG(isound), Win,Response_samprate,Kth_Neigh{isound}, Kth_Neigh_JK{isound}, pl);
%        [Trials{isound},PSTH{isound},PSTH_GaussFiltered{isound},JackKnife_GaussFiltered{isound},HwidthSpikes_PSTH{isound},HwidthSpikes_JK{isound},~,~,~] = spikeTimes_psth_gaussfilter_cal(1, EndIndex, Samprate(isound),response,Rate_BG(isound), Win,Response_samprate,Kth_Neigh{isound}, Kth_Neigh_JK{isound}, pl);
        [Trials{isound},PSTH{isound},PSTH_KDE_Filtered{isound},JackKnife_KDE_Filtered{isound},Spike_array{isound},~,~,~] = spikeTimes_psth_KDE_filter_cal(1, EndIndex, Samprate(isound),response,Rate_BG(isound), Win,Response_samprate, pl);
        
        %% Store other infos on section
        TDT_wavfiles{isound}=response.tdt_wavfile;
        Original_wavfiles{isound}=response.original_wavfile;
        if strcmp(response.stim_type,'call')
            ESex{isound}=response.stim_source_sex;
            Eage{isound}=response.callerAge;
            Erelated{isound}=response.stim_source;
            if ~isempty(OldWav)
                VocBank_Wavfiles{isound}=OldWav{VocBank_idx, 3};
            end
        else
            ESex{isound}='NaN';
            Eage{isound}='NaN';
            Erelated{isound}='NaN';
            if ~isempty(OldWav)
                VocBank_Wavfiles{isound}=OldWav{VocBank_idx, 3};
            end
        end
    
    end
    
    %% Plots if requested
    % Plot sound pressure waveform, spectrogram and psth isolated
        
    if pl>0
        for isound=1:nfiles
            figure(5)
            % Parameters for the Spectrogram
        nstd = 6;
        fband = 50;
        twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
        winLength = fix(twindow*Samprate(isound)/1000.0);  % Window length in number of points
        winLength = fix(winLength/2)*2;            % Enforce even window length
        increment = fix(0.001*Samprate(isound));           % Sampling rate of spectrogram in number of points - set at 1 kHz
            [s, to, fo, ~] = GaussianSpectrum([SectionWave{isound} ;zeros(Win*Samprate(isound)-length(SectionWave{isound}),1)], increment, winLength, Samprate(isound));
            subplot(3,1,1);
            DBNOISE = 40;
            logB = 20*log10(abs(s));
            maxB = max(max(logB));
            minB = maxB-DBNOISE;
            imagesc(to,fo,logB);          % to is in seconds
            axis xy;
            caxis('manual');
            caxis([minB maxB]);
            v_axis(1) = 0;
            v_axis(2) = Win;
            v_axis(3)= 0;
            v_axis(4)= 12000;
            axis(v_axis);
            ylabel('Frequency kHz');
            xlabel('Time (s)');
            cmap = spec_cmap();
            colormap(cmap);

            subplot(3,1,2);
            cla;
            %             wind1 = hanning(31)/sum(hanning(31));   % 31 ms smoothing
            %             smpsth = conv(psth,wind1);
            plot(1:length(PSTH{isound}),PSTH{isound}*1000, '-b');
            hold on
            plot(1:length(PSTH_KDE_Filtered{isound}),PSTH_KDE_Filtered{isound}*1000,'-r');
            legend('PSTH', 'Gaussian filtered PSTH')
            %axis([v_axis(1) v_axis(2) 0 1000*max(smpsth)]);
            ylabel('Rate (spikes/s)');
            xlabel('Time (ms)');
            hold off

            subplot(3,1,3)
            cla;
            wave=[SectionWave{isound} ;zeros(Win*Samprate(isound)-length(SectionWave{isound}),1)];
            plot(wave)
            pause(1)
        end
    end
    
    %% Get ready indices of unique sets of JK PSTH for subsequent use in Semantic_NeuroInfo_Poisson_savio.m for stims of interests
    Nb_bootstrap = 100;
    UCat = {'Ag', 'Be','DC','Di','LT', 'Ne','Te','Th','song'};
    IndVoc = zeros(nfiles,1);
    for cc=1:length(UCat)
        IndVoc = IndVoc + strcmp(VocType, UCat{cc});
    end
    GoodVocInd = find(IndVoc);
    NbGoodVoc = length(GoodVocInd);
    MinNbTrials = min(Ntrials(GoodVocInd));
    SetIndices_JK = cell(Nb_bootstrap,1);
    parfor bb=1:Nb_bootstrap
        LocalSet = nan(MinNbTrials, NbGoodVoc);
        for vv=1:NbGoodVoc
            LocalSet(:,vv) = randperm(Ntrials(GoodVocInd(vv)),MinNbTrials);
        end
        SetIndices_JK{bb} = LocalSet;
    end
    
    if pl>0
        figure(6)
        plot(cell2mat(PSTH_KDE_Filtered)');
        xlabel('Time (ms)')
        ylabel('spike rate spike/ms')
        title('KDE filtered spike rates of all the stimuli')
    end
    
    %% Fill-in output structure
    Res.VocType=VocType; % this is the type of vocalization (e.g. distance call DC, Nest call Ne, Aggressive call Ag...)
    Res.VocBank_wavfiles=VocBank_Wavfiles; % name of the wav file of the vocalization bank to which this section responded
    Res.Original_wavfiles=Original_wavfiles; % The real stim is a combination of 1 or 3 calls or 2.5s song. This is the original name of the wav file JEE constructed with the vocalization from the vocalization bank.
    Res.TDT_wavfiles=TDT_wavfiles; % The name of same previous wav stim given by TDT (stim1, stim2.... stim136...)
    Res.WavIndices = WavIndices; % begining and end indices of each section within the TDT_wavfiles, just to be able to retrieve the wavform if needed
    Res.SectionWave=SectionWave;
    Res.SectionLength = SectionLength; % duration of each section in ms
    Res.ESex=ESex; % Sex of the emitter of the vocalization
    Res.Eage=Eage;  % Age of the emitter of the vocalization
    Res.Erelated=Erelated; % Relation of the emitter to the subject (familiar, unfamiliar, self)
    Res.Trials=Trials; % Contains the spike arrival times in ms from the begining of the section and not in ms from the begining of the stim as in h5 files!!!
    Res.JackKnife_KDE_Filtered = JackKnife_KDE_Filtered;% Contains the spike rate of the ntrials JackKnifes sampled at Response_samprate from the begining of the section
%    Res.JKinput_Trialsfiltered = JKinput_Trialsfiltered;% Contains the spike rate of the ntrials convolved with time varying width estimated with the Jackknife estimated sampled at Response_samprate from the begining of the section
%    Res.Trials_GaussFiltered = Trials_GaussFiltered;
    Res.Response_samprate = Response_samprate; % Sampling rate of the neural responses in Hz
    Res.PSTH=PSTH;
    Res.Spectro=Spectro;
    Res.Spectroto=Spectroto;
    Res.Spectrofo=Spectrofo;
    Res.PSTH_KDE_Filtered=PSTH_KDE_Filtered;% Contains the spike rate (Gaussian filtered) calculated with all trials and sampled at Response_samprate from the begining of the section
%    Res.Kth_Neigh = Kth_Neigh;
%    Res.Kth_Neigh_JK = Kth_Neigh_JK;
    Res.Spike_array = Spike_array;
%    Res.HwidthSpikes_PSTH = HwidthSpikes_PSTH;
%    Res.HwidthSpikes_JK = HwidthSpikes_JK;
    Res.SetIndices_JK = SetIndices_JK;
    
    %% Get ready saving files and directories
    %OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatSavio');
    
    if Savio
        OutputDir='/global/scratch/jelie/MatFiles/FirstVoc1sMat';
    elseif Me
        OutputDir='/users/elie/Documents/CODE/data/matfile/FirstVoc1sMat';
    else
        OutputDir=fullfile('/auto','tdrive','julie','k6','julie','matfile','FirstVoc1sMat');
    end
    filename=fullfile(OutputDir,['FirstVoc1s_' Res.Site '.mat']);
    
    save(filename, '-struct', 'Res');
    fprintf('saved data under: %s\n', filename);
    clear duration Res VocType TDT_wavfiles Cut_orders Original_wavfiles WavIndices SectionWave ESex Eage Erelated Trials PSTH MeanRate StdRate Spectro Section_cat Section_zscore SectionLength Section_tvalue Section_pvalue sections_good_zscores
else
    fprintf(1, 'Warning: could not find stimType %s in h5 file %s\n', stimType, h5Path);
    fprintf(1, '\tAvailable options are:\n');
    for ic=1:nclasses
        fprintf(1,'%s\n', unit.classes{ic});
    end
    clear duration unit pthreshold
end
end
