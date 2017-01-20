function [filename]=Matfile_construct_F1s_JEE(h5Path, OldWav, pl)
%% This function cut tdt stim to isolate the begining of the stimulus which
...might contains several vocalizations in the first second

if nargin<2
    OldWav=[]; %If you don't have access to the name of the wave in the vocalization banck, don't try to look for them
end
if nargin<3
    pl=0; %set to 0 for no graph; 1 if you want to see final results per stim; 2 if you to see both steps graph
end


%%  I choose 1000 ms as a window to cut the voc of each stim
... to values obtain with soundduration_invest.m
Win=1;

%% Read input data from data base
unit = read_unit_h5file(h5Path, 'r');

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
    Trials_GaussFiltered = cell(nfiles,1);
    PSTH=cell(nfiles,1);
    PSTH_GaussFiltered=cell(nfiles,1);
    HwidthSpikes = cell(nfiles,1);
    Rate_BG=nan(nfiles,1);
    Spectro=cell(nfiles,1);
    

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
        if ismac()
            [~, username] = system('who am i');
            if strcmp(strtok(username), 'frederictheunissen')
                if strncmp('/auto/fdata/solveig',stim_name, 19)
                    stim_name = strcat('/Users/frederictheunissen/Documents/Data/solveig', stim_name(20:end));
                elseif strncmp('/auto/fdata/julie',stim_name, 17)
                    stim_name = strcat('/Users/frederictheunissen/Documents/Data/Julie', stim_name(18:end));
                end
            elseif strcmp(strtok(username), 'elie')
                if strncmp('/auto/fdata/solveig',stim_name, 19)
                    stim_name = strcat('/Users/frederictheunissen/Documents/Data/solveig', stim_name(20:end));
                elseif strncmp('/auto/fdata/julie',stim_name, 17)
                    stim_name = strcat('/Users/elie/Documents/CODE/data', stim_name(18:end));
                end
            end
        else
            stim_name = ['/auto/tdrive/' stim_name(strfind(stim_name,'fdata'):end)];
        end
        [sound_in, samprate] = audioread(stim_name);
        if pl>0
            fprintf(1,'Plot the sound isolated and the neural response\n');
            on=ones(size(sound_in));
            if length(sound_in)>Win*samprate
                on(Win*samprate+1 :end)=0;%Only keep the first response
            end
            on=on*0.8;
            plot_sound_spike_selection_h5(response, sound_in, samprate, on);
            pause(1);
        end
        
        % Save the wave file and wave indices
        if length(sound_in)<Win*samprate
            EndIndex_sound = length(sound_in);
        else
            EndIndex_sound=Win*samprate;
        end
        WavIndices{isound} = [1; EndIndex_sound];
        SectionWave{isound} = sound_in(WavIndices{isound}(1): WavIndices{isound}(2));
        SectionLength(isound) = length(SectionWave{isound})/samprate*1000;
        
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
                EndIndex = EndIndex_sound + Time_needed*samprate;
            else
                EndIndex = EndIndex_sound + post_response_time*samprate;
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
        winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
        winLength = fix(winLength/2)*2;            % Enforce even window length
        increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
        %calculate spectro
        [s, to, fo, ~] = GaussianSpectrum(SectionWave{isound}, increment, winLength, samprate);
        
        %reshape and store spectro
        D=length(to);
        F=length(fo);
        Spectro{isound}=reshape(abs(s),1,F*D);
        Spectroto= to;
        Spectrofo= fo;

        %% Isolate spikes that relate to the section and...
        ...calculate average (psth) for this section.
        [Trials{isound},Trials_GaussFiltered{isound},PSTH{isound},PSTH_GaussFiltered{isound},~,~,~,HwidthSpikes{isound}] = spikeTimes_psth_gaussfilter_cal(1, EndIndex, samprate,response,Rate_BG(isound), Win, pl);

        %% Plot sound pressure waveform, spectrogram and psth isolated
        if pl>0
            figure(5)
            [s, to, fo, ~] = GaussianSpectrum([SectionWave{isound} ;zeros(Win*samprate-length(SectionWave{isound}),1)], increment, winLength, samprate);
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
            plot(1:length(PSTH_GaussFiltered{isound}),PSTH_GaussFiltered{isound}*1000,'-r');
            legend('PSTH', 'Gaussian filtered PSTH')
            %axis([v_axis(1) v_axis(2) 0 1000*max(smpsth)]);
            ylabel('Rate (spikes/s)');
            xlabel('Time (ms)');
            hold off

            subplot(3,1,3)
            cla;
            wave=[SectionWave{isound} ;zeros(Win*samprate-length(SectionWave{isound}),1)];
            plot(wave)
            pause(1)
        end
        
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
    
    if pl>0
        figure(6)
        plot(cell2mat(PSTH_GaussFiltered)');
        xlabel('Time (ms)')
        ylabel('spike rate spike/ms')
        title('Gaussian filtered spike rates of all the stimuli')
    end
        
    
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
    Res.Trials_GaussFiltered = Trials_GaussFiltered;% Contains the spike rate per ms in bins of 1ms from the begining of the section and not in ms from the begining of the section
    Res.PSTH=PSTH;
    Res.Spectro=Spectro;
    Res.Spectroto=Spectroto;
    Res.Spectrofo=Spectrofo;
    Res.PSTH_GaussFiltered=PSTH_GaussFiltered;
    Res.HwidthSpikes = HwidthSpikes;

    if ismac()
            [~, username] = system('who am i');
            if strcmp(strtok(username), 'frederictheunissen')
                if strncmp('/auto/fdata/solveig',stim_name, 19)
                elseif strncmp('/auto/fdata/julie',stim_name, 17)
                    filename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['FirstVoc1s_' Res.Site '.mat']);
                end
            elseif strcmp(strtok(username), 'elie')
                filename = fullfile('/Users','elie','Documents','CODE','data','matfile','FirstVoc1sMat',['FirstVoc1s_' Res.Site '.mat']);
            end
    else
        filename=fullfile('/auto','tdrive','julie','k6','julie','matfile',Res.subject,['FirstVoc1s_' Res.Site '.mat']);
    end
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
    