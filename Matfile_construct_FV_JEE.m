function [filename]=Matfile_construct_FV_JEE(h5Path,  pl)
%% This function cut tdt stim to isolate vocalizations as much as possible...
...and just take the first vocalization of each stim that it places in ...
    ...a window of xxxms

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
    StimDir = '/global/scratch/jelie/Stims';
elseif ismac()
    Savio=0;
    Me = 1;
    fprintf(1, 'We are on my Mac Book Pro!\n')
    addpath(genpath('/Users/elie/Documents/CODE/SingleUnitModels'));
    addpath(genpath('/Users/elie/Documents/CODE/density_estimation'));
    addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'));
    addpath(genpath('/Users/elie/Documents/CODE/STRFLab/trunk'));
    addpath(genpath('/Users/elie/Documents/CODE/tlab/src'));
    StimDir = '/Users/elie/Documents/CODE/data/Stims/';
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

%%  I choose 600ms as a window to put the first voc of each stim in according...
... to values obtain with soundduration_invest.m

duration=0.02; % duration in ms of min silence period
if nargin<2
    pl=0; %set to 0 for no graph; 1 if you want to see final results per stim; 2 if you to see both steps graph
end

Win=0.6;

%% Parameters for the estiamtion of spike rate
%1 000 corresponds to 1 ms precision, TDT sampling rate is little above 24kHz so that's good enough 
Response_samprate = 1000;
% Number of time bins of the response with choosen sampling rate. For 1bin=1ms response_samprate=1000
Ntimebins = round(Response_samprate*Win);
% Setting the vector of points in ms at which an output value will be given
% for the filtered spike patterns
Tin = (1:Ntimebins)-0.5;

%% Read input data from data base
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
    %% Retrieve name of the stimuli in the vocalization database
    [OldWav]=relate_sound_files(Name, StimDir);
    
    %% Get responses
    responses = unit.class_responses.(prot);

    %% Cutting procedure based on zero values to calculate intensity values of each cut
    % This is the number of sound files played
    nfiles = length(responses);
    
    % set up output cell arrays
    intensity = cell(nfiles,1);
    dur = cell(nfiles,1);
    voctype = cell(nfiles,1);
    Sound = cell(nfiles,1);
    histDur=[];
    
    % Find on which machine we are and where are stored the data
    [SoundPath] =fileparts(responses{1}.original_wavfile);
    SoundPathParts = textscan(SoundPath, '%s', 'Delimiter', '\\/');
    DataDir = '';
    k = 2;
    while ~strcmp(SoundPathParts{1}{k},'Stims')
        DataDir = strcat(DataDir, '/', SoundPathParts{1}{k});
        k=k+1;
    end
    if ismac()
            [~, username] = system('who am i');
            if strcmp(strtok(username), 'frederictheunissen')
                if strcmp('/auto/fdata/julie',DataDir)
                    DataDir='/Users/frederictheunissen/Documents/Data/Julie';
                end
            elseif strcmp(strtok(username), 'elie')
                if strcmp('/auto/fdata/julie',DataDir)
                    DataDir='/Users/elie/Documents/CODE/data';
                    
                end
            end
    end
    
       %% Configure Parallel computing
    if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
        MyParPool = parpool(str2num(getenv('SLURM_CPUS_ON_NODE')),'IdleTimeout', Inf);
        %MyParPool = parpool(20,'IdleTimeout', Inf);
        system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
        [~,JobID] = system('echo $SLURM_JOB_ID');
        parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];
    end
    
    %Loop through files, detect sound periods
   parfor isound = 1:nfiles
        response=responses{isound};
        stim_name1=response.tdt_wavfile;
        stim_name1 = strcat(DataDir, stim_name1(18:end));
        [sound_in, samprate1] = audioread(stim_name1);
        iintensity=[];
        idur=[];
        iSound = zeros(0);
    
       % Read the callid.
        if strcmp(response.stim_type, 'song')
            vocid=response.stim_type;
        elseif strcmp(response.stim_type, 'call')
            vocid=response.callid;
        else
            vocid=response.stim_type;
        end
    
        %if strcmp(callId,'Te')||strcmp(callId,'Th')||strcmp(callId,'DC')||strcmp(callId,'LT')
            % Divide the sound into areas of non-zero time
            ind = find(sound_in);
            on = zeros(size(sound_in));
            sound_p1 = zeros(size(sound_in));
            sound_p1(1:end-1) = sound_in(2:end);
            sound_m1 = zeros(size(sound_in));
            sound_m1(2:end) = sound_in(1:end-1);
            sound_diff = sound_p1-sound_m1;
            ind2 = find(sound_diff);
            
            % Sounds is "on" when value is not zero and derivative is not zero
            on(ind) = 0.08;    % It is set at 0.8 to make a nice graphic.
            on(ind2) = 0.08;
            if pl==2        
                plot_sound_spike_selection_h5(response, sound_in, samprate1, on);
                pause(1);
            end
    
            % Find beggining and end
            laston = 0;
            for i=1:length(on)
                if (laston == 0) && (on(i)~=0)
                    begInd = i;
                elseif (laston~=0) && (on(i) == 0) || (laston~=0) && (i==length(on))
                    endInd1 = i;
                    durInd = endInd1 - begInd;
                
                    
                    if (durInd > samprate1*duration)    % Only examine stimulus intervals of certain length
                    % calculate duration, peak frequency and intensity
                    durVal = (endInd1-begInd)/samprate1;
                    intVal = std(sound_in(begInd:endInd1));
                    %[Pxx,f] = pwelch(sound_in(begInd:endInd), nw, fix(nw/2) , nw, samprate);
                    %Pmax = max(Pxx);
                    %fVal = f(Pxx==Pmax);
            
                     % Stuff values
                    iintensity = [iintensity intVal];
                    histDur = [histDur durVal];
                    idur = [idur durVal];
                    iSound = [iSound [begInd ; endInd1]];
                    % fprintf(1, 'Beg= %.2f End= %.2f Int= %.2f Dur= %.2f Freq= %.0f RD= %.2f\n', begTime, endTime, 20*log10(intVal), durVal, fVal, rateDiffVal);
                    end
                end
                laston = on(i);
            end
            voctype{isound} = vocid;       
            intensity{isound}=iintensity;
            dur{isound}=idur;
            Sound{isound}=iSound;
    end
    if pl
        figure(2)
        hist(histDur,60);
    end



    %%  Define the length of the section
    response=responses{1};
    stim_name=response.tdt_wavfile;
    stim_name = strcat(DataDir, stim_name(18:end));
    [~, samprate] = audioread(stim_name);
    Lsection=ceil(Win*samprate);



    %% Now run a second cutting system based on the intensity of the signal to get nicer cuts and enlarged the saved window around the vocalization
    % Create the structure that will contain the results
    Res = struct();
    %'subject',{},'Site', {},'VocType', {}, 'tdt_wavfiles', {},'original_wavfiles',{}, 'EndBegIndices', {}, 'Cut_orders', {}, 'sex', {}, 'age', {}, 'related', {}, 'trials', {}, 'PSTH',{}, 'spectro', {}, 'section_cat', {});

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
    Trials_BG=cell(nfiles,1);
    PSTH=cell(nfiles,1);
    KDE_Rate=cell(nfiles,1);
    PSTH_BG=cell(nfiles,1);
    Spectro=cell(nfiles,1);
    Section_cat=cell(nfiles,1);
    Spectroto = cell(nfiles,1);
    Spectrofo = cell(nfiles,1);
    
  

    for isound = 1:nfiles
        fprintf('sound %d/%d\n', isound, nfiles)
        response=responses{isound};
        stim_name=response.tdt_wavfile;
        stim_number=str2double(response.number);
        if ~isempty(OldWav)
            VocBank_idx= find(cell2mat(OldWav(:,1))==stim_number);
        end
    
        % Read the stim wave files on the cluster on a local mac machine.
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
    
        %Find Silences longer than 60ms between vocalizations
        iintensity=intensity{isound};
        %%%%%%%%%%%%%%%Increase this a little more??? to enable separation
        %%%%%%%%%%%%%%%of nest calls?
        Thresh=0.35*mean(iintensity); %silence is defined as 35% of average intensity calculated above
        yy = 1;
        Duration = ceil(0.01*samprate);% calls have to be spaced by at least 10ms to be separated. 10ms tail is added after the end sensus stricto of each sound and before the begining
        Sil = zeros(0); %Initialize the matrix that will contain start and end indices of silence slots
        while (yy>0) && (yy<length(sound_in))
            if abs(sound_in(yy))<Thresh
                Beg = yy;
                while abs(sound_in(yy))<Thresh && (yy<length(sound_in))
                    yy = yy + 1;
                end
                End = yy;
                if End-Beg>(Duration)
                    Sil = [Sil [Beg ; End]];
                end
            end
            yy = yy + 1;
        end
        [~,c]=size(Sil);
        on=ones(size(sound_in));
        if c~=0
            for cc=1:c
                on((Sil(1,cc)+Duration):1:(Sil(2,cc)-Duration))=0;% indicate zones of silences (0) and zones around calls that should be kept (60ms after the end and 10ms before the begining)
            end
            if Sil(2,c)==length(on)%if the wavfile ends by a silence
                on((Sil(2,c)-Duration):1:end)=0;% then make sure the silence was not supress at the end by the previous statement
            end
        end
        if on(end)==1
            on= [on;ones(Duration,1)];%in case the stim end by a vocalization, add some space after (10ms)
        end
        on=on*0.8;
        if pl>0        
            plot_sound_spike_selection_h5(response, sound_in, samprate, on);
            pause;
        end
    
        %% Isolate spikes that relate to before the 1rst vocalization and calculate avg background spike rate
        
        ntrials=length(response.trials);
        bgdur=zeros(1,ntrials);
        spike_count_bg = zeros(1,ntrials);
        ntimebins = round(1000*Win); % Number of time bins in ms of the background psth window
        psthbg = zeros(1, ntimebins);
        BG_Trials_temp = cell(ntrials,1);
        for it=1:ntrials
            trial = response.trials{it};
            bgdur(it)=trial.pre_time;
            spike_time_trials = trial.spikeTimes;
            
            %calculate spike rate for he whoe background period
            spike_time1=(spike_time_trials<0);
            spike_time2=(spike_time_trials>=(-bgdur(it)));
            SpikesIndices=find(spike_time1 .* spike_time2);
            spike_count_bg(it)=length(SpikesIndices);
            %bg_mean_rate = mean(spike_count_bg)./(bgdur);
            %bg_std_rate = std(spike_count_bg)./(bgdur);
            spike_rate_bg = spike_count_bg./bgdur;
            
            %calculate a psth background for a 600ms window (Win) before
            %the stim
            spike_time3=(spike_time_trials>=(-Win));
            SpikesIndicesPSTHbg=find(spike_time1 .* spike_time3);
            Spike_time_PSTHbg=spike_time_trials(SpikesIndicesPSTHbg);
            BG_Trials_temp {it} = (Spike_time_PSTHbg - (-Win))*1000; % here spike times arrival are referred to the begining of the 600ms window instead of the begining of the stim and in ms instead of s
            ns=length(Spike_time_PSTHbg);
            spike_array = zeros(1, ntimebins);
            for is=1:ns
                time_ind = ceil(Spike_time_PSTHbg(is)*1000-(-Win*1000));
                if (time_ind < 1 || time_ind > ntimebins)
                    fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d ntimebins = %d\n', response.number, it, time_ind, ntimebins);
                    continue;
                end
                spike_array(time_ind) = spike_array(time_ind) +1;
            end
            psthbg = psthbg + spike_array;
        end
        psthbg = psthbg./ntrials;
        
        
    
        %% Find beggining and end of two first vocalization-section
        WavIndices_temp=cell(2,1);
        nsections = 0;
        laston = 0;
        begInd=0;
        i=0;
        while (nsections<2) && (i<length(on))
            i = i + 1;
            %for i=1:length(on)
                if (laston == 0) && (on(i)~=0)
                    begInd = i;
                    
                elseif (laston~=0) && (on(i) == 0) || (laston~=0) && (i==length(on))
                    endInd = i;
                    durInd = endInd - begInd;
                    if begInd==1 %in case the stim start by a vocalization, the 10ms before begining of sound has not been added here yet and should be done below
                        durInd=durInd + Duration;
                    end

                    % calculate duration in sec and save the begining and end
                    % indices of the section if longer than 50ms (30ms of
                    % sound)
                    dur{isound}=(durInd)/samprate;
                    if durInd>ceil(0.05*samprate)%section needs to be longer than 50ms (20ms of silence + 30ms of sound) 
                        nsections=nsections+1;%increment the section number for that sound
                        if endInd>length(sound_in) % the stim end by a vocalization ensure that we keep the tail of psth by extending the wavsection
                            WavIndices_temp{nsections}=[begInd length(sound_in)];
                        else
                            WavIndices_temp{nsections}=[begInd endInd];
                        end
                    end
                end
                laston=on(i);
            %end
            
        end
        
        %% Isolate first vocalization section within a window of size Lsection
        PSTH_BG{isound}=psthbg; %save the psth of the correspoonding background calculated earlier
        Trials_BG{isound}=BG_Trials_temp; %save the spike arrival times of the background window on which is calculated the previous psth
        
        % Find how much silence there is between the first and second
        % vocalization sections or the end of the stim
        if nsections==2
            Sil = WavIndices_temp{nsections}(1) - WavIndices_temp{nsections-1}(2);
        elseif i==length(on)%section ends at the end of the stim
            Sil = 0; %no wave without sound available after the sectio
        else
            fprintf(1,'pb with the loop ln 341\n');
        end
        
        %define the begining and end of the section and making sure that...
        ...if the stim start by a vocalization there is still the 10ms of...
            ...silence before 
        if nsections==1 %une seule section découpée
            begInd= WavIndices_temp{nsections}(1);
            endInd= WavIndices_temp{nsections}(2);
        elseif nsections==2 %2 sections découpées
            begInd=WavIndices_temp{nsections-1}(1);
            endInd=WavIndices_temp{nsections-1}(2);
        else
            fprintf('problem with section cutting, the code cut %d when he should cut only the first 2\n', nsections);
        end
            
        if begInd==1 %the stim start by a vocalization
            prezeros=zeros(Duration,1);%padd with zero before for 10ms because was not done earlier
            if pl>0
                fprintf('the stim start by a vocalization, index of stim begining: %d\n', begInd')
            end
            begInd=-(length(prezeros)-1);%-1 because there is the zero index!
        end
            
        durInd =  endInd - begInd;
        
        
        %depending on the extracts length either enlarge the sound extract or cut
        if durInd<=Lsection
            % Extend wav file as much as possible around the vocalization
            % to include surrounding silence periods
            Needzeros = Lsection-durInd;
            if Needzeros<=Sil
                endInd_spike=endInd + Needzeros;
                begInd_spike=begInd;
            elseif Sil==0
                endInd_spike=endInd + ceil(Needzeros/2);
                begInd_spike=begInd - floor(Needzeros/2);
            else
                endInd_spike=endInd + Sil;
                begInd_spike=begInd - (Needzeros-Sil);
            end
            % write WavIndices with the new limits and save the section
            % of sound completed with needed zeros before the stim
            if begInd_spike<=0 %window begins before the start of the stim
                WavIndices{isound}(1)=1;
                prezeros=zeros(abs(begInd_spike),1);
            else
                WavIndices{isound}(1)=begInd_spike;
                prezeros=[];
            end
            if Sil==0 %window stops after the end of the stim
                WavIndices{isound}(2)=length(sound_in);
                postzeros=zeros((endInd_spike - length(sound_in)),1);
            else
                WavIndices{isound}(2)=endInd_spike;
                postzeros=[];
            end
             sec=[prezeros; sound_in(WavIndices{isound}(1):WavIndices{isound}(2)); postzeros];
            
            
            if length(sec)~=Lsection
                fprintf(1,'pb withduration of section, the sound length is %d when it should be %d\n', length(sec), Lsection);
            end
            
            SectionWave{isound} = sec;
            SectionLength(isound) = length(sec)/samprate*1000;%durée de l'extrait en ms ici
            Section_cat{isound}='full';
        else
            %cut wav file and store
            begInd_spike=begInd;
            endInd_spike=begInd + Lsection;
            if begInd_spike<=0 %window begins before the start of the stim
                WavIndices{isound}=[1 endInd_spike];
                prezeros=zeros(abs(begInd_spike),1);
                sec=[prezeros; sound_in(1:endInd_spike)];
            else
                WavIndices{isound}=[begInd_spike endInd_spike];
                sec=sound_in(begInd_spike : endInd_spike);
            end
            SectionWave{isound} = sec;
            SectionLength(isound) = length(sec)/samprate*1000;%durée de l'extrait en ms ici
            Section_cat{isound}='cut';
        end
                    
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
        VectorS=reshape(abs(s),1,F*D);
        Spectroto{isound}= to;
        Spectrofo{isound}= fo;
        Spectro{isound}=VectorS;

        %% Isolate spikes that relate to the section and...
        ...calculate average (psth) for this section.
        [Section_Trials,psth,~, ~, ~, ~, ~] = spikeTimes_psth_cal(begInd_spike, endInd_spike, samprate,response,spike_rate_bg, Win);
        Trials{isound}=Section_Trials;
        PSTH{isound}=psth;
        
        % calculate Kernel density estimate of the rate
        Ntrials=length(Trials{isound});
        Spike_count = nan(Ntrials,1);
        Spikes_Times = [];
       for tt=1:Ntrials
           Spikes_Times = [Spikes_Times Trials{isound}{tt}];
           Spike_count(tt) = length(Trials{isound}{tt});
       end
       
        %Calculate the kernel  density estimate of the rate
        [y,t,~,~,~,~,~] = ssvkernel(Spikes_Times,Tin);
        % check that the input Tin was correctly used
        if sum(Tin == t)~=length(Tin)
            error('WARNING: the kernel density estimation is using a different set of observation time points to return the filtered spike pattern\nThe # of time points used by ssvkernel is %d when we are asking for %d.\n', length(t),length(Tin));
        end
        
        % y is a density function that sums to 1
        % multiplying by the total number of spikes gives the number of expecting spike per time bin (here 1 ms) for all 10 Trials
        % dividing by the number of trials give the expected number of spikes per time bin for a trial
        % multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
        KDE_Rate{isound}=y* sum(Spike_count) / Ntrials * Response_samprate/1000;

        %% Plot sound pressure waveform, spectrogram and psth isolated
        if pl>0
            figure(3)
            
            subplot(3,1,1);
            DBNOISE = 40;
            logB = 20*log10(abs(s));
            maxB = max(max(logB));
            minB = maxB-DBNOISE;            
            imagesc(to,fo,logB);          % to is in seconds
            axis xy;
            caxis('manual');
            caxis([minB maxB]); 
            v_axis = nan(4,1);
            v_axis(1) = 0;
            v_axis(2) = SectionLength(isound)/1000;
            v_axis(3)= 0; 
            v_axis(4)= 12000;
            axis(v_axis);                                
            ylabel('Frequency kHz');
            cmap = spec_cmap();
            colormap(cmap);

            subplot(3,1,2);
            cla;
            wind1 = hanning(31)/sum(hanning(31));   % 31 ms smoothing
            smpsth = conv(psth,wind1);
            plot(((1:length(psth))/1000),smpsth(16:length(smpsth)-15)*1000);
            axis([v_axis(1) v_axis(2) 0 1000*max(smpsth)]);
            ylabel('Rate (spikes/s)');
            xlabel('Time (s)');

            subplot(3,1,3)
            cla;
            wave=SectionWave{isound};
            plot(wave)
            pause
        end
        
        %% Store other infos on section
        VocType{isound}=voctype{isound};
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
    

    Res.VocType=VocType(1:nfiles); % this is the type of vocalization (e.g. distance call DC, Nest call Ne, Aggressive call Ag...)
    Res.Section_cat=Section_cat(1:nfiles); % identify whether the section contain an entire vocalization (full) a portion of a vocalization (cut) just the end of a vocalization (cut_end)
    Res.VocBank_wavfiles=VocBank_Wavfiles(1:nfiles); % name of the wav file of the vocalization bank to which this section responded
    Res.Original_wavfiles=Original_wavfiles(1:nfiles); % The real stim is a combination of 1 or 3 calls or 2.5s song. This is the original name of the wav file JEE constructed with the vocalization from the vocalization bank.
    Res.TDT_wavfiles=TDT_wavfiles(1:nfiles); % The name of same previous wav stim given by TDT (stim1, stim2.... stim136...)
    Res.WavIndices = WavIndices(1:nfiles); % begining and end indices of each section within the TDT_wavfiles, just to be able to retrieve the wavform if needed
    Res.SectionWave=SectionWave(1:nfiles);
    Res.SectionLength = SectionLength(1:nfiles); % duration of each section in ms
    Res.ESex=ESex(1:nfiles); % Sex of the emitter of the vocalization
    Res.Eage=Eage(1:nfiles);  % Age of the emitter of the vocalization
    Res.Erelated=Erelated(1:nfiles); % Relation of the emitter to the subject (familiar, unfamiliar, self)
    Res.Trials=Trials(1:nfiles); % Contains the spike arrival times in ms from the begining of the section and not in ms from the begining of the stim as in h5 files!!!
    Res.PSTH=PSTH(1:nfiles);
    Res.KDE_Rate=KDE_Rate(1:nfiles); % Kernel density estimate of the rate for each section
    Res.Spectro=Spectro(1:nfiles);
    Res.Spectroto=Spectroto;
    Res.Spectrofo=Spectrofo;
    Res.VocDuration=dur;
    Res.PSTH_BG=PSTH_BG(1:nfiles);
    Res.Trials_BG=Trials_BG(1:nfiles);
    
%% Saving data
     if Savio
        OutputDir='/global/scratch/jelie/MatFiles/FirstVocMat';
        system(sprintf('mv /global/home/users/jelie/JobExecutableFiles/JobToDoSavio/ExJob%s* /global/home/users/jelie/JobExecutableFiles/JobToDoSavio/JobDoneSavio/', Res.Site))
    elseif Me
        OutputDir='/users/elie/Documents/CODE/data/matfile/FirstVocMat';
    else
        OutputDir=fullfile('/auto','tdrive','julie','k6','julie','matfile','FirstVocMat');
    end
    filename=fullfile(OutputDir,['FirstVoc_' Res.Site '.mat']);
    
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
    