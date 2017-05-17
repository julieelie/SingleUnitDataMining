function [Spike_ArrivalT_Trials,Psth, Psth_KDEfiltered, JK_KDEfiltered, Spike_array, Stim_mean_rate, Spike_std_rate, Spike_rate_stim] = spikeTimes_psth_KDE_filter_cal(BegInd, EndInd, Stim_samprate,Response,Spike_rate_bg, Win,Response_samprate, Fig)
%% This function calculates the PSTH, isolate arrival times of spikes...
...in reference to the begining of the section (begInd).
    % here spike arrival times in Spike_ArrivalT_Trials are given in ms
    
%% Treating input parameters
if nargin<8
    Fig=0;
end

if nargin<7
    Response_samprate = 1000; %1 000 corresponds to 1 ms precision, TDT sampling rate is little above 24kHz so that's good enough 
end
if Response_samprate~=1000
    error('This code is written for a sapling rate of the response at 1000Hz, otherwise the Kernel density estimation is incorrect\n')
end

%% Extracting important variables
% Expected size of the response vector given de response sampling rate:
% here the output Psth_KDEfiltered_all is completed with avg background rate if needed to match a size of Win s.
Resdur = Win*Response_samprate; 

% Duration of the stimulus in seconds
Stimdur=(EndInd+1-BegInd)/Stim_samprate;

% Number of time bins of the response with choosen sampling rate. For 1bin=1ms response_samprate=1000
Ntimebins = round(Response_samprate*Stimdur); 

% Number of stimulus presentations
Ntrials=length(Response.trials);

%% Isolate spike arrival times, construct spike array, count the number of spikes per trial
Spike_ArrivalT_Trials=cell(Ntrials,1);
Spike_count = zeros(1,Ntrials);
Spike_array = zeros(Ntrials, Ntimebins);

for it=1:Ntrials
    Trial = Response.trials{it};
    spike_time1=(Trial.spikeTimes>=(BegInd/Stim_samprate));
    spike_time2=(Trial.spikeTimes<=(EndInd/Stim_samprate));
    SpikesIndices=find(spike_time1 .* spike_time2);
    Spike_count(it)=length(SpikesIndices);
    Spike_time_trials_section=Trial.spikeTimes(SpikesIndices);
    % Save spike arrival times, referred to the begining of the section instead of the begining of the stim, and in ms instead of s 
    Spike_ArrivalT_Trials{it}=(Spike_time_trials_section+ (1-BegInd)/Stim_samprate).*1000; 
    for is=1:Spike_count(it)
        time_ind = floor(Spike_ArrivalT_Trials{it}(is)*Response_samprate/1000);
        if (time_ind < 1 || time_ind > Ntimebins)
            fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d ntimebins = %d\n', Response.number, it, time_ind, Ntimebins);
            continue;
        end
        Spike_array(it,time_ind) = Spike_array(it,time_ind) +1;
    end
end
Psth = sum(Spike_array,1)./Ntrials;

%% Obtain a concatenated version of spike arrival times for full trial and JK PSTHs
Spike_ArrivalT_AllT = nan(sum(Spike_count),1);
Spike_ArrivalT_JK = cell(Ntrials,1);
ss=0; % spike count for the All-trial set of spikes
for it=1:Ntrials
    % All-trial PSTH
    Spike_ArrivalT_AllT((ss+1):(ss+Spike_count(it))) = Spike_ArrivalT_Trials{it};
    ss = ss + Spike_count(it);
    
    % JK PSTH
    Local_trials = 1:Ntrials;
    Local_trials(it) = [];
    Spike_ArrivalT_JK{it} = nan(sum(Spike_count(Local_trials)),1);
    ss2=0; % spike count for that JK(it) set of spikes
    for it2=1:(Ntrials-1)
        it3 = Local_trials(it2);
        Spike_ArrivalT_JK{it}((ss2+1):(ss2+Spike_count(it3))) = Spike_ArrivalT_Trials{it3};
        ss2 = ss2 + Spike_count(it3);
    end
end

%% Calculate Kernel density estimation of the spike trains using the function from Hideaki Shimazaki and Shigeru Shinomoto 2010
% Kernel Bandwidth Optimization in Spike Rate Estimation 
% Journal of Computational Neuroscience 2010
% http://dx.doi.org/10.1007/s10827-009-0180-4

% Setting the vector of points in ms at which an output value will be given
% for the filtered spike patterns
Tin = (1:Ntimebins)-0.5;

% Calculating the filered spike pattern from all trials
[y,t,~,~,~,~,~] = ssvkernel(Spike_ArrivalT_AllT,Tin);
% check that the input Tin was correctly used
if sum(Tin == t)~=length(Tin)
    error('WARNING: line 88 the kernel density estimation is using a different set of observation time points to return the filtered spike pattern\nThe # of time points used by ssvkernel is %d when we are asking for %d.\n', length(t),length(Tin));
end
% y is a density function that sums to 1
% multiplying by the total number of spikes gives the number of expecting spike per time bin (here 1 ms) for all 10 Trials
% dividing by the number of trials give the expected number os spikes per time bin for a trial
% multiplying by the response sampling rate in kHz gives the expected spike rate to one stimulus presentation in spike/ms
Psth_KDEfiltered =  y * sum(Spike_count(Local_trials)) / Ntrials * Response_samprate/1000;

% Calculating the filered spike pattern from Ntrials-1 trials to obtain JK
% (Jackknife) filtered spike patterns
JK_KDEfiltered = nan(Ntrials,Ntimebins);
for it = 1:Ntrials
    Local_trials = 1:Ntrials;
    Local_trials(it) = [];
    [y,t,~,~,~,~,~] = ssvkernel(Spike_ArrivalT_JK{it},Tin);
    if sum(Tin == t)~=length(Tin)
        error('WARNING: line 105 the kernel density estimation for trial %d is using a different set of observation time points to return the filtered spike pattern\nThe # of time points used by ssvkernel is %d when we are asking for %d.\n', it, length(t),length(Tin));
    end
    JK_KDEfiltered(it,:) = y * sum(Spike_count(Local_trials)) / (Ntrials-1) * Response_samprate/1000;
end

%% Plot figures if requested
if Fig
    figure(6)
    subplot(3,1,1)
    shadedErrorBar(Tin,Psth_KDEfiltered, std(JK_KDEfiltered,1),{'Color', 'k','LineStyle','-', 'LineWidth',1},1)
    hold on
    for tt=1:Ntrials
        Local_trial = find(Spike_array(tt,:));
        for ss=1:length(Local_trial)
            hold on
            plot([Local_trial(ss) Local_trial(ss)], [-0.05-0.02*tt -0.03-0.02*tt], 'LineWidth', 2, 'Color', 'k')
        end
    end
    hold off
    xlabel('Time (ms)')
    ylabel('Kernel density estimates of spike pattern')
    
    % visually verify that the optimal kernel width is investigated by the
    % range of value tested
    subplot(3,1,2)
    [~,~,~,W,C,~,~] = sskernel(Spike_ArrivalT_AllT,Tin);
    plot(W, C,'Color', 'r', 'LineWidth', 2);
    xlabel('Kernel Width (# points)');
    ylabel('Cost function');
    YLim = get(gca,'YLim');
    subplot(3,1,3)
    plot(W, C,'Color', 'r', 'LineWidth', 2);
    xlabel('Kernel Width (# points)');
    ylabel('Cost function');
    xlim([0 20])
    ylim([min(YLim) min(YLim)+(max(YLim)-min(YLim))/4])
    pause(2);
end

%% elongate the smoothed psth vector if the section was shorter than Windur ms
if (Ntimebins/Response_samprate)<Win
    Psth_KDEfiltered=[Psth_KDEfiltered ones(1,Resdur-Ntimebins)*Spike_rate_bg/Response_samprate];
    JK_KDEfiltered = [JK_KDEfiltered repmat(Spike_rate_bg/Response_samprate,Ntrials,Resdur-Ntimebins)];
    if Fig
        figure(7)
        shadedErrorBar([],Psth_KDEfiltered, std(JK_KDEfiltered,1),{'Color', 'k','LineStyle','-', 'LineWidth',1},1)
        hold on
        for tt=1:Ntrials
            Local_trial = find(Spike_array(tt,:));
            for ss=1:length(Local_trial)
                hold on
                plot([Local_trial(ss) Local_trial(ss)], [-0.05-0.02*tt -0.03-0.02*tt], 'LineWidth', 2, 'Color', 'k')
            end
        end
        hold off
        xlabel('Time (ms)')
        ylabel('Kernel density estimates of spike pattern')
        pause(1);
    end
elseif (Ntimebins/Response_samprate)>Win
    error('WARNING line 157 length of psth = %d ms, larger than the expected response size of %d ms\n', Ntimebins, Resdur);
end

Stim_mean_rate = mean(Spike_count)./Stimdur;
Spike_std_rate = std(Spike_count)./Stimdur;
Spike_rate_stim = Spike_count./Stimdur;

end