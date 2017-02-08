function [Section_Trials,Spikesfiltered,Psth,GaussFiltered_psth, Stim_mean_rate, Spike_std_rate, Spike_rate_stim, HwidthSpikes] = spikeTimes_psth_gaussfilter_cal(begInd, endInd, stim_samprate,response,spike_rate_bg, Win,response_samprate, JackKnifeSwitch,Fig)
% here Section_Trials is now transformed from s to ms
%% This function calculates the PSTH, isolate arrival times of spikes...
...in reference to the begining of the section (begInd).
if nargin<8
    JackKnifeSwitch=1;
end

if nargin<9
    Fig=0;
end

if nargin<7
    response_samprate = 10000; %10 000 corresponds to 0.1 ms precision TDT sampling rate is little above 24kHz so that's good enough 
end
Resdur=Win*response_samprate; %expected size of the response vector given de response sampling rate: here the output gaussFiltered_psth is completed with avg background rate if needed to match a size of Win s. 
stimdur=(endInd+1-begInd)/stim_samprate;
ntimebins = round(response_samprate*(endInd+1-begInd)/stim_samprate); % Number of time bins with choosen sampling rate. For 1bin=1ms response_samprate=1000

%t=(1:ntimebins)./1000.0;
ntrials=length(response.trials);
Section_Trials=cell(ntrials,1);
spike_count = zeros(1,ntrials);
Spike_array = zeros(ntrials, ntimebins);
for it=1:ntrials
    trial = response.trials{it};
    spike_time_trials = trial.spikeTimes;
    spike_time1=(spike_time_trials>=(begInd/stim_samprate));
    spike_time2=(spike_time_trials<=(endInd/stim_samprate));
    SpikesIndices=find(spike_time1 .* spike_time2);
    spike_count(it)=length(SpikesIndices);
    spike_time_trials_section=spike_time_trials(SpikesIndices);
    Section_Trials{it}=(spike_time_trials_section+ (1-begInd)/stim_samprate).*1000; % here spike times arrival are referred to the begining of the section instead of the begining of the stim and in ms instead of s 
    for is=1:spike_count(it)
        time_ind = floor(Section_Trials{it}(is)*response_samprate/1000);
        if (time_ind < 1 || time_ind > ntimebins)
            fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d ntimebins = %d\n', response.number, it, time_ind, ntimebins);
            continue;
        end
        Spike_array(it,time_ind) = Spike_array(it,time_ind) +1;
    end
end
Psth = sum(Spike_array,1)./ntrials;

% Calculated the smoothed psth (convolve each spike with a gaussian which
% width is twice the distance to the closest spike considering all trials)
alpha_param=3;
if JackKnifeSwitch % ntrials gaussian filtered spike patterns are calculated using ntrials-1 trials and the gaussian filtered spike pattern using all trials is calcuated using a convolution on the concatenation of all trials
    Spike_array_JK = nan(ntrials+1,ntimebins);
    for it = 1:ntrials
        Spike_array_local = Spike_array;
        Spike_array_local(it,:)=[];
        Spike_array_JK(it,:) = sum(Spike_array_local,1);
    end
    Spike_array_JK(ntrials+1,:) = sum(Spike_array,1);
    [Spikesfiltered,~, HwidthSpikes] = gauss_filter_varying_window(Spike_array_JK,alpha_param, ceil([repmat(ntrials-1,ntrials,1); ntrials]./2), Fig);
    GaussFiltered_psth = Spikesfiltered(ntrials+1,:)./ntrials;% This contain the gaussian filtered spike pattern calculated with all trials.
    Spikesfiltered = Spikesfiltered(1:ntrials,:)./(ntrials-1);% This contain the ntrials gaussian filtered spike patterns calculated with ntrials-1 trials.
else % each trial is convolved with time varying gaussians and trials are averaged after convolution to obtain the PSTH
    [Spikesfiltered,~, HwidthSpikes] = gauss_filter_varying_window(Spike_array,alpha_param, Fig);
    GaussFiltered_psth = mean(Spikesfiltered,1);
end
if Fig
    figure(5)
    subplot(1,2,1)
    plot(GaussFiltered_psth,'LineWidth',2, 'Color','g')
    hold on
    for jj=1:ntrials
        plot(Spikesfiltered(jj,:), 'Color','k')
        hold on
        if jj==1
            plot(mean(Spikesfiltered,1), 'LineWidth',2,'Color','r')%this is just for the legend
            hold on
            legend('Actual spike rate','individual jackknife', 'Average jackknife spike rate')
        end
    end
    plot(mean(Spikesfiltered,1), 'LineWidth',2,'Color','r')
    hold off
    xlabel('Time bins')
    ylabel('Gaussian filtered spike patterns')
end
% elongate the smoothed psth vector if the section was shorter than Windur ms
if (ntimebins/response_samprate)<Win
    GaussFiltered_psth=[GaussFiltered_psth repmat(spike_rate_bg/response_samprate,1,Resdur-ntimebins)]; 
    Spikesfiltered = [Spikesfiltered repmat(spike_rate_bg/response_samprate,ntrials,Resdur-ntimebins)];
elseif (ntimebins/response_samprate)>Win
    fprintf(1,'Warning length of psth = %d ms, larger than %d ms\n', ntimebins, Resdur);
end
if Fig
    figure(5)
    subplot(1,2,2)
    plot(GaussFiltered_psth,'LineWidth',2, 'Color','g')
    hold on
    for jj=1:ntrials
        plot(Spikesfiltered(jj,:), 'Color','k')
        hold on
        if jj==1
            plot(mean(Spikesfiltered,1), 'LineWidth',2,'Color','r')%this is just for the legend
            hold on
            legend('Actual spike rate','individual jackknife', 'Average jackknife spike rate')
        end
    end
    plot(mean(Spikesfiltered,1), 'LineWidth',2,'Color','r')
    hold off
    xlabel('Time bins')
    ylabel('Gaussian filtered spike patterns')
end
Stim_mean_rate = mean(spike_count)./stimdur;
Spike_std_rate = std(spike_count)./stimdur;
Spike_rate_stim = spike_count./stimdur;

end