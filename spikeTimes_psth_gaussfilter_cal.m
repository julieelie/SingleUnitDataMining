function [Section_Trials,Spikesfiltered,Psth,GaussFiltered_psth, Stim_mean_rate, Spike_std_rate, Spike_rate_stim, HwidthSpikes] = spikeTimes_psth_gaussfilter_cal(begInd, endInd, samprate,response,spike_rate_bg, Win, Fig)
% here Section_Trials is now transformed from s to ms
%% This function calculates the PSTH, isolate arrival times of spikes...
...in reference to the begining of the section (begInd).
if nargin<7
    Fig=0;
end
Windur=Win*1000; %here the output gaussFiltered_psth is completed with avg background rate if needed to match a size of Windur ms
stimdur=(endInd+1-begInd)/samprate;
ntimebins = round(1000*(endInd+1-begInd)/samprate); % Number of time bins in ms

%t=(1:ntimebins)./1000.0;
ntrials=length(response.trials);
Section_Trials=cell(ntrials,1);
spike_count = zeros(1,ntrials);
Spike_array = zeros(ntrials, ntimebins);
for it=1:ntrials
    trial = response.trials{it};
    spike_time_trials = trial.spikeTimes;
    spike_time1=(spike_time_trials>=(begInd/samprate));
    spike_time2=(spike_time_trials<=(endInd/samprate));
    SpikesIndices=find(spike_time1 .* spike_time2);
    spike_count(it)=length(SpikesIndices);
    spike_time_trials_section=spike_time_trials(SpikesIndices);
    Section_Trials{it}=(spike_time_trials_section+ (1-begInd)/samprate).*1000; % here spike times arrival are referred to the begining of the section instead of the begining of the stim and in ms instead of s 
    for is=1:spike_count(it)
        time_ind = ceil(Section_Trials{it}(is));
        if (time_ind < 1 || time_ind > ntimebins)
            fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d ntimebins = %d\n', response.number, it, time_ind, ntimebins);
            continue;
        end
        Spike_array(it,time_ind) = Spike_array(it,time_ind) +1;
    end
end
Psth = sum(Spike_array,1)./ntrials;

% Calculated the smoothed psth (convolve each spike with a gaussian which
% width is twice the distance to the closest spike
alpha_param=3;
[Spikesfiltered,~, HwidthSpikes] = gauss_filter_varying_window(Spike_array,alpha_param, Fig);
GaussFiltered_psth = mean(Spikesfiltered,1);

% elongate the smoothed psth vector if the section was shorter than Windur ms
if ntimebins<Windur
    GaussFiltered_psth=[GaussFiltered_psth repmat(spike_rate_bg/1000,1,Windur-ntimebins)]; 
    Spikesfiltered = [Spikesfiltered repmat(spike_rate_bg/1000,ntrials,Windur-ntimebins)];
elseif ntimebins>Windur
    fprintf(1,'Warning length of psth = %d ms, larger than %d ms\n', ntimebins, Windur);
end

Stim_mean_rate = mean(spike_count)./stimdur;
Spike_std_rate = std(spike_count)./stimdur;
Spike_rate_stim = spike_count./stimdur;

end