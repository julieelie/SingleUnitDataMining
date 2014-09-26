
%This scripte retrieve the duration of the first cut in each stim and
%create and save a histogramm table
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/classreg'));
addpath(genpath('/auto/k5/matlab714/toolbox/econ'));
cd /auto/k6/julie/matfile
input_dir=pwd;
Subjects = dir(input_dir);
NVocSet = 0;

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        fprintf(1,'%s\n', Indiv);
        matfiles=dir(fullfile(Idir,'WholeVoc_*.mat'));
        lm=length(matfiles);
        
        % select only one file per site
        Site_Ind=zeros(lm,1);
        site=1;
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, sprintf('Site%d',site)))
                fprintf(1,'file %s has been selected for site %d\n', matfiles(ff).name, site);
                site=site+1;
                Site_Ind(ff)=1;
            end
        end
        FinalIndices=find(Site_Ind);
        
        LS = length(FinalIndices);
        NVocSet=NVocSet + LS;
    end
end
HistDur=zeros(NVocSet*150,1);
Vocatype=cell(NVocSet*150,1);
vv=0;
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        fprintf('%s\n', Indiv);        
        matfiles=dir(fullfile(Idir,'WholeVoc_*.mat'));
        lm=length(matfiles);
        
        % select only one file per site
        Site_Ind=zeros(lm,1);
        site=1;
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, sprintf('Site%d',site)))
                fprintf(1,'file %s has been selected for site %d\n', matfiles(ff).name, site);
                site=site+1;
                Site_Ind(ff)=1;
            end
        end
        FinalIndices=find(Site_Ind);
        
        LS = length(FinalIndices);
        
        for hh=1:LS
            MatfilePath=fullfile(Idir, matfiles(FinalIndices(hh)).name);
            fprintf(1,'Retrieving duration of first vocalization of each stim for %s\n', MatfilePath);
            Res=load(MatfilePath);
            Dur=Res.VocDuration;
            VocType=Res.VocType;
            DD=length(Dur);
            for dd=1:DD
                vv=vv+1;
                HistDur(vv)=Dur{dd}(1);
                Vocatype{vv}=VocType{dd};
            end
            fprintf(1,'done making calculus on %s\n', MatfilePath);
            clear MatfilePath
        end
    end
end
HistDur=HistDur(1:vv);
VocaType=Vocatype(1:vv);
filename=fullfile('/auto','k6','julie','matfile','HistoDuration1Voc.mat');
save(filename, 'HistDur','VocaType');
exit

figure(1)
hist(Hist.HistDur, 250);
figure(2)
for vv = 1:length(VT)
subplot(3,3,vv)
ind = find(strcmp(Hist.VocaType, VT{vv}));
hist(Hist.HistDur(ind))
title(sprintf('histogram of %s calls duration',VT{vv}));
fprintf('Min %s: %fms\n', VT{vv}, 1000*min(Hist.HistDur(ind)));
fprintf('Max %s: %fms\n',VT{vv}, 1000*max(Hist.HistDur(ind)));
fprintf('Mean %s: %fms\n',VT{vv}, 1000*mean(Hist.HistDur(ind)));
fprintf('Median %s: %fms\n',VT{vv}, 1000*median(Hist.HistDur(ind)));
fprintf('95%% quantile %s: %fms\n\n',VT{vv}, 1000.*quantile(Hist.HistDur(ind), 0.95));
end

qq = (0.1 : 0.05 : 1);
figure(3)
subplot(1,2,1)
for jj=1:length(qq)
plot(qq(jj), 1000*quantile(Hist.HistDur, qq(jj)), 'ko')
hold on
end
hline(250, 'r')
hline(600, 'g')
hold off
subplot(1,2,2)
for jj=1:length(qq)
plot(1000*quantile(Hist.HistDur, qq(jj)), qq(jj), 'ko')
hold on
end
vline(250, 'r')
vline(600, 'g')
hold off


figure(4)
for vv = 1:length(VT)
    subplot(3,3,vv)
    ind = find(strcmp(Hist.VocaType, VT{vv}));
    for jj=1:length(qq)
        plot(1000/quantile(1./Hist.HistDur(ind), qq(jj)),qq(jj), 'ko')
        hold on
        axis([0 2500 0 1])
    end
    hold off
end

figure(5)
ColOrd = get(gca,'ColorOrder');
cc=0;
for vv = 1:length(VT)
    ind = find(strcmp(Hist.VocaType, VT{vv}));
    cc=cc+1;
    %col = ColOrd(cc,:);
    for jj=1:length(qq)
        plot(1000/quantile(1./Hist.HistDur(ind), qq(jj)),qq(jj), 'o')
        hold on
        axis([0 2500 0 1])
    end
    hold on
end
    

