function CalMatfilesScript_JEE(UT,AT,OW)

addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/classreg'));
addpath(genpath('/auto/k5/matlab714/toolbox/econ'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/matfile

if nargin == 0
    UT='ALL';
end
if nargin == 1
    AT = 'R'; % 'R' for rate analysis running calls_selectivity_LRI_dprime set to 'CM' to run calculus of confusion matrices or to 'MD' to run modelizations of data
end
if nargin == 2
    OW = 1;%OW is the overwriting switch can be 1 (matfiles of units that were already processed will be overwritten) or 0
end

ExistingFiles=0;
input_dir=pwd;
Subjects = dir(input_dir);

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        if strcmp(AT, 'R') || strcmp(AT, 'MD')
            matfiles=dir(fullfile(Idir,'WholeVoc*.mat'));
        elseif strcmp(AT, 'CM')||strcmp(AT,'CMV')||strcmp(AT,'FVL')
            matfiles=dir(fullfile(Idir,'FirstVoc*.mat'));
        elseif strcmp(AT, 'RD') || strcmp(AT, 'RDP')
            matfiles=dir(fullfile(Idir,'ConfMat*.mat'));
        elseif strcmp(AT, 'RDPV')||strcmp(AT, 'RDPVn')
            matfiles=dir(fullfile(Idir,'ConfVoi*.mat'));
        end
        lm=length(matfiles);
        SS_Ind=zeros(lm,1);
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        if strcmp(UT, 'SS')
            Indices=find(SS_Ind);
        elseif strcmp(UT, 'S')
            Indices=find(~SS_Ind);
        elseif strcmp(UT, 'ALL')
            Indices=1:lm;
        end
        LM=length(Indices);
        
        %retrieve files already computed
        if strcmp(AT, 'R')
            Calmatfiles = dir(fullfile(Idir, 'LRI_DP*.mat'));
            CM = length(Calmatfiles);
        elseif strcmp(AT, 'CM')
            Calmatfiles = dir(fullfile(Idir, 'ConfMat*.mat'));
            CM = length(Calmatfiles);
        elseif strcmp(AT, 'CMV')
            Calmatfiles = dir(fullfile(Idir, 'ConfVoi*.mat'));
            CM = length(Calmatfiles);
        elseif strcmp(AT,'MD')
            Calmatfiles = dir(fullfile(Idir, 'Models*.mat'));
            CM = length(Calmatfiles);
        elseif strcmp(AT, 'RD')
            Calmatfiles = dir(fullfile(Idir, 'RandMat*.mat'));
            CM = length(Calmatfiles);
        elseif strcmp(AT, 'RDP')
            Calmatfiles = dir(fullfile(Idir, 'RandPMat*.mat'));
            CM = length(Calmatfiles);
        elseif strcmp(AT, 'RDPV')||strcmp(AT, 'RDPVn')
            Calmatfiles = dir(fullfile(Idir, 'RandPVoi*.mat'));
            CM = length(Calmatfiles);
        elseif strcmp(AT, 'FVL')
            Calmatfiles = dir(fullfile(Idir, 'ConfMat*.mat'));
            CM = length(Calmatfiles);
            
        end
        
        
        
        fprintf('Creating Calfile for %s\n',Indiv);
        for hh=1:LM
            matfile=matfiles(Indices(hh));
            MatfilePath=fullfile(Idir, matfile.name);
            
            %Make sure the file does not already exists
            FE=0;
            for cc = 1:CM
                if strcmp(Calmatfiles(cc).name(end-27:end), matfile.name(end-27:end))
                    fprintf('%s has already a Matfile\n', matfile.name);
                    if OW==0
                        FE=1;%switch to file exist only if OverWriting is not requested
                    end
                    ExistingFiles=ExistingFiles+1;
                end
            end
            
            if FE==0
                if strcmp(AT, 'R')
                    fprintf('Calculating dprimes, Index of selectivity for %s\n', MatfilePath);
                    [Calfilename]=calls_selectivity_LRI_dprime(MatfilePath);
                elseif strcmp(AT, 'CM')
                    fprintf('Calculating confusion matrices for %s\n', MatfilePath);
                    [Calfilename]=calls_selectivity_ConfusionMatrix(MatfilePath);
                elseif strcmp(AT,'FVL')
                    fprintf('Calculating nb cut vs intact stims for %s\n', MatfilePath);
                    [Calfilename]=nb_cutStim_cal(MatfilePath);
                elseif strcmp(AT, 'CMV')
                    fprintf('Calculating Voice confusion matrices for %s\n', MatfilePath);
                    [Calfilename]=Voice_selectivity_ConfusionMatrix(MatfilePath);
                elseif strcmp(AT, 'MD')
                    fprintf('Calculating Models for %s\n', MatfilePath);
                    [Calfilename]=SpectroSemantic_Neuro_model(MatfilePath);
                elseif strcmp(AT, 'RD')
                    fprintf('Calculating random confusion matrices for %s\n', MatfilePath);
                    [Calfilename]=RandomConfusionMatrices_Cal(MatfilePath);
                elseif strcmp(AT, 'RDP')
                    fprintf('Calculating random Predictions confusion matrices for %s\n', MatfilePath);
                    [Calfilename]=RandomPredictConfusionMatrices_Cal(MatfilePath);
                elseif strcmp(AT, 'RDPV')
                    fprintf('Calculating random Predictions Voice confusion matrices for %s\n', MatfilePath);
                    [Calfilename]=RandomPredictVoiceConfusionMatrices_Cal(MatfilePath);
                elseif strcmp(AT, 'RDPVn')
                    fprintf('adding stim names in random Predictions Voice confusion matrices for %s\n', MatfilePath);
                    [Calfilename]=AddVoicenames(MatfilePath);
                end
            
                fprintf('done making calculus on %s\nData save under %s\n', MatfilePath, Calfilename);
                clear MatfilePath
            end
        end
    end
end

if OW==1
    sprintf('%d matfiles already existed and has been overwritten', ExistingFiles);
elseif OW==0
    sprintf('%d matfiles already existed and were kept unchanged', ExistingFiles);
end
exit
