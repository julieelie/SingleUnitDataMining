% dataDir = '/Users/ellenzippi/Desktop/TheunissenLab/data/matfiles/'; 
dataDir = '/auto/tdrive/julie/k6/julie/h5/';
% matfileDir = '/Users/ellenzippi/Desktop/TheunissenLab/data/FirstVocMat/'; 
matfileDir = '/auto/tdrive/julie/k6/julie/matfile/FirstVocMat/';
% addpath(genpath('/Users/ellenzippi/Documents/MATLAB/tlab/src'));
addpath(genpath('/auto/fhome/julie/Code/SingleUnitDataMining'));
% dataFiles = dir(fullfile(dataDir)); %gets all files in struct
matfiles = dir(fullfile(matfileDir, '*.mat')); 
OldSite = 'None'; 
for k = 1:length(matfiles)
    matfileName = matfiles(k).name; 
    matfilePath = fullfile(matfileDir, matfileName);
    Res = load(matfilePath); 
    h5Name = strcat(matfileName(10:end-4), '.h5'); 
    h5Path = fullfile(dataDir, Res.subject, h5Name); 
    Inde = strfind(h5Name, 'e');
    CurrentSite = h5Name(1:Inde(2)-2); 
    if strcmp(OldSite, CurrentSite) == 0 %This is a new site!
        fprintf(1,'This is a new site\nMatfile=%s\nh5file=%s\n',matfilePath, h5Path) % Test the wrapper 
        [VocTypeHD] = VocTypeHD_Def_4Matfiles(h5Path, matfilePath);
        % Update the name of the oldsite
        OldSite = CurrentSite;
    elseif strcmp(OldSite, CurrentSite) ==1
        fprintf('Same site, copying VocTypeHD on matfile = %s\n', matfilePath) 
        Res.VocTypeHD = VocTypeHD; 
        save(matfilePath, '-struct', 'Res');
    end
end