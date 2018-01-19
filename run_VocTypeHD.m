% dataDir = '/Users/ellenzippi/Desktop/TheunissenLab/data/matfiles/'; 
dataDir = '/auto/tdrive/julie/k6/julie/h5/';
% matfileDir = '/Users/ellenzippi/Desktop/TheunissenLab/data/FirstVocMat/'; 
matfileDir = '/auto/tdrive/julie/k6/julie/matfile/FirstVocMat/';
% addpath(genpath('/Users/ellenzippi/Documents/MATLAB/tlab/src'));
addpath(genpath('/auto/fhome/julie/Code/SingleUnitDataMining'));
dataFiles = dir(fullfile(dataDir)); %gets all files in struct
OldSite = 'None'; 
for k = 1:length(dataFiles)
  if length(dataFiles(k).name)==11
      h5Dir = fullfile(dataDir, dataFiles(k).name);
      h5Files = dir(fullfile(h5Dir, '*ss*.h5')); 
      for j = 1:length(h5Files)
          h5Name = h5Files(j).name; 
          h5Path = fullfile(h5Dir, h5Name);
          Inde = strfind(h5Name, 'e');
          CurrentSite = h5Name(1:Inde(2)-1);
          matfileName = strcat('FirstVoc_', h5Name(1:end-3), '.mat'); 
          matfilePath = fullfile(matfileDir, matfileName);  
          %unit = read_unit_h5file(h5Path, 'r');
%           if exist(matfilePath) == 2 % This is just because I didn't have all the corresponding files (prevents errors for testing)
%               disp(sprintf('Current Unit: %s, Past Unit: %s', unit.site, site)); 
              if strcmp(OldSite, CurrentSite) == 0 %This is a new site!
%                   site = unit.site; 
%                   disp('The code would run here...') % Test the wrapper 

                  [VocTypeHD] = VocTypeHD_Def_4Matfiles(h5Path, matfilePath);
                  % Update the name of the oldsite
                  OldSite = CurrentSite;
              elseif strcmp(OldSite, CurrentSite) ==1
                  disp('Same site, copying VocTypeHD...') 
                  Res=load(matfilePath);
                  Res.VocTypeHD = VocTypeHD; 
                  save(matfilePath, '-struct', 'Res');
              end
%           end
      end
  end
end