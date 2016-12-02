% Read WholeVoc file
MatfilePath='WholeVoc_Site2_L1100R1450_e14_s0_ss1.mat';
Res = load(MatfilePath);
for ss=1:length(Res.Spectro)
    figure(100)
    STIM=reshape(Res.Spectro{ss}, length(Res.Spectrofo{ss}), length(Res.Spectroto{ss}));
    valc = max(abs(max(max(STIM))), abs(min(min(STIM))));
    if valc==0
        imagesc(Res.Spectroto{ss}*1000,Res.Spectrofo{ss}, STIM)
    else
        imagesc(Res.Spectroto{ss}*1000, Res.Spectrofo{ss}, STIM,[-valc valc])
    end
    xlabel('Time (ms)')
    ylabel('Frequencies')
    title(sprintf('Stim %d/%d',ss,length(Res.Spectro)))
    pause(1)
end
