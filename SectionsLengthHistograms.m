load('HistoDuration1Voc.mat')
% histogram of all the calls
subplot(2,5,1)
hist(HistDur,30)

% histogram per call category
IndicesAg=find(strcmp(VocaType, 'Ag'));
subplot(2,5,2)
hist(HistDur(IndicesAg))

IndicesBe=find(strcmp(VocaType, 'Be'));
subplot(2,5,3)
hist(HistDur(IndicesBe))

IndicesDC=find(strcmp(VocaType, 'DC'));
subplot(2,5,4)
hist(HistDur(IndicesDC))

IndicesDi=find(strcmp(VocaType, 'Di'));
subplot(2,5,5)
hist(HistDur(IndicesDi))

IndicesLT=find(strcmp(VocaType, 'LT'));
subplot(2,5,6)
hist(HistDur(IndicesLT))

IndicesNe=find(strcmp(VocaType, 'Ne'));
subplot(2,5,7)
hist(HistDur(IndicesNe))

IndicesTe=find(strcmp(VocaType, 'Te'));
subplot(2,5,8)
hist(HistDur(IndicesTe))

IndicesTh=find(strcmp(VocaType, 'Th'));
subplot(2,5,9)
hist(HistDur(IndicesTh))

IndicesSo=find(strcmp(VocaType, 'song'));
subplot(2,5,10)
hist(HistDur(IndicesSo))