function output = signalSynchCheck(testSignalFile, mesuredFile, searchW, fixInd)
startTic = tic;
[xr,fsr] = audioread(testSignalFile);
[xm,fsm] = audioread(mesuredFile);
info = audioinfo(testSignalFile);
xr = sum(xr,2);
sigInfo = split(info.Comment);
lUnit = str2double(sigInfo{2});
lSignal = length(xr);
tmpCover = fftfilt(ones(lSignal,1), abs(xm));
[~, maxIdx] = max(tmpCover);
topLoc = maxIdx - lSignal;
lBias = -searchW:searchW;
nRepetition = str2double(sigInfo{3});
avDiff = zeros(nRepetition - 1,length(lBias));
baseIdx = 1000:round(0.8*lUnit);
for ii = 1:nRepetition - 1
    for jj = 1:length(lBias)
        avDiff(ii, jj) = std(xm((ii-1)*lUnit+topLoc + baseIdx) ...
            - xm((ii-1)*lUnit+topLoc + baseIdx + lBias(jj) + lUnit));
    end
end
stdSig = std(xr);
stdMes = std(xm);
[minStd, minIdx] = min(mean(avDiff));
idShift = lBias(minIdx);
if  idShift == 0 && minStd/sqrt(stdSig*stdMes) < 0.1 && fsr == fsm
    isInSync = true;
else
    isInSync = false;
    if fixInd
        ttm = (1:length(xm))'/fsr;
        ttmFix = ttm*lUnit/(lUnit-idShift);
        output.xFix = interp1(ttm, xm, ttmFix, 'linear','extrap');
    end
end
output.isInSync = isInSync;
output.syncLevel = minStd/sqrt(2*stdSig*stdMes);
output.minLocation = idShift;
output.elapsedTime = toc(startTic);
end