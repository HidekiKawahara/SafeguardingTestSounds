function output = analyzeSafeguardedMeasurement(testFile, mesFile)
startTic = tic;
fileInfoStr = audioinfo(testFile);
output.fileInfoStr = fileInfoStr;
if isempty(fileInfoStr.Comment) || ~contains(fileInfoStr.Comment, "DFTA5SP")
    disp("These files are not properly prepared.")
    output = [];
    return
end
sigInfo = split(fileInfoStr.Comment);
lSegment = str2double(sigInfo{2});
nRepetition = str2double(sigInfo{3});
nLoops = str2double(sigInfo{4});
[xs,fss] = audioread(testFile);
mm = sum(abs(xs));
if mm(1) > mm(2)
    xs = xs(:, 1);
else
    xs = xs(:, 2);
end
syncOut = signalSynchCheck(testFile, mesFile, 5, 0);
output.isInSync = syncOut.isInSync;
if syncOut.isInSync
    [xm,~] = audioread(mesFile);
else
    syncRev = signalSynchCheck(testFile, mesFile, 100, 1);
    xm = syncRev.xFix;
end
mesfileInfoStr = audioinfo(mesFile);
mesComItem = split(mesfileInfoStr.Comment);
msk = ones(length(xs),1);
[~, maxIdx] = max(fftfilt(msk, abs(xm)));
fsmBias = maxIdx - length(xs) - 800;
flatSpl = 20*log10(std(xm(maxIdx - length(xs) + (1:length(xs))))) ...
    + str2double(mesComItem{2});
headRoom = round(lSegment*0.7);
silent = xm(1:fsmBias);
unitsilentL = [silent;silent(end-1:-1:2)];
nMulti = ceil(lSegment/length(unitsilentL));
silentL = zeros(nMulti*length(unitsilentL),1);
for ii = 1:nMulti
    silentL((ii-1)*length(unitsilentL) + (1:length(unitsilentL))) = ...
        unitsilentL;
end
baseIdx = 1:lSegment;
rawAvLTI = zeros(lSegment, nLoops);
rawAvLTIW = zeros(lSegment, nLoops);
xMesF = zeros(lSegment, 1);
xRefF = zeros(lSegment, 1);
avsdRand = zeros(lSegment, 1);
avsdRandW = zeros(lSegment, 1);
xMesFPw = zeros(lSegment, nLoops);
for jj = 1:nLoops
    tmpXrefF = zeros(lSegment, 1);
    tmpXmesF = zeros(lSegment, 1);
    headL = (jj-1)*nRepetition*lSegment;
    for ii = 1:nRepetition - 1
        xRef = xs(headL+(ii-1)*lSegment + baseIdx + headRoom);
        xMes = xm(headL+(ii-1)*lSegment + baseIdx + headRoom + fsmBias);
        tmpXrefF = tmpXrefF + fft(xRef);
        tmpXmesF = tmpXmesF + fft(xMes);
    end
    rawAvLTI(:, jj) = tmpXmesF./tmpXrefF;
    xMesFPw(:, jj) = abs(tmpXmesF/(nRepetition - 1)).^2;
    rawAvLTIW(:, jj) = rawAvLTI(:, jj) .* xMesFPw(:, jj);
    sdRand = zeros(lSegment, 1);
    for ii = 1:nRepetition - 1
        xRef = xs(headL+(ii-1)*lSegment + baseIdx + headRoom);
        xMes = xm(headL+(ii-1)*lSegment + baseIdx + headRoom + fsmBias);
        sdRand = sdRand + abs(fft(xMes)./fft(xRef) - rawAvLTI(:, jj)).^2;
    end
    sdRand = sdRand/(nRepetition - 2);
    xMesF = xMesF + tmpXmesF/(nRepetition - 1);
    xRefF = xRefF + tmpXrefF/(nRepetition - 1);
    avsdRand = avsdRand + sdRand;
    avsdRandW = avsdRandW + sdRand.*xMesFPw(:, jj);
end
avLTI = mean(rawAvLTI, 2);
avLTIW = mean(rawAvLTIW, 2)./mean(xMesFPw, 2);
sdSigDep = zeros(lSegment, 1);
sdSigDepW = zeros(lSegment, 1);
for jj = 1:nLoops
    sdSigDep = sdSigDep + abs(rawAvLTI(:, jj) - avLTI) .^2;
    sdSigDepW = sdSigDepW ...
        + (abs(rawAvLTI(:, jj) - avLTI) .^2) .*xMesFPw(:, jj);
end
if nLoops >1
    %sdSigDep = sdSigDep/(nLoops-1);
    sdSigDepW = sdSigDepW./sum(xMesFPw, 2)*nLoops/(nLoops-1);
else
    %sdSigDep = sdSigDep + NaN;
    sdSigDepW = sdSigDepW + NaN;
end
%avsdRand = avsdRand/nLoops;
avsdRandW = avsdRandW ./ sum(xMesFPw ,2);
output.TesfFile = testFile;
output.fAxis = (0:lSegment-1)'/lSegment*fss;
output.tAxis = (0:lSegment-1)/fss;
output.xMesF = xMesF/nLoops;
output.xRefF = xRefF/nLoops;
output.silentF = fft(silentL(baseIdx));
output.averageLTI = avLTI;
output.averageImpulse = real(ifft(avLTI,'symmetric'));
iResp = output.averageImpulse;
output.smoothedPwrResp = fftfilt(hanning(round(0.001*fss)), abs(iResp).^2);
output.smoothedPwrRespEnd = ...
    mean(output.smoothedPwrResp(end:-1:end-round(0.4*fss)));
output.squaredRand = sdRand;
%output.smoothedLTI = constRBWsmoothingTr(abs(avLTI).^2, fss, 6);
%output.smoothedRand = constRBWsmoothingTr(avsdRand, fss, 6);
%output.smoothedSigDep = constRBWsmoothingTr(sdSigDep, fss, 6);
output.smoothedSilentF = constRBWsmoothingTr(abs(output.silentF./xRefF).^2, fss, 6);
output.smoothedLTIW = constRBWsmoothingTr(abs(avLTIW).^2, fss, 6);
output.smoothedRandW = constRBWsmoothingTr(avsdRandW, fss, 6);
output.smoothedSigDepW = constRBWsmoothingTr(sdSigDepW, fss, 6);
output.flatSpl = flatSpl;
output.elaspedTime = toc(startTic);
end