function output = decayMap(iResp, fs)
startTic = tic;
fl = 1000/32*2^(-5/7);
fh = fs/2;
nInOct = 7*2;
mag = 4;%13.31;
fBankStr = sixTermFBdesign(fs, fl, fh, nInOct, mag);
fcList = fBankStr.fcList;
nChannels = length(fcList);
maxHalfL = fBankStr.fbank.filter(1).halfL;
nData = length(iResp);
rawMap = zeros(nData+maxHalfL*2+1, nChannels);
baseIdx = (1:nData);
for ii = 1:nChannels
    w = fBankStr.fbank.filter(ii).w;
    fOut = fftfilt(w, iResp);
    halfL = fBankStr.fbank.filter(ii).halfL;
    rawMap(baseIdx - halfL + maxHalfL, ii) = abs(fOut).^2;
end
[~, maxIdxIresp] = max(abs(iResp));
%txBuffer = ((1:nData+maxHalfL*2+1)'-maxIdxIresp-maxHalfL)/fs;
txBuffer = ((1:nData+maxHalfL*2+1)'- maxIdxIresp-maxHalfL)/fs;
%tDelMap = -0.02:0.001:1;
tDelMap = txBuffer(txBuffer>-0.001);
borderp = 0.0005;
tDelMap(tDelMap > borderp) = borderp*exp(tDelMap(tDelMap > borderp)/borderp)/exp(1);
tDelMap = tDelMap(tDelMap<txBuffer(end));
pwrDecayMap = zeros(length(tDelMap)-1, nChannels);
for ii = 1:length(tDelMap)-2
    pwrDecayMap(ii, :) = ...
        max(rawMap(txBuffer > tDelMap(ii) & txBuffer <= tDelMap(ii+2), :));
end
for ii = 1:nChannels
    pwrDecayMap(:, ii) = pwrDecayMap(:, ii) / max(pwrDecayMap(:, ii))*sqrt(fcList(ii));
end
output.spRangePowerResp = sum(rawMap(:, fcList>300 & fcList<3400),2);
output.spRangePowerResp = output.spRangePowerResp /max(output.spRangePowerResp(:));
output.txBuffer = txBuffer;
output.tDelMap = tDelMap;
output.pwrDecayMap = pwrDecayMap/max(pwrDecayMap(:));
output.pwrDecayMapMedian = median(output.pwrDecayMap(:));
output.rawMap = rawMap;
output.fBankStr = fBankStr;
output.fs = fs;
output.elapsedTime = toc(startTic);
end