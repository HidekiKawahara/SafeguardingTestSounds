function output = analyzeSafeguardedMusicTest(safeGSignal, mesSignal, lUnit, nRepetition, fs)
%%
startTic = tic;
lTotal = length(safeGSignal);
%%
locChecker = fftfilt(ones(lTotal,1), abs(mesSignal));
[~, mxp] = max(locChecker);
initWithMargin = mxp - lTotal - round(0.1*fs);
%%
xMes = mesSignal(initWithMargin + (1:lTotal));
marGinL = round(lUnit * 0.6);
nUnit = lTotal/(nRepetition*lUnit);
%%
fx = (0:lUnit-1)/lUnit*fs;
averageGain = zeros(lUnit,1);
pwErr = zeros(lUnit,1);
averageUnitGainRes = zeros(lUnit,nUnit);
for ii = 1:nUnit
    averageUnitGain = zeros(lUnit,1);
    for jj = 1:nRepetition-1
        initialLoc = marGinL + (ii-1)*lUnit*nRepetition + (jj-1)*nUnit;
        segRef = safeGSignal(initialLoc + (1:lUnit));
        segMes = xMes(initialLoc + (1:lUnit));
        tmpUnitGain = fft(segMes)./fft(segRef);
        averageUnitGain = averageUnitGain + tmpUnitGain;
    end
    averageUnitGain = averageUnitGain/(nRepetition-1);
    averageGain = averageGain + averageUnitGain;
    pwUnitErr = zeros(lUnit,1);
    for jj = 1:nRepetition-1
        initialLoc = marGinL + (ii-1)*lUnit*nRepetition + (jj-1)*nUnit;
        segRef = safeGSignal(initialLoc + (1:lUnit));
        segMes = xMes(initialLoc + (1:lUnit));
        tmpUnitGain = fft(segMes)./fft(segRef);
        pwUnitErr = pwUnitErr + abs(averageUnitGain-tmpUnitGain).^2;
    end
    pwUnitErr = pwUnitErr/(nRepetition-2);
    pwErr = pwErr + pwUnitErr;
    averageUnitGainRes(:, ii) = averageUnitGain;
end
averageGain = averageGain/nUnit;
pwErr = pwErr/nUnit;
sigDepErr = zeros(lUnit,1);
for ii = 1:nUnit
    sigDepErr = sigDepErr + abs(averageGain - averageUnitGainRes(:, ii)) .^2;
end
sigDepErr = sigDepErr/(nUnit - 1);
impulseResponse = real(ifft(averageGain));grid on;
%%
output.fftl = lUnit;
output.timeAxis = (1:lUnit)'/fs;
output.frequencyAxis = fx;
output.samplingFrequency = fs;
output.averageLTIspec = averageGain;
output.impulseResponse = impulseResponse;
output.randErrorPower = pwErr;
output.signalDependentErrorPower = sigDepErr;
output.elapsedTime = toc(startTic);
end