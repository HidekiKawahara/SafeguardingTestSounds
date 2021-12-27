function output = makeSafeguardedMonoMusic(loopUnitList, nRepetition,thresDb)
nUnit = length(loopUnitList);
loopStrct = struct;
unitAttributes = zeros(nUnit, 3);
for ii = 1:nUnit
    [x, fs] = audioread(loopUnitList{ii});
    loopStrct.unit(ii).signal = x;
    loopStrct.unit(ii).fs = fs;
    unitAttributes(ii, :) = [size(x) fs];
end
if sum(sum(abs(diff(unitAttributes)))) ~= 0
    disp(unitAttributes);
    disp('All unit has to have the same length, size, and sampling frequency.');
    output = [];
    return;
end
if unitAttributes(1, 2) == 2
    disp('Stereo unit is converted to monaural.');
    for ii = 1:nUnit
        x = loopStrct.unit(ii).signal;
        loopStrct.unit(ii).signal = (x(:,1) + x(:,2))/2;
    end
end
lUnit = unitAttributes(1,1);
y = zeros(nRepetition*lUnit*nUnit, 1);
for ii = 1:nUnit
    xf = fft(loopStrct.unit(ii).signal);
    xfFix = xf;
    thLevel = mean(abs(xf)) * 10^(thresDb/20);
    xfFix(abs(xf)<thLevel) = xfFix(abs(xf)<thLevel) ./ abs(xfFix(abs(xf)<thLevel)) * thLevel;
    xFix = real(ifft(xfFix));
    for jj = 1:nRepetition
        y((ii-1)*lUnit*nRepetition + (jj-1)*lUnit + (1:lUnit)) = ...
            xFix;
    end
end
output.loopStrcture = loopStrct;
output.safeGuardedSignal = y;
output.samplingFrequency = fs;
output.unitAttributes = unitAttributes;
end