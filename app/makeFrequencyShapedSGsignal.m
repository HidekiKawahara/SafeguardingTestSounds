function output = makeFrequencyShapedSGsignal(x, fs, thLevel, octBw)
if size(x,2) > 1
    x = x(:,1) + x(:,2);
end
lUnit = length(x);
fftl = lUnit;
xF = fft(x);
fx = (0:fftl-1)/fftl*fs;
fx(fx>fs/2) = fx(fx>fs/2)-fs;
fxU = abs(fx) * 2^(1/2/octBw);
fxL = abs(fx) * 2^(-1/2/octBw);
bwF = fxU-fxL;
bwF(1) = bwF(2);
cumsumAbsXF = cumsum(abs(xF)*fx(2));
meanAbsF = (interp1(fx, cumsumAbsXF, fxU, 'linear','extrap') ...
    - interp1(fx, cumsumAbsXF, fxL, 'linear','extrap')) ./ bwF;
xFSG = xF;
modTh = 10^(thLevel/20);
modCount = 0;
for ii = 1:fftl
    if abs(xF(ii)) < meanAbsF(ii)*modTh
        xFSG(ii) = modTh*meanAbsF(ii) * xF(ii) / abs(xF(ii));
        modCount = modCount + 1;
    end
end
xSG = real(ifft(xFSG,'symmetric'));
output.fx = fx;
output.fftl = fftl;
output.modCount = modCount;
output.meanAbsF = meanAbsF;
output.xF = xF;
output.xSG = xSG;
end