function output = musicTSPmeasurementRep(x, fs, nRep, chId)
% LICENSE: refer to LICENSE in this folder

%nRep = 4;
startTic = tic;
switch size(x, 2)
    case 2
        x = x(:,1) + x(:, 2);
    case 1
    otherwise
        output = [];
        disp('one or two channel is acceptable');
        return;
end
%fftl = 2^ceil(log2(length(x)));
fftl = length(x);
fx = (0:fftl-1)'/fftl*fs;
xf = fft(x, fftl);
xfFix = xf;
thLevel = max(abs(xf)) * 0.003;
xfFix(abs(xf)<thLevel) = xfFix(abs(xf)<thLevel) ./ abs(xfFix(abs(xf)<thLevel)) * thLevel;
xFix = real(ifft(xfFix));
xRef = zeros(fftl * (nRep + 2), 1);
for ii = 1:nRep + 2
   xRef((ii-1)*fftl + (1:fftl)) =  xFix;
end
switch chId
    case 'L'
        xRefSt = [xRef xRef*0];
    case 'R'
        xRefSt = [xRef*0 xRef];
end
player = audioplayer(xRefSt/max(abs(xRef))*0.8, fs, 24);
recorder = audiorecorder(fs, 24, 1, -1);
record(recorder);
playblocking(player);
xMes = getaudiodata(recorder);
xMes = xMes(:,1);
respF = zeros(fftl, nRep);
impResp = zeros(fftl, nRep);
meanRefF = zeros(fftl, 1);
meanMesF = zeros(fftl, 1);
meanRespF = zeros(fftl, 1);
meanImpResp = zeros(fftl, 1);
for ii = 1:nRep
    xRefSpec = fft(xRef(ii*fftl + (1:fftl)))/max(abs(xRef))*0.8;
    xMesSpec = fft(xMes(ii*fftl + (1:fftl)));
    respF(:, ii) = xMesSpec(:) ./ xRefSpec(:);
    impResp(:, ii) = real(ifft(respF(:, ii)));
    meanRespF = meanRespF + respF(:, ii);
    meanImpResp = meanImpResp + impResp(:, ii);
    meanRefF = meanRefF + xRefSpec;
    meanMesF = meanMesF + xMesSpec;
end
meanRespF = meanRespF/nRep;
meanImpResp = meanImpResp/nRep;
meanRefF = meanRefF/nRep;
meanMesF = meanMesF/nRep;
meanErrPwr = zeros(fftl, 1);
for ii = 1:nRep
    meanErrPwr = meanErrPwr + abs(respF(:, ii) - meanRespF) .^2;
end
meanErrPwr = meanErrPwr / (nRep - 1);
fu = fx * 2^(1/6);
fl = fx * 2^(-1/6);
fw = (fu-fl);
fw(1) = 1;
cumTtlPw = cumsum(fx(2)*abs(meanMesF).*abs(meanRefF));
cumPw = cumsum(fx(2)*abs(meanRespF).^2);
cumPwEr = cumsum(fx(2)*meanErrPwr);
smoothedResponse = (interp1(fx(:), cumPw, fu(:), 'linear','extrap') ...
    - interp1(fx(:), cumPw, fl(:), 'linear','extrap')) ./ fw(:);
smoothedError = (interp1(fx(:), cumPwEr, fu(:), 'linear','extrap') ...
    - interp1(fx(:), cumPwEr, fl(:), 'linear','extrap')) ./ fw(:);
output.meanImpResp = meanImpResp;
output.respF = respF;
output.meanRespF = meanRespF;
output.xfFix = xfFix;
output.xMes = xMes;
output.meanErrPwr = meanErrPwr;
output.smoothedResponse = smoothedResponse;
output.smoothedError = smoothedError;
%output.smoothedResponseW = smoothedResponseW;
%output.smoothedErrorW = smoothedErrorW;
output.cumTtlPw = cumTtlPw;
output.fAxis = fx;
output.tAxis = (1:fftl)'/fs;
output.elapsedTime = toc(startTic);
end