%% Music TSP prototype test for AST Fostex measured at 10cm, 3rdTime
%% 
% 
%% Clear work space

close all
clear variables
%% Load three loops

[x1, fs] = audioread("holiday.wav");
[x2, ~] = audioread("laboratory.wav");
[x3, ~] = audioread("legend.wav");
[x4, ~] = audioread("space.wav");
%% Single shot impulse response measurement

if 1 == 2% recording disabled
recorder = audiorecorder(fs, 24, 1, -1);
recordblocking(recorder,length(x1)*2/fs);
xSilent = getaudiodata(recorder);
output = musicTSPmeasurementRep(x1, fs, 4, 'L');
end
%% Calibration tone recording, static multiple sinusoids 94.0dB at microphone

if 1 == 2
recorder = audiorecorder(fs, 24, 1, -1);
recordblocking(recorder,3);
xStatNoise = getaudiodata(recorder);
end
%% Save snapshot then comment out

%save snap12021Dec170610.mat
%save snap12021Dec171630.mat
%%
load snap12021Dec171630.mat
%% Calibration level check

aWeightFilter = weightingFilter("A-weighting",fs);
xStatNoiseW = aWeightFilter(xStatNoise);
measuredLevel = 20*log10(std(xStatNoiseW(fs + (1:fs))));
calibrationLevel = 94.0 - measuredLevel;
disp("Calibration level:" + num2str(calibrationLevel) + " dB");
xSilentW = aWeightFilter(xSilent);
bgLevel = 20*log10(std(xSilentW(fs + (1:fs)))) + calibrationLevel;
disp("A-weighted backgroud noise level:" + num2str(bgLevel) + " dB");
xMesW = aWeightFilter(output.xMes);
mesLevel = 20*log10(std(xMesW(fs + (1:fs)))) + calibrationLevel;
disp("A-weighted signal level:" + num2str(mesLevel) + " dB");
%% Display results

disp(output)
%%
fx = output.fAxis;
tx = output.tAxis;
xSeg1 = x1;
xSeg1F = fft(xSeg1);
lSig = length(xSeg1);
refDb = 10*log10(abs(xSeg1F).^2/length(xSeg1)^2);
refFixDb = 10*log10(abs(output.xfFix).^2/length(xSeg1)^2);
figure;
set(gcf,'position', [3501         642         700         240])
plot(sort(refDb), (1:lSig)/lSig, 'LineWidth', 2);
hold all
plot(sort(refFixDb), (1:lSig)/lSig, 'LineWidth', 2);
set(gca, 'FontSize', 16, "LineWidth", 2)
xline(min(refDb(:,1)),'k-', 'LineWidth', 2, "Label",['L-ch: ' num2str(min(refDb(:,1)),'%6.2f') ' dB'], 'FontSize', 16 ...
    ,'LabelHorizontalAlignment','left')
xline(min(refDb(:,2)),'k-', 'LineWidth', 2, "Label",['R-ch: ' num2str(min(refDb(:,2)),'%6.2f') ' dB'], 'FontSize', 16)
xline(max(refDb(:,1)),'k-', 'LineWidth', 2, "Label",['L: ' num2str(max(refDb(:,1)),'%6.2f') ' dB'], 'FontSize', 16)
xline(max(refDb(:,2)),'k-', 'LineWidth', 2, "Label", ...
    ['R: ' num2str(max(refDb(:,2)),'%6.2f') ' dB'], 'FontSize', 16 ,'LabelHorizontalAlignment','left')
text(-115, 0.3, 'original','fontsize',16)
text(-71, 0.3, 'safe-guarded','fontsize',16)
xlabel('level (dB rel. MSB)')
ylabel('cumulative probability')
grid on;
print -depsc specLevelDistribution.eps
%%
figure;
set(gcf,'position', [3501         642         700         240])
semilogx(fx, refDb, "LineWidth", 2);grid on;
hold all
semilogx(fx, refFixDb, "LineWidth", 2);grid on;
set(gca, 'FontSize', 16, "LineWidth", 2)
axis([20 fs/2 -140 -20])
ylabel('level (dB rel. MSB)')
xlabel('freqiemcu (Hz)')
legend('original L-ch', 'original R-ch', 'safe-guarded mono', 'Location','southwest', 'FontSize', 16)
print -depsc safeGuardedSpectrum.eps
%%
figure;
set(gcf,'position', [3501         642         700         300])
xMesSeg1 = output.xMes(lSig+(1:lSig));
mesDb = 10*log10(abs(fft(xMesSeg1)).^2/length(xSeg1)^2)+calibrationLevel+3;
xSilentSeg = xSilent((1:lSig));
silentDb = 10*log10(abs(fft(xSilentSeg)).^2/length(xSeg1)^2)+calibrationLevel+3;
semilogx(fx, mesDb, "LineWidth", 2);grid on
hold all
semilogx(fx, silentDb, "LineWidth", 2);grid on
set(gca, 'FontSize', 16, "LineWidth", 2)
axis([20 fs/2 -20 90])
ylabel({' ';'SPL (dB rel. 20 \mu Pa)'})
xlabel('freqiemcu (Hz)')
text(1000, 77, 'loudspeaker output','fontsize',16)
text(1000, 10, 'background noise','fontsize',16)
print -depsc measuredSigAndBg.eps
%%
figure;
set(gcf,'position', [3501         642         700         300])
xMesSeg1 = output.xMes(lSig+(1:lSig));
mesDb = 10*log10(abs(fft(xMesSeg1)).^2/length(xSeg1)^2);
xSilentSeg = xSilent((1:lSig));
silentDb = 10*log10(abs(fft(xSilentSeg)).^2/length(xSeg1)^2);
refErrDb = 10*log10(output.meanErrPwr);
semilogx(fx, mesDb - refFixDb, "LineWidth", 2);grid on
hold all
semilogx(fx, refErrDb, "LineWidth", 2);grid on
semilogx(fx, silentDb - refFixDb, "LineWidth", 2);grid on
axis([20 fs/2 -80 0])
set(gca, 'FontSize', 16, "LineWidth", 2)
ylabel({' ';'gain (dB)'})
xlabel('freqiemcu (Hz)')
text(1000, -20, 'LTI-response','fontsize',16)
text(2000, -45, 'random-TV response','fontsize',16)
text(1000, -75, 'background noise effect','fontsize',16)
print -depsc rawResponse.eps
%%

figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fx, 20*log10(output.meanRespF)); grid on;
hold all
semilogx(fx, 10*log10(output.meanErrPwr))
axis([20 fs/2 -80 20])
xlabel('frequency (Hz)')
ylabel('relative level (dB)')
%%
tmpPw = 10 .^ ((silentDb - refFixDb)/10);
bgEffectDb = constRBWsmoothing(tmpPw, fs, 3);
figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fx, 10*log10(output.smoothedResponse), "LineWidth", 2); grid on;
hold all
semilogx(fx, 10*log10(output.smoothedError), "LineWidth", 2)
semilogx(fx, 10*log10(bgEffectDb), "LineWidth", 2)
set(gca, 'FontSize', 16, "LineWidth", 2)
%%
axis([20 fs/2 -80 0])
xlabel('frequency (Hz)')
ylabel({' ';'gain (dB)'})
text(1000, -12, 'LTI-response','fontsize',16)
text(1000, -30, 'random-TV response','fontsize',16)
text(1000, -55, 'background noise effect','fontsize',16)
print -depsc smoothedOneTypResponse.eps
%%

tmpPw = 10 .^ ((silentDb - refFixDb)/10);
bgEffectDb = constRBWsmoothing(tmpPw, fs, 3);
sampleFxIdx = round(50*2.^(0:1/3:log2(fs/2/50))/fx(2))+1;
fixOneThird = -10*log10(length(sampleFxIdx));
meanPw = mean(output.smoothedResponse(sampleFxIdx));
mesLevel = 20*log10(std(aWeightFilter(output.xMes))) + calibrationLevel;
figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fx, 10*log10(output.smoothedResponse/meanPw)+mesLevel+fixOneThird, "LineWidth", 2); grid on;
hold all
semilogx(fx, 10*log10(output.smoothedError/meanPw)+mesLevel+fixOneThird, "LineWidth", 2)
semilogx(fx, 10*log10(bgEffectDb/meanPw)+mesLevel+fixOneThird, "LineWidth", 2)
set(gca, 'FontSize', 16, "LineWidth", 2);
%axis([20 fs/2 [-80 20] + mesLevel+fixOneThird]);
axis([20 fs/2 0 100]);
xlabel('frequency (Hz)')
ylabel({' ';'virtual 1/3 oct SPL (dB)'})
text(1000, 90, 'LTI-response','fontsize',16)
text(2000, 56, 'random-TV response','fontsize',16)
text(2000, 32, 'background noise effect','fontsize',16)
title("SPL (A-weight): " + num2str(mesLevel,'%5.1f') + " dB")
print -depsc smoothedOneTypResponseV3L.eps

%% Non-linear response separation test

%segTopLocs = [19.7, 23.6, 27.5, 30.2, 35.4, 39.3, 42.8, 99.3];
if 1 == 2% disable not to destroy
mesSet = struct;
for ii = 1:4
    switch ii
        case 1
            segX = x1;
        case 2
            segX = x2;
        case 3
            segX = x3;
        case 4
            segX = x4;
    end
    output = musicTSPmeasurementRep(segX, fs, 4, 'L');
    mesSet.rec(ii) = output;
end
end
%%
%load mesSet2021dec140252.mat
%save mesSet2021dec170625.mat mesSet
%save mesSet2021dec170639.mat mesSet
%save mesSet2021dec170823.mat mesSet
%save mesSet2021dec171640.mat mesSet
%save mesSet2021dec171653.mat mesSet
load mesSet2021dec171653.mat
%%
sigLevelRMS = 0;
for ii = 1:4
    figure;
    set(gcf,'position', [3501         642         700         300])
    semilogx(fx, 10*log10(mesSet.rec(ii).smoothedResponse), "LineWidth", 2); 
    grid on;
    hold all;
    semilogx(fx, 10*log10(mesSet.rec(ii).smoothedError), "LineWidth", 2); 
    set(gca, 'FontSize', 16, "LineWidth", 2)
    axis([20 fs/2 -60 20])
    xlabel('frequency (Hz)')
    ylabel({' ';'gain (dB)'})
    sigLevel = 20*log10(std(aWeightFilter(mesSet.rec(ii).xMes))) + calibrationLevel;
    sigLevelRMS = sigLevelRMS + 10^(sigLevel/10);
    switch ii
        case 1
            title(['holiday.wav  SPL(A):' num2str(sigLevel,'%5.1f') ' dB' ])
        case 2
            title(['laboratory.wav  SPL(A):' num2str(sigLevel,'%5.1f') ' dB' ])
        case 3
            title(['legend.wav  SPL(A):' num2str(sigLevel,'%5.1f') ' dB' ])
        case 4
            title(['space.wav  SPL(A):' num2str(sigLevel,'%5.1f') ' dB' ])
    end
end
averageSigLevel = 10*log10(sigLevelRMS/4);
disp("average SPL(A):" + num2str(averageSigLevel,'%5.1f') + " dB")
%%
if 1 == 2
angg = (1:lSig)/lSig*2*pi;
ddd = zeros(length(angg),1);
for ii = 1:length(angg)
    ddd(ii) = std(abs(mesSet.rec(1).meanRespF*exp(1i*angg(ii)) - mesSet.rec(4).meanRespF));
end
end
%%

hsLTIresp = zeros(lSig,1);
erVTrespSq = zeros(lSig,1);
for ii = 1:4
    if ii ==4
        hsLTIresp = hsLTIresp + mesSet.rec(ii).meanRespF.*exp(1i*unitx*1024);
    else
    hsLTIresp = hsLTIresp + mesSet.rec(ii).meanRespF;
    end
    erVTrespSq = erVTrespSq + mesSet.rec(ii).meanErrPwr;
end
hsLTIresp = hsLTIresp/4;
erVTrespSq = erVTrespSq/4;
ernLTIrespSq = zeros(lSig,1);
for ii = 1:4
    if ii ==4
        ernLTIrespSq = ernLTIrespSq + abs(mesSet.rec(ii).meanRespF.*exp(1i*unitx*1024)-hsLTIresp).^2;
    else
    ernLTIrespSq = ernLTIrespSq + abs(mesSet.rec(ii).meanRespF-hsLTIresp).^2;
    end
end
ernLTIrespSq = ernLTIrespSq/(length(segTopLocs)-1);
%%
figure;
semilogx(fx, 20*log10(abs(hsLTIresp))); grid on
hold all
semilogx(fx, 10*log10(abs(erVTrespSq))); grid on
semilogx(fx, 10*log10(abs(ernLTIrespSq))); grid on
axis([20 fs/2 -60 20])
xlabel('frequency (Hz)')
ylabel('relative level (dB)')
%%
figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fx, 10*log10(constRBWsmoothing(abs(hsLTIresp).^2, fs, 3)), 'LineWidth', 2); grid on
hold all
semilogx(fx, 10*log10(constRBWsmoothing(erVTrespSq, fs, 3)), 'LineWidth', 2); grid on
semilogx(fx, 10*log10(constRBWsmoothing(ernLTIrespSq, fs, 3)), 'LineWidth', 2); grid on
set(gca, 'FontSize', 15, 'LineWidth', 2)
axis([20 fs/2 -60 20])
xlabel('frequency (Hz)')
ylabel('gain (dB)')
outFname = 'fostex10cm85db.eps';
%outFname = 'fostex10cm95db.eps';
switch outFname
    case 'fostex10cm85db.eps'
        text(1000, -4, 'LTI response', 'FontSize', 15)
        text(1000, -30, 'signal dependent response', 'FontSize', 15)
        text(1000, -56, 'random-TV response', 'FontSize', 15)
        print -depsc fostex10cm85db.eps
    case 'fostex10cm95db.eps'
        text(1000, 4, 'LTI response', 'FontSize', 15)
        text(1000, -18, 'signal dependent response', 'FontSize', 15)
        text(1000, -50, 'random-TV response', 'FontSize', 15)
        print -depsc fostex10cm95db.eps
end
%%
sampleFxIdx = round(50*2.^(0:1/3:log2(fs/2/50))/fx(2))+1;
fixOneThird = -10*log10(length(sampleFxIdx));
smoothedLTIrespSq = constRBWsmoothing(abs(hsLTIresp).^2,fs,3);
meanPw = mean(smoothedLTIrespSq(sampleFxIdx));
%mesLevel = 20*log10(std(aWeightFilter(output.xMes))) + calibrationLevel;
fixLv = averageSigLevel+fixOneThird;
figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fx, 10*log10(constRBWsmoothing(abs(hsLTIresp).^2, fs, 3))+fixLv, 'LineWidth', 2); grid on
hold all
semilogx(fx, 10*log10(constRBWsmoothing(erVTrespSq, fs, 3))+fixLv, 'LineWidth', 2); grid on
semilogx(fx, 10*log10(constRBWsmoothing(ernLTIrespSq, fs, 3))+fixLv, 'LineWidth', 2); grid on
set(gca, 'FontSize', 15, 'LineWidth', 2)
axis([20 fs/2 0 100])
xlabel('frequency (Hz)')
ylabel('virtual 1/3 oct SPL (dB)')
outFname = 'fostex10cm85db.eps';
%outFname = 'fostex10cm95db.eps';
switch outFname
    case 'fostex10cm85db.eps'
        text(1000, 70, 'LTI response', 'FontSize', 15)
        text(1000, 44, 'signal dependent response', 'FontSize', 15)
        text(1000, 15, 'random-TV response', 'FontSize', 15)
        title("SPL (A-weight): " + num2str(averageSigLevel,'%5.1f') + " dB")
        print -depsc fostex10cm85dbV3L.eps
    case 'fostex10cm95db.eps'
        text(1000, 85, 'LTI response', 'FontSize', 15)
        text(1000, 62, 'signal dependent response', 'FontSize', 15)
        text(1000, 27, 'random-TV response', 'FontSize', 15)
        title("SPL (A-weight): " + num2str(averageSigLevel,'%5.1f') + " dB")
        print -depsc fostex10cm95dbV3L.eps
end
%% Measurent results saved then commented out

%save mesSet2021dec140252.mat mesSet
%% Simulation using white noise test signal and observation noise

lTestSig = 100000;
fxs = (0:lTestSig-1)/lTestSig*fs;
xTestOrg = randn(lTestSig, 1);
xBgNoise = randn(lTestSig, 1);
%% SNR 40dB without low level clipping

xRefF = fft(xTestOrg);
xNoise = fft(xBgNoise);
xMes = xTestOrg + xBgNoise/sqrt(10000);
xMesF = fft(xMes);
Hraw = xMesF ./ xRefF;
figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fxs, 10*log10(abs(xRefF).^2/lTestSig^2));grid on;
hold all
semilogx(fxs, 10*log10(abs(xNoise).^2/lTestSig^2/10000));grid on;
axis([20 fs/2 -140 0])

figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fxs, 20*log10(abs(Hraw)));grid on;
set(gca, 'xlim', [20 fs/2])
%% SNR 40dB with 0dB from average absolute value

xRefFixF = xRefF;
thetaL = mean(abs(xRefF))/sqrt(1);
xRefFixF(abs(xRefFixF) < thetaL) = ...
    thetaL * xRefFixF(abs(xRefFixF) < thetaL) ./ abs(xRefFixF(abs(xRefFixF) < thetaL));
xRefFix = real(ifft(xRefFixF));
xMesS = xRefFix + xBgNoise*std(xRefFix)/100;
xMesSF = fft(xMesS);
HrawS = xMesSF ./ xRefFixF;

figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fxs, 10*log10(abs(xRefFixF).^2/lTestSig^2));grid on;
hold all
semilogx(fxs, 10*log10(abs(xNoise*std(xRefFix)).^2/lTestSig^2/10000));grid on;
semilogx(fxs, 10*log10(abs(xNoise).^2/lTestSig^2/10000));grid on;
axis([20 fs/2 -140 0])

figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fxs, 20*log10(abs(HrawS)));grid on;
set(gca, 'xlim', [20 fs/2])
%% Figures for AST

fx = fxs;
figure;
set(gcf,'position', [3501         642         700         300])
semilogx(fx, 10*log10(abs(xRefF).^2/lTestSig^2),'LineWidth',2);grid on;
hold all
semilogx(fx, 10*log10(abs(xRefFixF).^2/lTestSig^2),'LineWidth',2);grid on;
semilogx(fx, 10*log10(abs(xNoise*std(xRefFix)).^2/lTestSig^2/10000),'LineWidth',2);grid on;
semilogx(fx, 10*log10(abs(xNoise).^2/lTestSig^2/10000),'LineWidth',2);grid on;
set(gca, 'FontSize', 15, 'LineWidth', 2)
axis([20 fs/2 -120 -30])
xlabel('frequency (Hz)')
ylabel('level (dB rel. MSB)')
text(3000, -49, 'X_s: safe','Color',[1 1 1], 'FontSize', 18)
text(8000, -62, 'X','Color',[1 1 1], 'FontSize', 18)
text(30, -78, 'adjusted additive noise','Color',[0 0 0], 'FontSize', 18)
text(2000, -95, 'original additive noise','Color',[1 1 1], 'FontSize', 18)
print -depsc whiteNoiseClipSpec.eps

figure;
set(gcf,'position', [3501         642         700         300], 'Color','none','InvertHardcopy','off')
semilogx(fx, 20*log10(abs(Hraw)),'LineWidth',2);grid on;
hold on
semilogx(fx, 20*log10(abs(HrawS)),'LineWidth',2);grid on;
set(gca, 'FontSize', 15, 'LineWidth', 2,'Color','none')
set(gca, 'xlim', [20 fs/2])
xlabel('frequency (Hz)')
ylabel('gain (dB)')
legend('raw H', 'safeguarded H','Location','northwest')
print -depsc whitenoiseRawAndSafeH.eps
%%

figure;
set(gcf,'position', [3501         642         700         220], 'Color','none','InvertHardcopy','off')
plot(sort(10*log10(abs(xRefF).^2/lTestSig^2)) , (1:lTestSig)/lTestSig,'LineWidth',2);
hold all
plot(sort(10*log10(abs(xRefFixF).^2/lTestSig^2)) , (1:lTestSig)/lTestSig,'LineWidth',2);
plot(sort(10*log10(abs(xNoise*std(xRefFix)).^2/lTestSig^2/10000)) , (1:lTestSig)/lTestSig,'LineWidth',2);
plot(sort(10*log10(abs(xNoise).^2/lTestSig^2/10000)) , (1:lTestSig)/lTestSig,'LineWidth',2);
set(gca, 'FontSize', 15, 'LineWidth', 2,'Color','none')
grid on;
xlabel('level (dB rel. MSB)')
ylabel('cumulative probability')
text(-48, 0.3, 'X_s: safe','Color',[0 0 0], 'FontSize', 18)
text(-60, 0.3, 'X','Color',[0 0 0], 'FontSize', 18)
text(-92, 0.3, 'adjusted','Color',[0 0 0], 'FontSize', 18)
text(-92, 0.22, 'additive noise','Color',[0 0 0], 'FontSize', 18)
text(-135, 0.3, 'original additive noise','Color',[0 0 0], 'FontSize', 18)
print -depsc whiteNoiseDistribution.eps
%% SNR and clipping level dependency

clipLevel = -50:2:20;
snrList = 10:10:50;
stdResult = zeros(length(clipLevel),length(snrList));
maxResult = zeros(length(clipLevel),length(snrList));
for jj = 1:length(snrList)
    for ii = 1:length(clipLevel)
        xRefFixF = xRefF;
        thetaL = mean(abs(xRefF))*10^(clipLevel(ii)/20);
        xRefFixF(abs(xRefFixF) < thetaL) = ...
            thetaL * xRefFixF(abs(xRefFixF) < thetaL) ./ abs(xRefFixF(abs(xRefFixF) < thetaL));
        xRefFix = real(ifft(xRefFixF));
        xMesS = xRefFix + xBgNoise*std(xRefFix)/10^(snrList(jj)/20);
        xMesSF = fft(xMesS);
        HrawS = xMesSF ./ xRefFixF;
        stdResult(ii,jj) = std(20*log10(abs(HrawS)));
        maxResult(ii,jj) = max(20*log10(abs(HrawS)));
    end
end
%%
figure;
set(gcf,'position', [3501         642         700         400])
semilogy(clipLevel,stdResult,'LineWidth',2);grid on;
set(gca, 'FontSize', 15, 'LineWidth', 2)
xlabel('clipping level (dB from the average original absolute spectrum)')
ylabel('standard deviation (dB)')
text(-40, 5, 'SNR 10 dB','Color',[0 0 0], 'FontSize', 18)
text(-40, 2, 'SNR 20 dB','Color',[0 0 0], 'FontSize', 18)
text(-40, 0.8, 'SNR 30 dB','Color',[0 0 0], 'FontSize', 18)
text(-40, 0.28, 'SNR 40 dB','Color',[0 0 0], 'FontSize', 18)
text(-40, 0.09, 'SNR 50 dB','Color',[0 0 0], 'FontSize', 18)
%%
figure;
set(gcf,'position', [3501         642         700         260])
semilogy(clipLevel,maxResult,'LineWidth',2);grid on;
set(gca, 'FontSize', 15, 'LineWidth', 2)
axis([-50 20 0.07 100])
xlabel('flooring level (dB from the average original absolute spectrum)')
ylabel('maximum deviation (dB)')
text(5, 10, 'SNR 10 dB','Color',[0 0 0], 'FontSize', 18)
text(5, 3.3, 'SNR 20 dB','Color',[0 0 0], 'FontSize', 18)
text(5, 1.2, 'SNR 30 dB','Color',[0 0 0], 'FontSize', 18)
text(5, 0.4, 'SNR 40 dB','Color',[0 0 0], 'FontSize', 18)
text(5, 0.13, 'SNR 50 dB','Color',[0 0 0], 'FontSize', 18)
print -depsc clipOnMaxDev.eps
%% Spectral smoothing effects on LTI distribution

lTestSig = 100000;
fxs = (0:lTestSig-1)/lTestSig*fs;
xTestOrg = randn(lTestSig, 1);
xBgNoise = randn(lTestSig, 1);
xRefFixF = xRefF;
thetaL = mean(abs(xRefF))/sqrt(1);
snr = 30;
xRefFixF(abs(xRefFixF) < thetaL) = ...
    thetaL * xRefFixF(abs(xRefFixF) < thetaL) ./ abs(xRefFixF(abs(xRefFixF) < thetaL));
xRefFix = real(ifft(xRefFixF));
xMesS = xRefFix + xBgNoise*std(xRefFix)/10^(snr/20);
xMesSF = fft(xMesS);
HrawS = xMesSF ./ xRefFixF;

%%
figure;
semilogx(fx, 20*log10(abs(HrawS)));grid on;
fxL = fxs * 2^(-1/6);
fxH = fxs * 2^(1/6);
fxL6 = fxs * 2^(-1/12);
fxH6 = fxs * 2^(1/12);
fw = (fxH-fxL);
fw(1) = fw(2);

%figure;
tic
nTest = 1000;
snrList = 10:10:50;
stdstdSmsFix = zeros(lTestSig, length(snrList));
stdstdSmsSmsFix = zeros(lTestSig, length(snrList));
for kk = 1:length(snrList)
    stdSmsFix = zeros(lTestSig, nTest);
    stdSmsSmsFix = zeros(lTestSig, nTest);
    snr = snrList(kk);
    for ii = 1:nTest
        xTestOrg = randn(lTestSig, 1);
        xBgNoise = randn(lTestSig, 1);
        xRefFixF = xRefF;
        thetaL = mean(abs(xRefF))/sqrt(1);
        xRefFixF(abs(xRefFixF) < thetaL) = ...
            thetaL * xRefFixF(abs(xRefFixF) < thetaL) ./ abs(xRefFixF(abs(xRefFixF) < thetaL));
        xRefFix = real(ifft(xRefFixF));
        xMesS = xRefFix + xBgNoise*std(xRefFix)/10^(snr/20);
        xMesSF = fft(xMesS);
        HrawS = xMesSF ./ xRefFixF;
        cumFixPwr = cumsum(abs(HrawS).^2*fxs(2));
        HbndPwr = (interp1(fxs,cumFixPwr,fxH,'linear','extrap') - interp1(fxs,cumFixPwr,fxL,'linear','extrap'))./fw;
        cumFixSmsPwr = cumsum(HbndPwr*fxs(2));
        HbndSmsPwr = (interp1(fxs,cumFixSmsPwr,fxH6,'linear','extrap') ...
            - interp1(fxs,cumFixSmsPwr,fxL6,'linear','extrap'))./fw;
        %semilogx(fx, 10*log10(HbndPwr));grid on;
        %hold all
        %drawnow
        stdSmsFix(:, ii) = 10*log10(HbndPwr);
        stdSmsSmsFix(:, ii) = 10*log10(HbndSmsPwr);
    end
    stdstdSmsFix(:, kk) = std(stdSmsFix');
    stdstdSmsSmsFix(:, kk) = std(stdSmsSmsFix');
end
toc
%%
figure;
set(gcf,'position', [3501         642         700         240])
loglog(fxs, stdstdSmsFix, 'LineWidth', 2);grid on;
set(gca, 'xlim', [20 fs/2])
set(gca, 'FontSize', 15, 'LineWidth', 2)
xlabel('frequency (Hz)')
ylabel('standard deviation (dB)')
text(1000, 0.13, 'SNR 10 dB','Color',[0 0 0], 'FontSize', 15)
text(1000, 0.03*1.3, 'SNR 20 dB','Color',[0 0 0], 'FontSize', 15)
text(1000, 0.01*1.3, 'SNR 30 dB','Color',[0 0 0], 'FontSize', 15)
text(1000, 0.003*1.3, 'SNR 40 dB','Color',[0 0 0], 'FontSize', 15)
text(1000, 0.001*1.3, 'SNR 50 dB','Color',[0 0 0], 'FontSize', 15)
print -depsc whiteNSmoothingSD.eps
%% Clipping effect on the original signal

clipLevel = -40:2:20;
waveDistortion = zeros(length(clipLevel), 1);
for ii = 1:length(clipLevel)
    xRefFixF = fft(xTestOrg);
    thetaL = mean(abs(xRefFixF))*10^(clipLevel(ii)/20);
    xRefFixF(abs(xRefFixF) < thetaL) = ...
        thetaL * xRefFixF(abs(xRefFixF) < thetaL) ./ abs(xRefFixF(abs(xRefFixF) < thetaL));
    xRefFix = real(ifft(xRefFixF));
    waveDistortion(ii) = std(xTestOrg-xRefFix);
end
%% 
% mdl = fitlm(clipLevel(5:end-8), 20*log10(waveDistortion(5:end-8)))
% 
% 
% 
% mdl = 
% 
% 
% 
% 
% 
% Linear regression model:
% 
% y ~ 1 + x1
% 
% 
% 
% Estimated Coefficients:
% 
% Estimate       SE        tStat       pValue  
% 
% ________    ________    _______    __________
% 
% 
% 
% (Intercept)    -10.321      0.16566    -62.305    1.6482e-21
% 
% x1               1.995     0.009319     214.08    1.3121e-30
% 
% 
% 
% 
% 
% Number of observations: 19, Error degrees of freedom: 17
% 
% Root Mean Squared Error: 0.445
% 
% R-squared: 1,  Adjusted R-Squared: 1
% 
% F-statistic vs. constant model: 4.58e+04, p-value = 1.31e-30


figure;
set(gcf,'position', [3501         642         700         300])
plot(clipLevel,20*log10(waveDistortion),'LineWidth',2);
grid on;
elvl = -10.321 + clipLevel*1.995;
hold all
plot(clipLevel,elvl,'LineWidth',2);
set(gca, 'FontSize', 15, 'LineWidth', 2)
xlabel('clipping level (dB)')
ylabel('noise level (dB)')
legend('observation', 'linear regression','Location','northwest', 'FontSize', 18)
text(-10, -60, 'y ~ -10.321 + 1.995 x', 'FontSize', 18)
print -depsc equivalentNoise.eps
%% Estimation of random response

clipLevel = -50:5:20;
snrList = 10:10:50;
nRepetition = 10;
HrawSRslt = zeros(lTestSig, nRepetition);
snrLevelEst = zeros(length(snrList), length(clipLevel));
for kk = 1:length(snrList)
    for ii = 1:length(clipLevel)
        averageResp = zeros(lTestSig, 1);
        for jj = 1:nRepetition
            xRefFixF = fft(xTestOrg);
            thetaL = mean(abs(xRefFixF))*10^(clipLevel(ii)/20);
            xRefFixF(abs(xRefFixF) < thetaL) = ...
                thetaL * xRefFixF(abs(xRefFixF) < thetaL) ./ abs(xRefFixF(abs(xRefFixF) < thetaL));
            xRefFix = real(ifft(xRefFixF));
            xMesS = xRefFix + std(xRefFix)*randn(lTestSig, 1)/10^(snrList(kk)/20);
            xMesSF = fft(xMesS);
            HrawS = xMesSF ./ xRefFixF;
            HrawSRslt(:, jj) = HrawS;
            averageResp = averageResp + HrawS;
        end
        averageResp = averageResp/nRepetition;
        sqNoise = 0;
        for jj = 1:nRepetition
            sqNoise = sqNoise + mean(abs(HrawSRslt(:, jj)-averageResp).^2);
        end
        sqNoiseLevel = sqNoise/(nRepetition -1);
        snrLevelEst(kk, ii) = sqNoise/(nRepetition -1);
    end
    %snrLevelEst(kk, :) = snrLevelEst(kk, :) / sqNoise / 10^(snrList(kk)/10);
end
%%
figure;
set(gcf,'position', [3501         642         700         230])
plot(clipLevel,10*log10(snrLevelEst), 'LineWidth', 2);grid on
set(gca, 'FontSize', 15, 'LineWidth', 2)
xlabel('flooring level (dB from the average original absolute spectrum)')
ylabel('level (dB)')
text(5, -7, 'SNR 10 dB','Color',[0 0 0], 'FontSize', 18)
text(5, -17, 'SNR 20 dB','Color',[0 0 0], 'FontSize', 18)
text(5, -27, 'SNR 30 dB','Color',[0 0 0], 'FontSize', 18)
text(5, -37, 'SNR 40 dB','Color',[0 0 0], 'FontSize', 18)
text(5, -47, 'SNR 50 dB','Color',[0 0 0], 'FontSize', 18)
axis([-50 20  -50 5])
print -depsc randomLevelMod.eps
%% Signal dependent deviation

nType = 4;
nRepetition = 3;
snrList = 10:10:50;
xRefSet = randn(lTestSig, nType);
averageResp = zeros(lTestSig, 1);
%nlAlpha = 0.0002*2.^(0:1/3:log2(0.2/0.0002));
alpp = 0.4;
inputLevel = -60:5:0;
clipLevel = 0;
sigDepLevel = zeros(length(inputLevel), length(snrList));
averagesqNoiseLevelRecord = zeros(length(inputLevel), length(snrList));
averagesqSignalLevelRecord = zeros(length(inputLevel), length(snrList));
for nn = 1:length(snrList)
    snr = snrList(nn);
    for mm = 1:length(inputLevel) 
        ampL = 10^(inputLevel(mm)/20);
        averageRespType = averageResp*0;
        sqNoiseLevelRecord = zeros(nType, 1);
        sqAverageLevelRecord = zeros(nType, 1);
        averageRespRecord = zeros(lTestSig, nType);
        averagesqNoiseLevel = 0;
        averagSignalLevelUnit = zeros(nType, 1);
        for kk = 1:nType
            xRefF = fft(xRefSet(:, kk));
            averageResp = averageResp*0;
            averagSignalLevel = 0;
            for jj = 1:nRepetition
                xRefFixF = xRefF;
                thetaL = mean(abs(xRefF))*10^(clipLevel/20);
                xRefFixF(abs(xRefFixF) < thetaL) = ...
                    thetaL * xRefFixF(abs(xRefFixF) < thetaL) ./ abs(xRefFixF(abs(xRefFixF) < thetaL));
                xRefFix = real(ifft(xRefFixF));
                xMesS = (exp(alpp*xRefFix*ampL) - 1)/alpp + ampL*randn(lTestSig, 1)/10^(snr/20);
                xMesSF = fft(xMesS);
                HrawS = xMesSF ./ xRefFixF;
                HrawSRslt(:, jj) = HrawS;
                averageResp = averageResp + HrawS;
            end
            averageResp = averageResp/nRepetition;
            averageRespType = averageRespType + averageResp;
            averageRespRecord(:, kk) = averageResp;
            sqNoise = 0;
            for jj = 1:nRepetition
                sqNoise = sqNoise + mean(abs(HrawSRslt(:, jj)-averageResp).^2);
                averagSignalLevel = averagSignalLevel + mean(abs(averageResp) .^2);
            end
            averagSignalLevel = averagSignalLevel/nRepetition;
            sqNoiseLevel = sqNoise/(nRepetition -1);
            sqNoiseLevelRecord(kk) = sqNoiseLevel;
            sqAverageLevelRecord(kk) = averagSignalLevel;
            averagesqNoiseLevel = averagesqNoiseLevel + sqNoiseLevel;
        end
        averageRespType = averageRespType/nType;
        averagesqNoiseLevel = averagesqNoiseLevel/nType;
        averagesqNoiseLevelRecord(mm, nn) = mean(sqNoiseLevelRecord);%averagesqNoiseLevel/nType;
        averagesqSignalLevelRecord(mm, nn) = mean(sqAverageLevelRecord);
        sqNoiseType = 0;
        for kk = 1:nType
            sqNoiseType = sqNoiseType + mean(abs(averageRespRecord(:, kk)-averageRespType).^2);
        end
        sqNoiseLevelType = sqNoiseType/(nType -1);
        sigDepLevel(mm,nn) = sqNoiseLevelType;
    end
end
%%
figure;
set(gcf,'position', [3501         642         700         300])
plot(inputLevel,10*log10(sigDepLevel), 'LineWidth', 2);grid on;
hold all;
semilogx(inputLevel, 10*log10(averagesqNoiseLevelRecord),'k--');
set(gca, 'FontSize', 15, 'LineWidth', 2)
xlabel('input level (dB)')
ylabel('output level (dB)')
text(-47, -61, '10','Color',[0 0 0], 'FontSize', 16)
text(-47, -71, '20','Color',[0 0 0], 'FontSize', 16)
text(-47, -81, '30','Color',[0 0 0], 'FontSize', 16)
text(-47, -91, '40','Color',[0 0 0], 'FontSize', 16)
text(-47, -101, '50','Color',[0 0 0], 'FontSize', 16)
text(-50, -45, 'SNR (dB)','Color',[0 0 0], 'FontSize', 16)
print -depsc nlOnDeviation.eps
%%
figure;
outLevel = 10*log10(averagesqSignalLevelRecord);
set(gcf,'position', [3501         642         700         300])
plot(inputLevel,10*log10(sigDepLevel)-outLevel, 'LineWidth', 3);grid on;
hold all;
semilogx(inputLevel, 10*log10(averagesqNoiseLevelRecord)-outLevel,'k--');
set(gca, 'FontSize', 15, 'LineWidth', 2)
xlabel('input level (dB)')
ylabel('relative level (dB rel. outout)')
%%

text(-57, -12, '10','Color',[0 0 0], 'FontSize', 16)
text(-57, -22, '20','Color',[0 0 0], 'FontSize', 16)
text(-57, -32, '30','Color',[0 0 0], 'FontSize', 16)
text(-57, -42, '40','Color',[0 0 0], 'FontSize', 16)
text(-57, -52, '50','Color',[0 0 0], 'FontSize', 16)
text(-59, -5, 'SNR (dB) for Sig.Dep.','Color',[0 0 0], 'FontSize', 16)
%print -depsc nlOnDeviation.eps
text(-5, -6, '10','Color',[0 0 0], 'FontSize', 16)
text(-5, -16, '20','Color',[0 0 0], 'FontSize', 16)
text(-5, -26, '30','Color',[0 0 0], 'FontSize', 16)
text(-5, -36, '40','Color',[0 0 0], 'FontSize', 16)
text(-5, -46, '50','Color',[0 0 0], 'FontSize', 16)
text(-16, -56, 'SNR (dB) for rand.','Color',[0 0 0], 'FontSize', 16)
%%
print -depsc nlOnDeviation.eps