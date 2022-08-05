%% Test script for smoothed transfer function

% Licence
% Copyright 2022 Hideki Kawahara
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Please execute this using cell mode

%% clock speed check
[xm,fsm] = audioread("/Volumes/SSPL-UT/acousticTest/loopMusic1/mesSig20220802T205507.wav");
[xr,fsr] = audioread("/Volumes/SSPL-UT/acousticTest/loopMusic1/testSig20220802T205507.wav");
signalInfo = audioinfo("/Volumes/SSPL-UT/acousticTest/loopMusic1/testSig20220802T205507.wav");
 %%
%[xm,fsm] = audioread("/Volumes/SSPL-UT/acousticTest/loopMusic1/mesSig20220805T010413.wav");
%[xr,fsr] = audioread("/Volumes/SSPL-UT/acousticTest/loopMusic1/testSig20220805T010413.wav");
%signalInfo = audioinfo("/Volumes/SSPL-UT/acousticTest/loopMusic1/testSig20220805T010413.wav");
%mesInfo = audioinfo("/Volumes/SSPL-UT/acousticTest/loopMusic1/mesSig20220805T010413.wav");

signalInfoCell = split(signalInfo.Comment);
unitL = str2double(signalInfoCell{2});
nTune = str2double(signalInfoCell{4});
nRepetition = str2double(signalInfoCell{3});
mesInfoCell = split(mesInfo.Comment);
calDB = str2double(mesInfoCell{2});
%%
revSeg = xm(fsr + (1:fsr));
revSeg = revSeg(end:-1:1);
tic
pos1 = fftfilt(revSeg,fftfilt(hanning(11),xm(fsr+(1:unitL))));
pos2 = fftfilt(revSeg,fftfilt(hanning(11),xm(fsr+unitL+(1:unitL))));
pos3 = fftfilt(revSeg,fftfilt(hanning(11),xm(fsr+2*unitL+(1:unitL))));
pos4 = fftfilt(revSeg,fftfilt(hanning(11),xm(fsr+3*unitL+(1:unitL))));
toc
%%
w1 = hanning(441);
w1 = w1/sum(w1);
w2 = hanning(21);
w2 = w2/sum(w2);
w = w1;
w(221+(-10:10)) = w(221+(-10:10))-w2;
w = -w;
%%
positions = zeros(nTune*nRepetition,4);
peakLocId = 0;
avSlope = zeros(nTune,1);
rawSlope = zeros(nTune*(nRepetition-1),1);
for iTuneId = 1:nTune
    revSeg = xm(3*fsr + (1:fsr) + (iTuneId-1)*4*unitL);
    revSeg = revSeg(end:-1:1);
    tunePos = zeros(nRepetition,1);
    for jj = 1:nRepetition
        xmSeg = xm(fsr+(1:unitL) + (iTuneId-1)*4*unitL + (jj-1)*unitL);
        tmpSig = fftfilt(revSeg,fftfilt(w,xmSeg));
        [~, maxIdx] = max(tmpSig);
        peakLocId = peakLocId + 1;
        positions(peakLocId,1) = maxIdx;
        positions(peakLocId,2:4) = [tmpSig(maxIdx-1) tmpSig(maxIdx) tmpSig(maxIdx+1)];
        r = [tmpSig(maxIdx-1) tmpSig(maxIdx) tmpSig(maxIdx+1)] - tmpSig(maxIdx);
        maxp = -(r(3)-r(1))/(r(1)+r(3))/2;
        positions(peakLocId,1) = maxIdx + maxp;
        tunePos(jj) = maxIdx + maxp;
    end
    avSlope(iTuneId) = mean(diff(tunePos))/unitL;
    rawSlope((iTuneId-1)*(nRepetition-1)+(1:nRepetition-1)) = ...
        diff(tunePos)/unitL;
end
%%
tfrm = (3+(0:15)*unitL/fsr);
figure;stem(tfrm, (positions(:,1)-positions(1,1))/fsr*1000,"LineWidth",2);
grid on;
set(gca,"FontSize",15,"LineWidth",2)
xlabel("frame location (s)")
ylabel("time axis deviation (ms)")
%%
calRaw = 10.0 ^(calDB/20);
segSilent = xm(1:45000);
rmsSilentSPL = 20*log10(std(segSilent)*calRaw);
segSignal = xm(10*fsr:80*fsr);
rmsSignalSPL = 20*log10(std(segSignal)*calRaw);
segTest = xr(10*fsr:80*fsr,1);
rmsTestSPL = 20*log10(std(segTest)*calRaw);
aFilter = weightingFilter("A-weighting",fsr);
rmsSilentAwSPL = 20*log10(std(aFilter(segSilent))*calRaw);
rmsSignalAwSPL = 20*log10(std(aFilter(segSignal))*calRaw);
%%
pwSilent = abs(fft(segSilent.*blackman(length(segSilent)))).^2;
pwSignal = abs(fft(segSignal.*blackman(length(segSignal)))).^2;
pwTest = abs(fft(segTest.*blackman(length(segTest)))).^2;
pwSilent = pwSilent/sum(pwSilent)*std(segSilent*calRaw)^2;
pwSignal = pwSignal/sum(pwSignal)*std(segSignal*calRaw)^2;
pwTest = pwTest/sum(pwTest)*std(segSignal*calRaw)^2;
fxSilent = (0:length(segSilent)-1)'/length(segSilent)*fsr;
fxSignal = (0:length(segSignal)-1)'/length(segSignal)*fsr;
cumPwSilent = cumsum(pwSilent);
cumPwSignal = cumsum(pwSignal);
cumPwTest = cumsum(pwTest);
fxSignalH = fxSignal*2^(1/6);
fxSignalL = fxSignal*2^(-1/6);
fxSilentH = fxSilent*2^(1/6);
fxSilentL = fxSilent*2^(-1/6);
smSilent = interp1(fxSilent,cumPwSilent,fxSignalH,"linear","extrap") ...
    - interp1(fxSilent,cumPwSilent,fxSignalL,"linear","extrap");
figure;
semilogx(fxSignal, 10*log10(smSilent),"LineWidth",2);grid on;
hold all
%%
smSignal = interp1(fxSignal,cumPwSignal,fxSignalH,"linear","extrap") ...
    - interp1(fxSignal,cumPwSignal,fxSignalL,"linear","extrap");
semilogx(fxSignal, 10*log10(smSignal),"LineWidth",2);grid on;
smTest = interp1(fxSignal,cumPwTest,fxSignalH,"linear","extrap") ...
    - interp1(fxSignal,cumPwTest,fxSignalL,"linear","extrap");
semilogx(fxSignal, 10*log10(smTest),"LineWidth",2);grid on;
set(gca,"FontSize",15,"LineWidth",2)
axis([20 fsr/2 0 70])
xlabel("frequency (Hz)")
ylabel("sound pressure level (dB rel. 20 \mu Pa)")
legend("backgrund","measured","test","Location","south","Orientation","horizontal")
%%
figure;
semilogx(fxSignal,10*log10(smSignal./smTest),"LineWidth",2);grid on;
hold all
semilogx(fxSignal,10*log10(smSilent./smTest),"LineWidth",2);grid on;
axis([20 fsr/2 -40 20])
xlabel("frequency (Hz)")
ylabel("relative gain (dB)")
legend("LTI response","noise response","Location","south","Orientation","horizontal")
%%
