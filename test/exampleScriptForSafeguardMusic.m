%% Example script for maki safeguarded music

loopUnit = {'holiday.wav', 'laboratory.wav', 'legend.wav', 'space.wav'};
nRepetition = 4;
%--------
theFormAvdB = 0;
tmp = makeSafeguardedMonoMusic(loopUnit,nRepetition,theFormAvdB);
fnameLch = 'safeguarded0dBL.wav';
fnameRch = 'safeguarded0dBR.wav';
fnameMono = 'safeguarded0dBMono.wav';
y = tmp.safeGuardedSignal;
y = 0.8 * y /max(abs(y));
fs = tmp.samplingFrequency;
audiowrite(fnameLch, [y y*0], fs, "BitsPerSample", 24);
audiowrite(fnameRch, [y*0 y], fs, "BitsPerSample", 24);
audiowrite(fnameMono, y, fs, "BitsPerSample", 24);

%--------
theFormAvdB = -10;
tmp = makeSafeguardedMonoMusic(loopUnit,nRepetition,theFormAvdB);
fnameLch = 'safeguardedM10dBL.wav';
fnameRch = 'safeguardedM10dBR.wav';
fnameMono = 'safeguardedM10dBMono.wav';
y = tmp.safeGuardedSignal;
y = 0.8 * y /max(abs(y));
fs = tmp.samplingFrequency;
audiowrite(fnameLch, [y y*0], fs, "BitsPerSample", 24);
audiowrite(fnameRch, [y*0 y], fs, "BitsPerSample", 24);
audiowrite(fnameMono, y, fs, "BitsPerSample", 24);
