function output = sixTermFBdesign(fs, fl, fh, nInOct, mag)
startTic = tic;
fcList = fl * 2.^(0:1/nInOct:log2(fh/fl));
nChannels = length(fcList);
fbank = struct;
for ii = 1:nChannels
    halfL = round(mag*fs/fcList(ii));
    ttLocal = (-halfL:halfL)'/fs;
    fbank.filter(ii).w = sixtermCos(halfL*2+1) .* exp(1i*2*pi*fcList(ii)*ttLocal);
    fbank.filter(ii).halfL = halfL;
    fbank.filter(ii).tAxis = ttLocal;
end
output.fcList = fcList;
output.fs = fs;
output.fbank = fbank;
output.nChannels = nChannels;
output.elapsedTime = toc(startTic);
end