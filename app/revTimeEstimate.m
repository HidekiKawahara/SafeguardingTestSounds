function output = revTimeEstimate(decayMapOutput)
%% Test script
tx = decayMapOutput.txBuffer;
dBResp = 10*log10(decayMapOutput.spRangePowerResp);
interCept = mean(dBResp(abs(tx - 0.01)<0.001));
rtList = 0.05:0.01:1;
%respMedian = 10*log10(decayMapOutput.pwrDecayMapMedian);
respMedian = median(dBResp);
effAv = zeros(length(rtList),1);
for ii = 1:length(rtList)
    y = interCept - 60/rtList(ii)*tx(tx>0.01);
    err = abs(y-dBResp(tx>0.01)).^2;
    effAv(ii) = mean(err(y>respMedian+6));
end
[~, minIdx] = min(effAv);
y = interCept - 60/rtList(minIdx)*tx(tx>0.01);
idxOrg = 1:length(tx);
txSel = tx(tx>0.01);
idxSel = idxOrg(tx>0.01);
txSel = txSel(y>respMedian+6);
idxSel = idxSel(y>respMedian+6);
if max(diff(idxSel)) > 1
    idxEnd = 1;
    for ii = 1:length(idxSel)-1
        if idxSel(ii+1)-idxSel(ii) >2
            break
        end
        idxEnd = idxEnd + 1;
    end
    txSel = txSel(1:idxEnd);
    idxSel = idxSel(1:idxEnd);
end
%dBRespSel = dBResp(tx>0.01);
dBRespSel = dBResp(idxSel);
H = [ones(length(dBRespSel),1) txSel];
r = (H'*H)\(H'*dBRespSel);
t60 = -60/r(2);
output.t60 = t60;
output.interCept = r(1);
output.slope = r(2);
output.dBRespSel = dBRespSel;
output.txSel = txSel;
