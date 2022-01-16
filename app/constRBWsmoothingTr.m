function smoothedPw = constRBWsmoothingTr(pw, fs, nInOctave)
fxL = length(pw);
fx = (0:fxL-1)'/fxL*fs;
fu = fx * 2^(1/nInOctave/2);
fl = fx * 2^(-1/nInOctave/2);
fw = (fu-fl);
fw(1) = 1;
cumPw = cumsum(fx(2)*pw);
smoothedPw = (interp1(fx(:), cumPw, fu(:), 'linear','extrap') ...
    - interp1(fx(:), cumPw, fl(:), 'linear','extrap')) ./ fw(:);
cumSmPw = cumsum(fx(2)*smoothedPw);
smoothedPw = (interp1(fx(:), cumSmPw, fu(:), 'linear','extrap') ...
    - interp1(fx(:), cumSmPw, fl(:), 'linear','extrap')) ./ fw(:);
end