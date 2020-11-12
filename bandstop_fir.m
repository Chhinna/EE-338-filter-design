F_samp = 260e3;

%Band Edge speifications
F_pl = 53.2e3;
F_sl = 57.2e3;
F_sh = 77.2e3;
F_ph = 81.2e3;

w_pl = F_pl*2*pi/F_samp;
w_sh = F_sh*2*pi/F_samp;
w_sl = F_sl*2*pi/F_samp;
w_ph = F_ph*2*pi/F_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end


N_min = ceil((A-7.95) / (2.285*0.0307*pi));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 36;

%Ideal bandstop impulse response of length "n"
w_c1 = (w_pl + w_sl)/2;
w_c2 = (w_ph + w_sh)/2;
bs_ideal =  ideal_lp(pi,n) -ideal_lp(w_c2,n) + ideal_lp(w_c1,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, F_samp);
plot(f,abs(H))
xlim([20000,130000])

yline(0.15,'-r','Magnitude = 0.15')
yline(0.85,'-r','Magnitude = 0.85')
yline(1.15,'-r','Magnitude = 1.15')
xline(53.2e3,'-m','f = 52.9kHz')
xline(57.2e3,'-m','f = 56.9kHz')
xline(77.2e3,'-m','f = 76.9kHz')
xline(81.2e3,'-m','f = 80.9kHz')

xlim([20000,120000])
ylim([0,1.5])
xlabel('Frequency (in Hz)')
ylabel('Magnitude')
title('Magnitude Plot')
hold on
grid on
