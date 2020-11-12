%Butterworth Analog LPF parameters
Wc = 1.08;              %cut-off frequency
N = 7;                  %order 

%poles of Butterworth polynomial of degree 7 in the open CLHP 
p1 = -0.0546 + 1i*0.9871;
p2 = -0.0564 - 1i*0.9871;
p3 = -0.2037 + 1i*0.2645;
p4 = -0.2037 - 1i*0.2645;
p5 = -0.1491 + 1i*0.7226;
p6 = -0.1491 - 1i*0.7226;


%Band Edge speifications
F_pl = 53.2e3;
F_sl = 57.2e3;
F_sh = 77.2e3;
F_ph = 81.2e3;

%Transformed Band Edge specs using Bilinear Transformation
F_samp = 260e3;         
W_pl = tan(F_pl/F_samp*pi);
W_sl = tan(F_sl/F_samp*pi); 
W_sh = tan(F_sh/F_samp*pi);
W_ph = tan(F_ph/F_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(W_pl*W_ph);
B = W_ph-W_pl;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 ],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;
disp(ds);
disp(ns);
%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz);                                           %frequency response
disp(dz);
disp(nz);
%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, F_samp);
plot(f,abs(H))

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