%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 6;

% Open CLHP Poles of the Chebyshev Polynomial of order 4
p1 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p2 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p3 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p4 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);        
p5 = -sin(5*pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(5*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p6 = -sin(5*pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(5*pi/(2*N))*cosh(asinh(1/epsilon)/N);       

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
n3 = [1 -p5-p6 p5*p6];
den = conv(n1,n2);             
den = conv(den,n3);
num = [den(7)*sqrt(1/(1+epsilon*epsilon))]; 
%Band Edge speifications
F_sl = 41.7e3;
F_pl = 45.7e3;
F_ph = 65.7e3;
F_sh = 69.7e3;

%Transformed Band Edge specs using Bilinear Transformation
F_samp = 330e3;
Ws1 = tan(F_sl/F_samp*pi);          
Wp1 = tan(F_pl/F_samp*pi);
Wp2 = tan(F_ph/F_samp*pi);
Ws2 = tan(F_sh/F_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(Wp1*Wp2);
B = Wp2-Wp1;

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BPF
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       %frequency response in dB

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, F_samp);
plot(f,abs(H))

yline(0.15,'-r','Magnitude = 0.15');
yline(0.85,'-r','Magnitude = 0.85');
yline(1,'-r','Magnitude = 1');
xline(41.7e3,'-m','f = 41.7kHz');
xline(45.7e3,'-m','f = 45.7kHz');
xline(65.7e3,'-m','f = 65.7kHz');
xline(69.7e3,'-m','f = 69.7kHz');

xlim([20000,130000]);
ylim([0,1.5]);
xlabel('Frequency (in Hz)');
ylabel('Magnitude');
title('Magnitude Plot');
hold on
grid on