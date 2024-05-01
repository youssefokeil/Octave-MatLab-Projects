
# clear command window
clc;
# remove all variables from workspace
clear
# close all figures
close all

%%%%%% Frequencies initialization %%%%%%
# carrier frequency initialization
fc=1e4; #10KHz
# message frequency initialization
fm=100; #100Hz

%%%%%% Time and frequency definition %%%%%%
# Sampling time is identified using fc so that we can sample carrier signal tooo
ts= 0.1/fc;
# T is the last value of t on graph, we're drawing 2 periods of the message signal
T=2/fm;
# Number of time samples
N=ceil(T/ts);
# time vector
t= 0:ts:((N-1)*ts)

%%%%%% Frequency domain %%%%%%
# df=fs/N
df=1/T;
# sampling frequency, inverse of sampling tinme
fs=1/ts;
# frequency vector
if(rem(N,2)==0) #even
  f= -(0.5*fs): df : (0.5*fs-df)
else #odd
  f= -(0.5*fs-0.5*df): df : (0.5*fs-0.5*df)
end

%%%%%% Baseband signal definition in Time %%%%%%
# m(t) definition
m= sin(2*pi*fm*t);

%%%%%% frequency components of baseband sigmal %%%%%%
# fast fourier transform of baseband signal
M=fftshift(fft(m)/N);

%%%%%% Carrier signal definition in time %%%%%%
# c(t)
c= cos(2*pi*fc*t);

%%%%%% Carrier signal definition in frequency %%%%%%
# C(f) using fft
C=fftshift(fft(c)/N);

%%%%%% DSB-SC Modulated signal definition in time %%%%%%
# s(t)
s= m.*c;

%%%%%% DSB-SC Modulated signal definition in frequency %%%%%%
# S(f) using fft
S=fftshift(fft(s)/N);

%%%%%% Coherent detector (receiver) time response %%%%%%
# r(t)
r=s.*c;
# plotting receiver signal with respect to time
figure(1);
plot(t,r);
xlabel("Time (sec)");
ylabel("r(t)")

%%%%%% Coherent detector (receiver) frequency response %%%%%%
# R(f)
R=fftshift(fft(r)/N);
# plotting receiver signal with respect to time
figure(2);
plot(f,abs(R));
xlabel("Frequency (Hz)");
ylabel("R(f)")

%%%%%% Filter (ideal LPF) frequency response %%%%%%
# H(f)
H=abs(f)<1.1*fm;
# plotting filter signal with respect to frequency
figure(3);
plot(f,H);
xlabel("Frequency (Hz)");
ylabel("H(f)")
box off;

%%%%%% demodulated signal after ideal LPF %%%%%%
# M_rec(f)
M_rec=H.*R;
# m_rec(t)
m_rec=ifft(ifftshift(M_rec)*N)
# plotting received signal after filter with respect to time
figure(4);
plot(t,m_rec,"-r");
# adding original song compared to received one
hold on;
plot(t,m,"--b");
xlabel("Time (sec)");
ylabel("m(t)")
box off;
