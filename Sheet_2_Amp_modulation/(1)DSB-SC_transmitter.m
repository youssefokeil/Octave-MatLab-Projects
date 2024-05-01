# Created by Octave 9.1.0, Tue Apr 30 14:26:57 2024 GMT <unknown@Okeil>

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
# plotting m(t) with t
figure(1);
plot(t,m);
xlabel("Time (sec)");
ylabel("m(t)");

%%%%%% frequency components of baseband sigmal %%%%%%
# fast fourier transform of baseband signal
M=fftshift(fft(m)/N);
# plotting frequency component of signal
figure(2);
plot(f, abs(M));
xlabel("Frequency (Hz)");
ylabel("M(f)");

%%%%%% Carrier signal definition in time %%%%%%
# c(t)
c= cos(2*pi*fc*t);
# plotting carrier signal with respect to time
figure(3);
plot(t,c);
xlabel("Time (sec)");
ylabel("c(t)");

%%%%%% Carrier signal definition in frequency %%%%%%
# C(f) using fft
C=fftshift(fft(c)/N);
# plotting carrier signal with respect to frequency
figure(4);
plot(f,abs(C));
xlabel("Frequency (Hz)");
ylabel("C(f)")

%%%%%% DSB-SC Modulated signal definition in time %%%%%%
# s(t)
s= m.*c;
# plotting modulated signal with respect to time
figure(5);
plot(t,s);
xlabel("Time (sec)");
ylabel("s(t)")

%%%%%% DSB-SC Modulated signal definition in frequency %%%%%%
# S(f) using fft
S=fftshift(fft(s)/N);
# plotting modulated signal with respect to frequency
figure(6);
plot(f,abs(S));
xlabel("Frequency (Hz)");
ylabel("S(f)")


# change fc, fm from above to 500 & 10 respectively and see the figures
