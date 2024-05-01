# clear command window
clc;
# remove all variables from workspace
clear
# close all figures
close all
# importing signal package
pkg load signal

%%%%%% Frequencies initialization %%%%%%
# carrier frequency initialization
fc=1e4; #10KHz
# message frequency initialization
fm=100; #100Hz
Am=10;
Ac=1;
kf=100;
beta=kf*Am/fm;

%%%%%% Time and frequency definition %%%%%%
# Sampling time is identified using fc so that we can sample carrier signal tooo
ts= 0.1/fc/10;
# T is the last value of t on graph, we're drawing 2 periods of the message signal
T=10/fm;
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
  f= -(0.5*fs): df : (0.5*fs-df);
else #odd
  f= -(0.5*fs-0.5*df): df : (0.5*fs-0.5*df);
end

%%%%%%% Baseband Signal %%%%%%%%
m= Am*sin(2*pi*fm*t);
# fast fourier transform of baseband signal
M=fftshift(fft(m)/N);

%%%%%%%% Modulated Signal %%%%%%%%%%%%%
# time domain, cumtrapz is integration function
s= Ac*cos(2*pi*fc*t+2*pi*kf*cumtrapz(t,m));
# frequency domain
S=fftshift(fft(s)/N);

%%%%%%%%%%% plotting FM modulated signal %%%%%%%%
# plotting modulated signal with respect to time
figure(1);
plot(t,s);
xlabel("Time (sec)");
ylabel("s(t)")

# plotting modulated signal with respect to frequency
figure(2);
# using stem instead of plot
stem(f,abs(S));
xlabel("Frequency (Hz)");
ylabel("|S(f)|")
box off;


%%%%% Receiver %%%%%%%%%%%%%%
%%%% starting with differentiator %%%%
g_rec=diff(s)./diff(t);
t_new=t(2:end);
%%%%%% Envelope detector %%%%%%%%%
g_rec=abs(hilbert(g_rec));
%%%%%%% Removing DC level (mean of signal) %%%%%%%%%
g_rec=g_rec-mean(g_rec);
%%%%%%% Changing time vector to remove transient response %%%%%%%%%%%%
index=t>1/fm & t<9/fm;
t_new=t_new(index);
g_rec=g_rec(index);
%%%% Plotting original signal %%%%%%%
figure(3)
plot(t,m/max(m),"-b");
hold on;
%%%%%% Plotting retrieved signal %%%%%%%%%%
plot(t_new,g_rec/max(g_rec),"--r");
hold on;
legend("Original Signal", "Retrieved Signal");
